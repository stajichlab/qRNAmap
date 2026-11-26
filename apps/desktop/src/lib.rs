use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::Command;

use serde::{Deserialize, Serialize};
use tauri::{WebviewUrl, WebviewWindowBuilder};

use qrnamap_core::{
    gene_counts_from_abundance, kallisto_quant, parse_tx2gene, zero_warning,
    Quantifier, RunConfig, Sample, Strandedness,
};

/// Information about native execution capabilities on this machine.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NativeCapabilities {
    pub available: bool,
    pub kallisto_path: Option<String>,
    pub kallisto_version: Option<String>,
    pub num_cpus: usize,
    pub temp_dir: String,
}

/// Request to prepare or verify a native index on disk.
#[derive(Debug, Clone, Deserialize)]
pub struct PrepareIndexRequest {
    pub index_path: Option<String>,
    pub index_bytes: Option<Vec<u8>>,
    pub index_gz_bytes: Option<Vec<u8>>,
    pub fasta_path: Option<String>,
    pub fasta_bytes: Option<Vec<u8>>,
}

/// Response with path and size of prepared native index.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrepareIndexResponse {
    pub index_path: String,
    pub index_size: u64,
}

/// Request to run native quantification on one sample.
#[derive(Debug, Clone, Deserialize)]
pub struct NativeQuantRequest {
    pub name: String,
    pub fastq1: Vec<String>,
    pub fastq2: Vec<String>,
    pub index_path: String,
    pub tx2gene: String,
    pub strandedness: Option<String>,
    pub threads: Option<usize>,
    pub fragment_mean: Option<f64>,
    pub fragment_sd: Option<f64>,
    pub out_dir: Option<String>,
}

/// Response returned from native quantification to frontend.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NativeQuantResponse {
    pub sample: String,
    pub n_processed: u64,
    pub n_pseudoaligned: u64,
    pub p_pseudoaligned: f64,
    pub whole_file: bool,
    pub gene_counts: BTreeMap<String, u64>,
    pub abundance_tsv: String,
    pub warning: String,
    pub args: Vec<String>,
}

/// Request to parse and resolve a local sample sheet from disk.
#[derive(Debug, Clone, Deserialize)]
pub struct ReadSheetRequest {
    pub path: String,
    pub base_dir: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResolvedSampleRow {
    pub sample: String,
    pub name: String,
    pub fastq1: Vec<String>,
    pub fastq2: Vec<String>,
    pub fastq1_raw: Vec<String>,
    pub fastq2_raw: Vec<String>,
    pub fastq1_exists: Vec<bool>,
    pub fastq2_exists: Vec<bool>,
    pub covariates: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadSheetResponse {
    pub path: String,
    pub dir: String,
    pub content: String,
    pub columns: Vec<(String, String)>,
    pub samples: Vec<ResolvedSampleRow>,
    pub all_files_exist: bool,
    pub warnings: Vec<String>,
}

/// Look for the `kallisto` binary on the system.
pub fn find_kallisto() -> Option<PathBuf> {
    if let Ok(custom) = std::env::var("QRNAMAP_KALLISTO") {
        let p = PathBuf::from(custom);
        if p.is_file() {
            return Some(p);
        }
    }

    if let Ok(path_var) = std::env::var("PATH") {
        for dir in std::env::split_paths(&path_var) {
            let candidate = dir.join("kallisto");
            if candidate.is_file() {
                return Some(candidate);
            }
        }
    }

    // Standard fallback locations on macOS and Linux
    let home = std::env::var("HOME").unwrap_or_default();
    let well_known = [
        "/opt/homebrew/bin/kallisto",
        "/usr/local/bin/kallisto",
        "/usr/bin/kallisto",
        "/bin/kallisto",
        "/opt/conda/bin/kallisto",
    ];

    for path in &well_known {
        let p = PathBuf::from(path);
        if p.is_file() {
            return Some(p);
        }
    }

    if !home.is_empty() {
        let user_paths = [
            format!("{home}/.cargo/bin/kallisto"),
            format!("{home}/.pixi/bin/kallisto"),
            format!("{home}/.local/bin/kallisto"),
            format!("{home}/miniforge3/bin/kallisto"),
            format!("{home}/miniconda3/bin/kallisto"),
            format!("{home}/anaconda3/bin/kallisto"),
            format!("{home}/projects/funannotate/.pixi/envs/default/bin/kallisto"),
        ];
        for path in &user_paths {
            let p = PathBuf::from(path);
            if p.is_file() {
                return Some(p);
            }
        }
    }

    None
}

/// Get kallisto version text.
pub fn get_kallisto_version(bin: &Path) -> Result<String, String> {
    let output = Command::new(bin)
        .arg("version")
        .output()
        .map_err(|e| format!("cannot run {}: {e}", bin.display()))?;

    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
    if !stdout.is_empty() {
        Ok(stdout)
    } else if !stderr.is_empty() {
        Ok(stderr)
    } else {
        Ok("kallisto (unknown version)".to_string())
    }
}

/// Expand `~` and relative paths.
pub fn resolve_path_str(path: &str, base_dir: Option<&Path>) -> PathBuf {
    let p = if path.starts_with("~/") || path == "~" {
        if let Ok(home) = std::env::var("HOME") {
            PathBuf::from(home).join(&path[2..])
        } else {
            PathBuf::from(path)
        }
    } else {
        PathBuf::from(path)
    };

    if p.is_absolute() {
        p
    } else if let Some(base) = base_dir {
        base.join(p)
    } else {
        std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")).join(p)
    }
}

// ---------------------------------------------------------------------------------------------
// Commands module to isolate Tauri macros and avoid namespace collisions in tests

pub mod commands {
    use super::*;

    #[tauri::command]
    pub fn open_url(url: String) -> Result<(), String> {
        if url.starts_with("http://") || url.starts_with("https://") {
            open::that(&url).map_err(|e| e.to_string())
        } else {
            Err("disallowed scheme".to_string())
        }
    }

    #[tauri::command]
    pub fn check_native_capabilities() -> NativeCapabilities {
        let kallisto = find_kallisto();
        let (kallisto_path, kallisto_version, available) = match kallisto {
            Some(ref p) => {
                let ver = get_kallisto_version(p).ok();
                (Some(p.to_string_lossy().to_string()), ver, true)
            }
            None => (None, None, false),
        };

        let num_cpus = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let temp_dir = std::env::temp_dir().to_string_lossy().to_string();

        NativeCapabilities {
            available,
            kallisto_path,
            kallisto_version,
            num_cpus,
            temp_dir,
        }
    }

    #[tauri::command]
    pub fn prepare_native_index(req: PrepareIndexRequest) -> Result<PrepareIndexResponse, String> {
        let scratch_dir = std::env::temp_dir().join("qrnamap_desktop");
        fs::create_dir_all(&scratch_dir).map_err(|e| format!("cannot create scratch dir: {e}"))?;

        if let Some(ref path_str) = req.index_path {
            let p = resolve_path_str(path_str, None);
            if p.is_file() {
                let size = fs::metadata(&p).map(|m| m.len()).unwrap_or(0);
                return Ok(PrepareIndexResponse {
                    index_path: p.to_string_lossy().to_string(),
                    index_size: size,
                });
            }
            return Err(format!("specified index path does not exist: {}", p.display()));
        }

        let target_idx = scratch_dir.join("transcripts.idx");

        if let Some(ref gz) = req.index_gz_bytes {
            let mut decoder = flate2::read::GzDecoder::new(&gz[..]);
            let mut out = File::create(&target_idx)
                .map_err(|e| format!("cannot create {}: {e}", target_idx.display()))?;
            std::io::copy(&mut decoder, &mut out)
                .map_err(|e| format!("cannot decompress index.gz: {e}"))?;
            let size = fs::metadata(&target_idx).map(|m| m.len()).unwrap_or(0);
            return Ok(PrepareIndexResponse {
                index_path: target_idx.to_string_lossy().to_string(),
                index_size: size,
            });
        }

        if let Some(ref bytes) = req.index_bytes {
            fs::write(&target_idx, bytes)
                .map_err(|e| format!("cannot write {}: {e}", target_idx.display()))?;
            let size = bytes.len() as u64;
            return Ok(PrepareIndexResponse {
                index_path: target_idx.to_string_lossy().to_string(),
                index_size: size,
            });
        }

        if req.fasta_path.is_some() || req.fasta_bytes.is_some() {
            let kallisto = find_kallisto()
                .ok_or_else(|| "kallisto binary not found to build index".to_string())?;

            let fasta_file = if let Some(ref path_str) = req.fasta_path {
                resolve_path_str(path_str, None)
            } else {
                let p = scratch_dir.join("transcripts.fa");
                fs::write(&p, req.fasta_bytes.as_ref().unwrap())
                    .map_err(|e| format!("cannot write transcripts.fa: {e}"))?;
                p
            };

            let output = Command::new(&kallisto)
                .arg("index")
                .arg("-i")
                .arg(&target_idx)
                .arg(&fasta_file)
                .output()
                .map_err(|e| format!("failed to start kallisto index: {e}"))?;

            if !output.status.success() {
                let stderr = String::from_utf8_lossy(&output.stderr);
                return Err(format!("kallisto index failed: {stderr}"));
            }

            let size = fs::metadata(&target_idx).map(|m| m.len()).unwrap_or(0);
            return Ok(PrepareIndexResponse {
                index_path: target_idx.to_string_lossy().to_string(),
                index_size: size,
            });
        }

        Err("no index or FASTA provided".to_string())
    }

    #[tauri::command]
    pub fn run_native_quant(req: NativeQuantRequest) -> Result<NativeQuantResponse, String> {
        let kallisto = find_kallisto()
            .ok_or_else(|| "kallisto binary not found on this system. Install kallisto (e.g. brew install kallisto).".to_string())?;

        let index_path = resolve_path_str(&req.index_path, None);
        if !index_path.is_file() {
            return Err(format!(
                "kallisto index does not exist at {}",
                index_path.display()
            ));
        }

        // Verify all FASTQ files exist
        let mut resolved_r1 = Vec::new();
        for f in &req.fastq1 {
            let p = resolve_path_str(f, None);
            if !p.is_file() {
                return Err(format!("mate 1 FASTQ file not found: {}", p.display()));
            }
            resolved_r1.push(p.to_string_lossy().to_string());
        }

        let mut resolved_r2 = Vec::new();
        for f in &req.fastq2 {
            let p = resolve_path_str(f, None);
            if !p.is_file() {
                return Err(format!("mate 2 FASTQ file not found: {}", p.display()));
            }
            resolved_r2.push(p.to_string_lossy().to_string());
        }

        if resolved_r1.is_empty() {
            return Err(format!("sample '{}' has no mate 1 reads", req.name));
        }
        if !resolved_r2.is_empty() && resolved_r1.len() != resolved_r2.len() {
            return Err(format!(
                "sample '{}' mate 1 ({}) and mate 2 ({}) lane counts differ",
                req.name,
                resolved_r1.len(),
                resolved_r2.len()
            ));
        }

        let out_dir = match req.out_dir {
            Some(ref d) => resolve_path_str(d, None).join(&req.name),
            None => std::env::temp_dir()
                .join("qrnamap_quant")
                .join(&req.name),
        };
        fs::create_dir_all(&out_dir).map_err(|e| format!("cannot create out dir {}: {e}", out_dir.display()))?;

        let strandedness = match req.strandedness.as_deref().unwrap_or("unstranded") {
            "forward" => Strandedness::Forward,
            "reverse" => Strandedness::Reverse,
            _ => Strandedness::Unstranded,
        };

        let threads = req.threads.unwrap_or_else(|| {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
        }).max(1) as u32;

        let fragment_len = if resolved_r2.is_empty() {
            let mean = req.fragment_mean.ok_or_else(|| {
                format!("single-end sample '{}' requires fragment length mean", req.name)
            })?;
            let sd = req.fragment_sd.ok_or_else(|| {
                format!("single-end sample '{}' requires fragment length SD", req.name)
            })?;
            if !(mean > 0.0 && sd > 0.0 && mean.is_finite() && sd.is_finite()) {
                return Err(format!(
                    "sample '{}' fragment length mean and SD must be finite positive numbers",
                    req.name
                ));
            }
            Some((mean, sd))
        } else {
            None
        };

        let sample = Sample {
            name: req.name.clone(),
            fastq1: resolved_r1,
            fastq2: resolved_r2,
            covariates: BTreeMap::new(),
        };

        let mut cfg = RunConfig::default();
        cfg.strandedness = strandedness;
        cfg.threads = threads;
        cfg.fragment_len = fragment_len;
        cfg.quantifier = Quantifier::Kallisto;

        let call = kallisto_quant(
            &sample,
            &index_path.to_string_lossy(),
            &out_dir.to_string_lossy(),
            &cfg,
        ).map_err(|e| format!("config error for sample '{}': {e}", req.name))?;

        // Run kallisto directly
        let mut cmd = Command::new(&kallisto);
        cmd.args(&call.args);
        cmd.current_dir(&out_dir);

        let output = cmd.output().map_err(|e| format!("failed to start kallisto: {e}"))?;

        let stderr = String::from_utf8_lossy(&output.stderr);
        if !output.status.success() {
            if stderr.contains("zero reads pseudoaligned") || stderr.contains("no reads pseudoaligned") {
                return Err(format!(
                    "sample '{}': no reads pseudoaligned to transcripts (kallisto exited 1). Check strandedness setting or reference.",
                    req.name
                ));
            }
            let tail = stderr.lines().rev().take(5).collect::<Vec<_>>();
            let mut tail_rev = tail;
            tail_rev.reverse();
            return Err(format!(
                "kallisto failed for sample '{}' (status {}):\n{}",
                req.name,
                output.status,
                tail_rev.join("\n")
            ));
        }

        let abundance_path = out_dir.join("abundance.tsv");
        let run_info_path = out_dir.join("run_info.json");

        if !abundance_path.is_file() {
            return Err(format!(
                "kallisto completed successfully but did not produce {}",
                abundance_path.display()
            ));
        }

        let abundance_tsv = fs::read_to_string(&abundance_path)
            .map_err(|e| format!("cannot read {}: {e}", abundance_path.display()))?;

        let run_info_text = fs::read_to_string(&run_info_path)
            .map_err(|e| format!("cannot read {}: {e}", run_info_path.display()))?;

        let run_info_val: serde_json::Value = serde_json::from_str(&run_info_text)
            .map_err(|e| format!("cannot parse run_info.json: {e}"))?;

        let n_processed = run_info_val.get("n_processed").and_then(|v| v.as_u64()).unwrap_or(0);
        let n_pseudoaligned = run_info_val.get("n_pseudoaligned").and_then(|v| v.as_u64()).unwrap_or(0);
        let p_pseudoaligned = run_info_val
            .get("p_pseudoaligned")
            .and_then(|v| v.as_f64())
            .map(|p| p / 100.0)
            .unwrap_or(0.0);

        // Parse tx2gene
        let tx2gene_map = if req.tx2gene.starts_with('/') || req.tx2gene.ends_with(".tsv") || req.tx2gene.ends_with(".txt") {
            let p = resolve_path_str(&req.tx2gene, None);
            if p.is_file() {
                let tsv_text = fs::read_to_string(&p)
                    .map_err(|e| format!("cannot read tx2gene file {}: {e}", p.display()))?;
                parse_tx2gene(&tsv_text).map_err(|e| format!("cannot parse tx2gene table: {e}"))?
            } else {
                parse_tx2gene(&req.tx2gene).map_err(|e| format!("cannot parse tx2gene table: {e}"))?
            }
        } else {
            parse_tx2gene(&req.tx2gene).map_err(|e| format!("cannot parse tx2gene table: {e}"))?
        };

        let gene_counts = gene_counts_from_abundance(&abundance_tsv, &tx2gene_map)
            .map_err(|e| format!("cannot calculate gene counts from abundance: {e}"))?;

        let mut warnings = Vec::new();
        if let Some(w) = zero_warning(&req.name, &gene_counts) {
            warnings.push(w);
        }
        if p_pseudoaligned < 0.5 {
            warnings.push(format!(
                "pseudoaligned fraction is low ({:.1}%). Check strandedness setting (currently '{}').",
                p_pseudoaligned * 100.0,
                req.strandedness.as_deref().unwrap_or("unstranded")
            ));
        }

        Ok(NativeQuantResponse {
            sample: req.name,
            n_processed,
            n_pseudoaligned,
            p_pseudoaligned,
            whole_file: true,
            gene_counts,
            abundance_tsv,
            warning: warnings.join(" "),
            args: call.args,
        })
    }

    #[tauri::command]
    pub fn read_local_sample_sheet(req: ReadSheetRequest) -> Result<ReadSheetResponse, String> {
        let sheet_path = resolve_path_str(&req.path, None);
        if !sheet_path.is_file() {
            return Err(format!("sample sheet not found at {}", sheet_path.display()));
        }

        let content = fs::read_to_string(&sheet_path)
            .map_err(|e| format!("cannot read sample sheet: {e}"))?;

        let parsed = qrnamap_core::SampleSheet::parse(&content)
            .map_err(|errs| format!("invalid sample sheet: {}", errs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join("; ")))?;

        let sheet_dir = sheet_path.parent().unwrap_or(Path::new("."));
        let base_dir = req.base_dir.as_deref().map(|b| resolve_path_str(b, None)).unwrap_or_else(|| sheet_dir.to_path_buf());

        let mut samples = Vec::new();
        let mut all_files_exist = true;
        let mut warnings = Vec::new();

        for s in &parsed.samples {
            let mut r1_paths = Vec::new();
            let mut r1_exists = Vec::new();
            for f in &s.fastq1 {
                let res = resolve_path_str(f, Some(&base_dir));
                let ex = res.is_file();
                if !ex {
                    all_files_exist = false;
                    warnings.push(format!("sample '{}' mate 1 missing: {}", s.name, res.display()));
                }
                r1_paths.push(res.to_string_lossy().to_string());
                r1_exists.push(ex);
            }

            let mut r2_paths = Vec::new();
            let mut r2_exists = Vec::new();
            for f in &s.fastq2 {
                let res = resolve_path_str(f, Some(&base_dir));
                let ex = res.is_file();
                if !ex {
                    all_files_exist = false;
                    warnings.push(format!("sample '{}' mate 2 missing: {}", s.name, res.display()));
                }
                r2_paths.push(res.to_string_lossy().to_string());
                r2_exists.push(ex);
            }

            samples.push(ResolvedSampleRow {
                sample: s.name.clone(),
                name: s.name.clone(),
                fastq1: r1_paths,
                fastq2: r2_paths,
                fastq1_raw: s.fastq1.clone(),
                fastq2_raw: s.fastq2.clone(),
                fastq1_exists: r1_exists,
                fastq2_exists: r2_exists,
                covariates: s.covariates.clone(),
            });
        }

        let columns = parsed
            .columns
            .into_iter()
            .map(|(name, typ)| {
                let t_str = match typ {
                    qrnamap_core::ColType::Factor => "factor".to_string(),
                    qrnamap_core::ColType::Numeric => "numeric".to_string(),
                };
                (name, t_str)
            })
            .collect();

        Ok(ReadSheetResponse {
            path: sheet_path.to_string_lossy().to_string(),
            dir: sheet_dir.to_string_lossy().to_string(),
            content,
            columns,
            samples,
            all_files_exist,
            warnings,
        })
    }

    #[tauri::command]
    pub fn pick_file_dialog(title: Option<String>, filter_name: Option<String>, filter_exts: Option<Vec<String>>) -> Result<Option<String>, String> {
        let mut dialog = rfd::FileDialog::new();
        if let Some(t) = title {
            dialog = dialog.set_title(&t);
        }
        if let (Some(name), Some(exts)) = (filter_name, filter_exts) {
            let str_refs: Vec<&str> = exts.iter().map(|s| s.as_str()).collect();
            dialog = dialog.add_filter(&name, &str_refs);
        }
        let res = dialog.pick_file().map(|p| p.to_string_lossy().to_string());
        Ok(res)
    }

    #[tauri::command]
    pub fn pick_directory_dialog(title: Option<String>) -> Result<Option<String>, String> {
        let mut dialog = rfd::FileDialog::new();
        if let Some(t) = title {
            dialog = dialog.set_title(&t);
        }
        let res = dialog.pick_folder().map(|p| p.to_string_lossy().to_string());
        Ok(res)
    }

    #[tauri::command]
    pub fn resolve_path(path: String, base_dir: Option<String>) -> Result<serde_json::Value, String> {
        let base = base_dir.as_deref().map(|b| resolve_path_str(b, None));
        let resolved = resolve_path_str(&path, base.as_deref());
        let exists = resolved.is_file() || resolved.is_dir();
        let is_file = resolved.is_file();
        let size = if is_file {
            fs::metadata(&resolved).map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };

        Ok(serde_json::json!({
            "input": path,
            "resolved": resolved.to_string_lossy(),
            "exists": exists,
            "isFile": is_file,
            "size": size,
        }))
    }
}

pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            commands::open_url,
            commands::check_native_capabilities,
            commands::prepare_native_index,
            commands::run_native_quant,
            commands::read_local_sample_sheet,
            commands::pick_file_dialog,
            commands::pick_directory_dialog,
            commands::resolve_path,
        ])
        .setup(|app| {
            #[allow(unused_variables)]
            let window = WebviewWindowBuilder::new(app, "main", WebviewUrl::default())
                .title("qRNAmap")
                .inner_size(1100.0, 820.0)
                .min_inner_size(800.0, 600.0)
                .resizable(true)
                .on_navigation(|url| {
                    if url.scheme() == "http" || url.scheme() == "https" {
                        let _ = open::that(url.as_str());
                        return false;
                    }
                    true
                })
                .on_new_window(|url, _features| {
                    if url.scheme() == "http" || url.scheme() == "https" {
                        let _ = open::that(url.as_str());
                    }
                    tauri::webview::NewWindowResponse::Deny
                })
                .build()?;

            #[cfg(debug_assertions)]
            {
                let _ = window.open_devtools();
            }
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running qRNAmap desktop application");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_kallisto() {
        let k = find_kallisto();
        assert!(k.is_some(), "kallisto should be found on this machine");
        let path = k.unwrap();
        assert!(path.is_file(), "kallisto path must be an existing file");
        let version = get_kallisto_version(&path).expect("should get version");
        assert!(version.to_lowercase().contains("kallisto"), "version string should contain 'kallisto'");
    }

    #[test]
    fn test_check_capabilities() {
        let caps = commands::check_native_capabilities();
        assert!(caps.available);
        assert!(caps.kallisto_path.is_some());
        assert!(caps.num_cpus >= 1);
    }

    #[test]
    fn test_resolve_path_str() {
        let resolved = resolve_path_str("samples.txt", Some(Path::new("/tmp")));
        assert_eq!(resolved, PathBuf::from("/tmp/samples.txt"));

        let abs = resolve_path_str("/absolute/path", Some(Path::new("/tmp")));
        assert_eq!(abs, PathBuf::from("/absolute/path"));
    }

    #[test]
    fn test_read_local_sample_sheet() {
        let sheet_path = "/Users/jstajich/projects/Dhingra/samples.txt";
        if Path::new(sheet_path).is_file() {
            let res = commands::read_local_sample_sheet(ReadSheetRequest {
                path: sheet_path.to_string(),
                base_dir: None,
            }).expect("should parse samples.txt");

            assert!(!res.samples.is_empty(), "should have parsed samples");
            assert_eq!(res.samples[0].sample, "AF293_Azole_2");
            assert_eq!(res.samples[0].name, "AF293_Azole_2");
            assert!(res.samples[0].fastq1[0].contains("LCS9590_Af293_Azole_2_Clean_Data1.fq.gz"));
            assert!(res.samples[0].fastq1_exists[0], "first sample fastq1 file should exist on disk");
            assert!(res.samples[0].fastq2_exists[0], "first sample fastq2 file should exist on disk");
        }
    }

    #[test]
    fn test_native_quant_flow() {
        if find_kallisto().is_none() {
            return;
        }

        let dir = tempfile::tempdir().unwrap();
        let fa_path = dir.path().join("transcripts.fa");
        let tx1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let tx2_seq = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
        let fa_content = format!(">tx1\n{}\n>tx2\n{}\n", tx1_seq, tx2_seq);
        fs::write(&fa_path, &fa_content).unwrap();

        // 1. Prepare native index
        let prep = commands::prepare_native_index(PrepareIndexRequest {
            index_path: None,
            index_bytes: None,
            index_gz_bytes: None,
            fasta_path: Some(fa_path.to_string_lossy().to_string()),
            fasta_bytes: None,
        }).expect("should build native kallisto index");

        assert!(Path::new(&prep.index_path).is_file());
        assert!(prep.index_size > 0);

        // 2. Create FASTQ files
        let fq1 = dir.path().join("sample_R1.fastq");
        let fq2 = dir.path().join("sample_R2.fastq");

        let read1 = &tx1_seq[..50];
        let qual = "I".repeat(50);
        let fq1_content = format!("@read1/1\n{}\n+\n{}\n", read1, qual);
        let read2 = &tx1_seq[50..100];
        let fq2_content = format!("@read1/2\n{}\n+\n{}\n", read2, qual);

        fs::write(&fq1, &fq1_content).unwrap();
        fs::write(&fq2, &fq2_content).unwrap();

        let tx2gene = "tx1\tgene1\ntx2\tgene2\n";

        // 3. Run native quant
        let quant_res = commands::run_native_quant(NativeQuantRequest {
            name: "test_sample".to_string(),
            fastq1: vec![fq1.to_string_lossy().to_string()],
            fastq2: vec![fq2.to_string_lossy().to_string()],
            index_path: prep.index_path,
            tx2gene: tx2gene.to_string(),
            strandedness: Some("unstranded".to_string()),
            fragment_mean: None,
            fragment_sd: None,
            threads: Some(1),
            out_dir: Some(dir.path().join("quant_out").to_string_lossy().to_string()),
        }).expect("should successfully run native quant");

        assert_eq!(quant_res.sample, "test_sample");
        assert_eq!(quant_res.n_processed, 1);
        assert_eq!(quant_res.n_pseudoaligned, 1);
        assert_eq!(quant_res.whole_file, true);
        assert!(quant_res.gene_counts.contains_key("gene1"));
        assert_eq!(*quant_res.gene_counts.get("gene1").unwrap(), 1);
    }

    #[test]
    fn test_dhingra_real_data_quant() {
        if find_kallisto().is_none() {
            return;
        }

        let dhingra_transcripts = Path::new("/Users/jstajich/projects/Dhingra/results/transcripts.fa.zst");
        let dhingra_tx2gene = Path::new("/Users/jstajich/projects/Dhingra/results/tx2gene.tsv.zst");
        let r1 = Path::new("/Users/jstajich/projects/Dhingra/input/LCS9590_Af293_Azole_2_Clean_Data1.fq.gz");
        let r2 = Path::new("/Users/jstajich/projects/Dhingra/input/LCS9590_Af293_Azole_2_Clean_Data2.fq.gz");

        if !dhingra_transcripts.is_file() || !dhingra_tx2gene.is_file() || !r1.is_file() || !r2.is_file() {
            return;
        }

        let dir = tempfile::tempdir().unwrap();
        let fa_decompressed = dir.path().join("transcripts.fa");
        let tx2gene_decompressed = dir.path().join("tx2gene.tsv");

        // Decompress transcripts and tx2gene using qrnamap_cli::io::open_input
        {
            let mut reader = qrnamap_cli::io::open_input(dhingra_transcripts).unwrap();
            let mut out = fs::File::create(&fa_decompressed).unwrap();
            std::io::copy(&mut reader, &mut out).unwrap();
        }
        {
            let mut reader = qrnamap_cli::io::open_input(dhingra_tx2gene).unwrap();
            let mut out = fs::File::create(&tx2gene_decompressed).unwrap();
            std::io::copy(&mut reader, &mut out).unwrap();
        }

        let prep = commands::prepare_native_index(PrepareIndexRequest {
            index_path: None,
            index_bytes: None,
            index_gz_bytes: None,
            fasta_path: Some(fa_decompressed.to_string_lossy().to_string()),
            fasta_bytes: None,
        }).expect("should build native index from Dhingra transcripts");

        assert!(Path::new(&prep.index_path).is_file());

        // Create a fast subset of 1,000 read pairs from the real Dhingra FASTQs for quick unit test
        let sub_r1 = dir.path().join("dhingra_sub_R1.fq.gz");
        let sub_r2 = dir.path().join("dhingra_sub_R2.fq.gz");

        let slice_fastq = |src: &Path, dst: &Path, max_lines: usize| {
            let reader = qrnamap_cli::io::open_input(src).unwrap();
            let out_f = fs::File::create(dst).unwrap();
            let mut out_gz = flate2::write::GzEncoder::new(out_f, flate2::Compression::fast());
            use std::io::BufRead;
            for (idx, line) in reader.lines().enumerate() {
                if idx >= max_lines {
                    break;
                }
                let l = line.unwrap();
                use std::io::Write;
                writeln!(out_gz, "{}", l).unwrap();
            }
            out_gz.finish().unwrap();
        };

        slice_fastq(r1, &sub_r1, 4000); // 1,000 read pairs
        slice_fastq(r2, &sub_r2, 4000);

        let quant = commands::run_native_quant(NativeQuantRequest {
            name: "AF293_Azole_2_sub".to_string(),
            fastq1: vec![sub_r1.to_string_lossy().to_string()],
            fastq2: vec![sub_r2.to_string_lossy().to_string()],
            index_path: prep.index_path,
            tx2gene: tx2gene_decompressed.to_string_lossy().to_string(),
            strandedness: Some("unstranded".to_string()),
            fragment_mean: None,
            fragment_sd: None,
            threads: Some(2),
            out_dir: Some(dir.path().join("dhingra_out").to_string_lossy().to_string()),
        }).expect("should quantify Dhingra sample");

        assert_eq!(quant.n_processed, 1000);
        assert!(quant.n_pseudoaligned > 500, "pseudoaligned: {}", quant.n_pseudoaligned);
        assert!(!quant.gene_counts.is_empty(), "gene counts should be non-empty");
    }
}
