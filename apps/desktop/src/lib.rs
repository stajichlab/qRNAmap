use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

use serde::{Deserialize, Serialize};
use tauri::{WebviewUrl, WebviewWindowBuilder};

use qrnamap_core::{
    gene_counts_from_abundance, kallisto_quant, parse_tx2gene, zero_warning,
    Quantifier, RunConfig, Sample, Strandedness,
};

/// The kallisto version of the page's WebAssembly build (site/vendor/kallisto/PROVENANCE.md). A page-built index is only
/// handed to a native kallisto of the same major.minor version: kallisto index files are not compatible across versions.
pub const PAGE_KALLISTO_VERSION: &str = "0.52.0";

/// The native kallisto quant that is running now (at most one), so `cancel_native_quant` can stop it.
static RUNNING_QUANT: Mutex<Option<Child>> = Mutex::new(None);
/// Set by `cancel_native_quant`, so a killed run is reported as cancelled (not as a kallisto failure) on every OS.
static CANCEL_REQUESTED: AtomicBool = AtomicBool::new(false);

/// Information about native execution capabilities on this machine.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NativeCapabilities {
    pub available: bool,
    pub kallisto_path: Option<String>,
    pub kallisto_version: Option<String>,
    pub num_cpus: usize,
    pub temp_dir: String,
    /// True when the native kallisto has the same major.minor version as the page build, so a page-built index can be used.
    pub index_compatible: bool,
    /// A sentence for the page when the native kallisto cannot use the page-built index (or cannot run).
    pub warning: Option<String>,
}

/// `major.minor` of a kallisto version text such as "kallisto, version 0.52.0".
pub fn major_minor(version_text: &str) -> Option<(u32, u32)> {
    let v = version_text.split_whitespace().last()?;
    let mut parts = v.trim_start_matches('v').split('.');
    let major = parts.next()?.parse().ok()?;
    let minor = parts.next()?.parse().ok()?;
    Some((major, minor))
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

/// Whether the native kallisto can load an index, and why not.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NativeIndexCheck {
    pub ok: bool,
    pub detail: String,
}

/// Whether R with DESeq2 is on this computer (for "Run DESeq2 (R) on this computer").
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Deseq2Check {
    pub available: bool,
    pub rscript: Option<String>,
    pub version: Option<String>,
    pub detail: String,
}

/// Run DESeq2 (through R, as `qrnamap de` does) on the page's gene count matrix and sample sheet.
#[derive(Debug, Clone, Deserialize)]
pub struct LocalDeRequest {
    pub counts_tsv: String,
    pub samples_tsv: String,
    pub design: String,
    pub contrast: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocalDeResponse {
    /// `de.results.csv` as DESeq2 wrote it (gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj; NA for missing).
    pub results_csv: String,
    pub deseq2_version: String,
    pub rscript: String,
    /// The folder with every output (results, normalized counts, VST, versions) and the inputs, kept for the user.
    pub out_dir: String,
}

/// Write the page's counts and sheet to `dir` and run DESeq2 there with `rscript` (`qrnamap de`'s code). The outputs stay
/// in `dir`.
pub fn run_local_deseq2(req: &LocalDeRequest, rscript: &Path, dir: &Path) -> Result<LocalDeResponse, String> {
    fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    let counts = dir.join("counts.tsv");
    let samples = dir.join("samples.tsv");
    fs::write(&counts, &req.counts_tsv).map_err(|e| format!("cannot write {}: {e}", counts.display()))?;
    fs::write(&samples, &req.samples_tsv).map_err(|e| format!("cannot write {}: {e}", samples.display()))?;
    let args = qrnamap_cli::de::DeArgs {
        counts,
        samples,
        design: req.design.clone(),
        contrast: req.contrast.clone(),
        shrink: None,
        out: dir.join("de"),
        tmpdir: Some(dir.join("tmp")),
        rscript: Some(rscript.to_path_buf()),
    };
    let out = qrnamap_cli::de::run_de(&args).map_err(|e| e.to_string())?;
    let results = dir.join("de.results.csv");
    let results_csv = fs::read_to_string(&results).map_err(|e| format!("cannot read {}: {e}", results.display()))?;
    Ok(LocalDeResponse {
        results_csv,
        deseq2_version: out.deseq2_version,
        rscript: rscript.to_string_lossy().to_string(),
        out_dir: dir.to_string_lossy().to_string(),
    })
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
        ];
        for path in &user_paths {
            let p = PathBuf::from(path);
            if p.is_file() {
                return Some(p);
            }
        }
    }

    conda_kallisto_candidates(&home).into_iter().find(|p| p.is_file())
}

/// kallisto in conda, mamba and micromamba installs: each root's `bin/` and then its environments (`envs/*/bin/`, by
/// name). An app started from a desktop launcher (Linux) or the Dock (macOS) does not get the shell's PATH, so an
/// environment that the user activates in a terminal is not on it.
pub fn conda_kallisto_candidates(home: &str) -> Vec<PathBuf> {
    let mut roots: Vec<PathBuf> = Vec::new();
    for var in ["CONDA_PREFIX", "MAMBA_ROOT_PREFIX", "CONDA_ROOT"] {
        if let Ok(v) = std::env::var(var) {
            if !v.is_empty() {
                roots.push(PathBuf::from(v));
            }
        }
    }
    if !home.is_empty() {
        for dir in ["micromamba", "miniforge3", "mambaforge", "miniconda3", "anaconda3", ".local/share/mamba", ".conda"] {
            roots.push(Path::new(home).join(dir));
        }
    }
    roots.push(PathBuf::from("/opt/conda"));
    roots.push(PathBuf::from("/opt/miniforge3"));

    let mut out = Vec::new();
    for root in roots {
        out.push(root.join("bin").join("kallisto"));
        if let Ok(entries) = fs::read_dir(root.join("envs")) {
            let mut envs: Vec<PathBuf> = entries.flatten().map(|e| e.path()).collect();
            envs.sort();
            out.extend(envs.into_iter().map(|env| env.join("bin").join("kallisto")));
        }
    }
    out
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

/// Write index bytes to `target`, gunzipping them first when `gz` is set.
pub fn write_index_bytes(target: &Path, bytes: &[u8], gz: bool) -> Result<(), String> {
    if gz {
        let mut decoder = flate2::read::GzDecoder::new(bytes);
        let mut out = File::create(target).map_err(|e| format!("cannot create {}: {e}", target.display()))?;
        std::io::copy(&mut decoder, &mut out).map_err(|e| format!("cannot decompress index.gz: {e}"))?;
    } else {
        fs::write(target, bytes).map_err(|e| format!("cannot write {}: {e}", target.display()))?;
    }
    Ok(())
}

/// Can this `kallisto` load `index`?
///
/// kallisto index files are not portable between builds of the same version (0.52.0): an index built by kallisto.wasm
/// (the page, and the pre-built indexes on S3) loads in macOS Homebrew kallisto but crashes the Linux conda build
/// (`std::bad_array_new_length`), and a Linux conda index is "Corrupted index; unitig not found" for the other two
/// (docs/browser/page-v1-measurements.md). So the desktop app checks before it uses one: a quant of one made-up read.
/// The index loaded when kallisto prints its "[quant] processed" line; the exit status alone does not tell, because a read
/// that does not pseudoalign ends with exit 1 too. `work` gets the read and the (unused) output.
pub fn index_loads(kallisto: &Path, index: &Path, work: &Path) -> Result<(), String> {
    fs::create_dir_all(work).map_err(|e| format!("cannot create {}: {e}", work.display()))?;
    let read = work.join("one-read.fastq");
    let seq = "ACGT".repeat(15);
    fs::write(&read, format!("@check\n{seq}\n+\n{}\n", "I".repeat(seq.len())))
        .map_err(|e| format!("cannot write {}: {e}", read.display()))?;
    let output = Command::new(kallisto)
        .args(["quant", "-t", "1", "--single", "-l", "200", "-s", "20", "-i"])
        .arg(index)
        .arg("-o")
        .arg(work.join("out"))
        .arg(&read)
        .stdin(Stdio::null())
        .output()
        .map_err(|e| format!("failed to start kallisto: {e}"))?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    if stderr.contains("[quant] processed") {
        return Ok(());
    }
    // After a failed load kallisto may print its usage text: quote the error line, not the last line.
    let lines: Vec<&str> = stderr.lines().map(str::trim).filter(|l| !l.is_empty()).collect();
    let last = lines
        .iter()
        .rev()
        .find(|l| l.contains("Error") || l.contains("what():") || l.contains("Corrupted"))
        .or(lines.last())
        .copied()
        .unwrap_or("");
    let how = match output.status.code() {
        Some(code) => format!("exit {code}"),
        None => "stopped by a signal".to_string(),
    };
    Err(format!("{} cannot read this index ({how}{}{last})", kallisto.display(), if last.is_empty() { "" } else { ": " }))
}

/// A fresh working folder for one `kallisto index` run, next to the index it writes. kallisto 0.50+ makes its temporary
/// folder `tmp` in the working directory, and fails with "could not create directory tmp" when that directory is not
/// writable (an app started from a Linux desktop launcher often runs in `/`) or when another build is using `tmp`.
/// The index and FASTA paths given to kallisto must be absolute. The folder is removed when the value is dropped.
pub fn index_work_dir(target: &Path) -> Result<tempfile::TempDir, String> {
    let parent = target.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|e| format!("cannot create {}: {e}", parent.display()))?;
    tempfile::Builder::new()
        .prefix(".kallisto-index-")
        .tempdir_in(parent)
        .map_err(|e| format!("cannot make a working folder in {}: {e}", parent.display()))
}

/// Where the desktop app keeps the native index of a genome: one file per genome key and kallisto version, so an index is
/// built once per machine. The key is the genome accession (or "upload"), cut to safe characters.
pub fn native_index_path(dir: &Path, key: &str, kallisto_version: &str) -> PathBuf {
    let safe = |s: &str| -> String {
        let t: String = s.chars().map(|c| if c.is_ascii_alphanumeric() || c == '.' || c == '-' { c } else { '_' }).collect();
        t.trim_matches(|c| c == '.' || c == '_').chars().take(80).collect()
    };
    let version = kallisto_version.rsplit(' ').next().unwrap_or(kallisto_version);
    let key = match safe(key) { k if k.is_empty() => "genome".to_string(), k => k };
    dir.join(format!("{key}.kallisto-{}.idx", safe(version)))
}

/// Build a native index of `fasta` at `target` with `kallisto index`, then check that it loads. The index is written to a
/// temporary name and renamed at the end, so a stopped build never leaves a file that looks finished.
pub fn build_native_index_file(kallisto: &Path, fasta: &Path, target: &Path, threads: usize) -> Result<u64, String> {
    let partial = target.with_extension("idx.partial");
    let cwd = index_work_dir(target)?;
    let output = Command::new(kallisto)
        .current_dir(cwd.path())
        .arg("index")
        .arg("-t")
        .arg(threads.max(1).to_string())
        .arg("-i")
        .arg(&partial)
        .arg(fasta)
        .stdin(Stdio::null())
        .output()
        .map_err(|e| format!("failed to start kallisto index: {e}"))?;
    if !output.status.success() || !partial.is_file() {
        let _ = fs::remove_file(&partial);
        let stderr = String::from_utf8_lossy(&output.stderr);
        let tail: Vec<&str> = stderr.lines().filter(|l| !l.trim().is_empty()).collect();
        return Err(format!("kallisto index failed: {}", tail[tail.len().saturating_sub(3)..].join(" ")));
    }
    fs::rename(&partial, target).map_err(|e| format!("cannot rename {}: {e}", partial.display()))?;
    let work = target.with_extension("check");
    let checked = index_loads(kallisto, target, &work);
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = checked {
        let _ = fs::remove_file(target);
        return Err(format!("the index built by this kallisto does not load in it: {e}"));
    }
    Ok(fs::metadata(target).map(|m| m.len()).unwrap_or(0))
}

/// The folder of the desktop app's scratch files (the page's index, native indexes, checks).
pub fn scratch_dir() -> Result<PathBuf, String> {
    let dir = std::env::temp_dir().join("qrnamap_desktop");
    fs::create_dir_all(&dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    Ok(dir)
}

/// Run `cmd` to the end, keeping its `Child` in `RUNNING_QUANT` so `cancel_native_quant` can kill it.
/// stdout and stderr are read on their own threads (a full pipe would otherwise stall kallisto).
/// Returns Err("Cancelled ...") when the process was killed through `cancel_native_quant`.
pub fn run_cancellable(mut cmd: Command) -> Result<std::process::Output, String> {
    use std::io::Read;
    cmd.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = cmd.spawn().map_err(|e| format!("failed to start kallisto: {e}"))?;
    let mut out_pipe = child.stdout.take();
    let mut err_pipe = child.stderr.take();
    let out_thread = std::thread::spawn(move || {
        let mut buf = Vec::new();
        if let Some(p) = out_pipe.as_mut() { let _ = p.read_to_end(&mut buf); }
        buf
    });
    let err_thread = std::thread::spawn(move || {
        let mut buf = Vec::new();
        if let Some(p) = err_pipe.as_mut() { let _ = p.read_to_end(&mut buf); }
        buf
    });
    {
        let mut guard = RUNNING_QUANT.lock().unwrap_or_else(|e| e.into_inner());
        if guard.is_some() {
            let _ = child.kill();
            let _ = child.wait();
            return Err("another native kallisto run is still going".to_string());
        }
        CANCEL_REQUESTED.store(false, Ordering::SeqCst);
        *guard = Some(child);
    }
    let status = loop {
        {
            let mut guard = RUNNING_QUANT.lock().unwrap_or_else(|e| e.into_inner());
            let child = guard.as_mut().expect("RUNNING_QUANT holds the child until it ends");
            match child.try_wait() {
                Ok(Some(status)) => { *guard = None; break status; }
                Ok(None) => {}
                Err(e) => { *guard = None; return Err(format!("cannot wait for kallisto: {e}")); }
            }
        }
        std::thread::sleep(std::time::Duration::from_millis(200));
    };
    let stdout = out_thread.join().unwrap_or_default();
    let stderr = err_thread.join().unwrap_or_default();
    if CANCEL_REQUESTED.swap(false, Ordering::SeqCst) {
        return Err("Cancelled: kallisto was stopped.".to_string());
    }
    Ok(std::process::Output { status, stdout, stderr })
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

        let page = major_minor(PAGE_KALLISTO_VERSION);
        let native = kallisto_version.as_deref().and_then(major_minor);
        let index_compatible = available && native.is_some() && native == page;
        let warning = if !available {
            Some(format!(
                "No kallisto binary was found (looked in QRNAMAP_KALLISTO, PATH, Homebrew, and conda, mamba, micromamba \
                 and pixi folders). Install kallisto 0.52 ({}) or set QRNAMAP_KALLISTO to its path.",
                if cfg!(target_os = "macos") {
                    "for example `brew install kallisto`"
                } else {
                    "for example `micromamba create -n kallisto -c conda-forge -c bioconda kallisto=0.52`"
                }
            ))
        } else if !index_compatible {
            Some(format!(
                "The native kallisto ({}) is not version {PAGE_KALLISTO_VERSION}, so it cannot read the index this app builds. \
                 Install kallisto 0.52 or set QRNAMAP_KALLISTO to a 0.52 binary.",
                kallisto_version.as_deref().unwrap_or("unknown version")
            ))
        } else {
            None
        };

        NativeCapabilities {
            available,
            kallisto_path,
            kallisto_version,
            num_cpus,
            temp_dir,
            index_compatible,
            warning,
        }
    }

    /// Write the kallisto index that the page built (or the gzip of a prebuilt one) to disk for native kallisto.
    /// The bytes come as the raw IPC body (a Uint8Array on the page side), not as a JSON number array: an index is
    /// 100 to 300 MB, and JSON would make it several times bigger and slow. Header `x-qrnamap-gzip: 1` marks a gzip body.
    #[tauri::command(async)]
    pub fn write_native_index(request: tauri::ipc::Request<'_>) -> Result<PrepareIndexResponse, String> {
        let tauri::ipc::InvokeBody::Raw(bytes) = request.body() else {
            return Err("write_native_index needs the index bytes as a raw (Uint8Array) body".to_string());
        };
        let gz = request
            .headers()
            .get("x-qrnamap-gzip")
            .and_then(|v| v.to_str().ok())
            == Some("1");
        let scratch_dir = std::env::temp_dir().join("qrnamap_desktop");
        fs::create_dir_all(&scratch_dir).map_err(|e| format!("cannot create scratch dir: {e}"))?;
        let target_idx = scratch_dir.join("transcripts.idx");
        write_index_bytes(&target_idx, bytes, gz)?;
        let size = fs::metadata(&target_idx).map(|m| m.len()).unwrap_or(0);
        Ok(PrepareIndexResponse { index_path: target_idx.to_string_lossy().to_string(), index_size: size })
    }

    /// Check that the native kallisto can load the index at `index_path` (see `index_loads`).
    #[tauri::command(async)]
    pub fn check_native_index(index_path: String) -> NativeIndexCheck {
        let Some(kallisto) = find_kallisto() else {
            return NativeIndexCheck { ok: false, detail: "no kallisto binary was found".to_string() };
        };
        let work = match scratch_dir() { Ok(d) => d.join("check"), Err(e) => return NativeIndexCheck { ok: false, detail: e } };
        let result = index_loads(&kallisto, &resolve_path_str(&index_path, None), &work);
        let _ = fs::remove_dir_all(&work);
        match result {
            Ok(()) => NativeIndexCheck { ok: true, detail: String::new() },
            Err(detail) => NativeIndexCheck { ok: false, detail },
        }
    }

    /// The native index this machine built before for genome `key` with the current kallisto, if it is there and loads.
    #[tauri::command(async)]
    pub fn cached_native_index(key: String) -> Option<PrepareIndexResponse> {
        let kallisto = find_kallisto()?;
        let version = get_kallisto_version(&kallisto).ok()?;
        let path = native_index_path(&scratch_dir().ok()?.join("native"), &key, &version);
        if !path.is_file() {
            return None;
        }
        let work = path.with_extension("check");
        let ok = index_loads(&kallisto, &path, &work).is_ok();
        let _ = fs::remove_dir_all(&work);
        if !ok {
            let _ = fs::remove_file(&path);
            return None;
        }
        let size = fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
        Some(PrepareIndexResponse { index_path: path.to_string_lossy().to_string(), index_size: size })
    }

    /// Build a native kallisto index from the transcript FASTA in the raw IPC body (a Uint8Array; header
    /// `x-qrnamap-key` names the genome, so the index is kept for the next time). Used when the native kallisto cannot read
    /// the page's kallisto.wasm index.
    #[tauri::command(async)]
    pub fn build_native_index(request: tauri::ipc::Request<'_>) -> Result<PrepareIndexResponse, String> {
        let tauri::ipc::InvokeBody::Raw(fasta) = request.body() else {
            return Err("build_native_index needs the transcript FASTA as a raw (Uint8Array) body".to_string());
        };
        if !fasta.starts_with(b">") {
            return Err("build_native_index: the body is not a FASTA file".to_string());
        }
        let key = request.headers().get("x-qrnamap-key").and_then(|v| v.to_str().ok()).unwrap_or("upload").to_string();
        let kallisto = find_kallisto().ok_or_else(|| "no kallisto binary was found".to_string())?;
        let version = get_kallisto_version(&kallisto)?;
        let dir = scratch_dir()?.join("native");
        fs::create_dir_all(&dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
        let target = native_index_path(&dir, &key, &version);
        let fasta_path = target.with_extension("transcripts.fa");
        fs::write(&fasta_path, fasta).map_err(|e| format!("cannot write {}: {e}", fasta_path.display()))?;
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(8);
        let built = build_native_index_file(&kallisto, &fasta_path, &target, threads);
        let _ = fs::remove_file(&fasta_path);
        let size = built?;
        Ok(PrepareIndexResponse { index_path: target.to_string_lossy().to_string(), index_size: size })
    }

    /// Is R with DESeq2 on this computer? Tries each Rscript found (`qrnamap_cli::de::find_deseq2`).
    #[tauri::command(async)]
    pub fn check_deseq2() -> Deseq2Check {
        match qrnamap_cli::de::find_deseq2() {
            Ok((r, v)) => Deseq2Check { available: true, rscript: Some(r.to_string_lossy().to_string()), version: Some(v), detail: String::new() },
            Err(detail) => Deseq2Check { available: false, rscript: None, version: None, detail },
        }
    }

    /// Run DESeq2 through R on the page's counts and sheet. The outputs go to a new folder under the app's scratch folder.
    #[tauri::command(async)]
    pub fn run_deseq2_local(req: LocalDeRequest) -> Result<LocalDeResponse, String> {
        let (rscript, _) = qrnamap_cli::de::find_deseq2().map_err(|e| format!("R with DESeq2 is not available: {e}"))?;
        let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
        let dir = scratch_dir()?.join(format!("deseq2-{stamp}"));
        run_local_deseq2(&req, &rscript, &dir)
    }

    /// Stop the native kallisto quant that is running now, if any. The waiting `run_native_quant` then returns an error
    /// that starts with "Cancelled".
    #[tauri::command]
    pub fn cancel_native_quant() -> bool {
        let mut guard = RUNNING_QUANT.lock().unwrap_or_else(|e| e.into_inner());
        match guard.as_mut() {
            Some(child) => {
                CANCEL_REQUESTED.store(true, Ordering::SeqCst);
                child.kill().is_ok()
            }
            None => false,
        }
    }

    #[tauri::command(async)]
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

            let target_idx = std::path::absolute(&target_idx).unwrap_or(target_idx);
            let fasta_file = std::path::absolute(&fasta_file).unwrap_or(fasta_file);
            let cwd = index_work_dir(&target_idx)?;
            let output = Command::new(&kallisto)
                .current_dir(cwd.path())
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

    // `async`: a full-depth run takes minutes; a plain (sync) command would run on the main thread and freeze the window.
    #[tauri::command(async)]
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

        let output = run_cancellable(cmd).map_err(|e| format!("sample '{}': {e}", req.name))?;

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

    #[tauri::command(async)]
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

    /// Pick several files at once (the FASTQ files of a sample sheet). Returns their absolute paths (empty when cancelled).
    #[tauri::command]
    pub fn pick_files_dialog(title: Option<String>, filter_name: Option<String>, filter_exts: Option<Vec<String>>) -> Result<Vec<String>, String> {
        let mut dialog = rfd::FileDialog::new();
        if let Some(t) = title {
            dialog = dialog.set_title(&t);
        }
        if let (Some(name), Some(exts)) = (filter_name, filter_exts) {
            let str_refs: Vec<&str> = exts.iter().map(|s| s.as_str()).collect();
            dialog = dialog.add_filter(&name, &str_refs);
        }
        Ok(dialog
            .pick_files()
            .unwrap_or_default()
            .into_iter()
            .map(|p| p.to_string_lossy().to_string())
            .collect())
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
            commands::pick_files_dialog,
            commands::write_native_index,
            commands::check_native_index,
            commands::cached_native_index,
            commands::build_native_index,
            commands::check_deseq2,
            commands::run_deseq2_local,
            commands::cancel_native_quant,
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

    /// Tests that start kallisto share RUNNING_QUANT (one native run at a time), so they take this lock. Tests that write a
    /// script and run it take it too: on Linux, exec of a file fails with "Text file busy" (ETXTBSY) while a child that
    /// another test thread is forking still holds the file open for writing.
    static NATIVE_RUN_LOCK: Mutex<()> = Mutex::new(());
    fn native_run_lock() -> std::sync::MutexGuard<'static, ()> {
        NATIVE_RUN_LOCK.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// The kallisto of this machine for the tests that run it. None skips the test (with a note), unless
    /// QRNAMAP_REQUIRE_TOOLS=1 (CI), where a missing kallisto is a failure (ruling R5, as crates/qrnamap-core/tests).
    fn kallisto_for_test() -> Option<PathBuf> {
        let k = find_kallisto();
        if k.is_none() {
            if std::env::var("QRNAMAP_REQUIRE_TOOLS").map(|v| v == "1").unwrap_or(false) {
                panic!("kallisto was not found and QRNAMAP_REQUIRE_TOOLS=1 (set QRNAMAP_KALLISTO or put it on PATH)");
            }
            eprintln!("SKIP: kallisto not found (set QRNAMAP_REQUIRE_TOOLS=1 to make this a failure)");
        }
        k
    }

    #[test]
    fn test_major_minor() {
        assert_eq!(major_minor("kallisto, version 0.52.0"), Some((0, 52)));
        assert_eq!(major_minor("0.46.1"), Some((0, 46)));
        assert_eq!(major_minor("kallisto (unknown version)"), None);
        assert_eq!(major_minor(PAGE_KALLISTO_VERSION), Some((0, 52)));
    }

    #[test]
    fn test_write_index_bytes_plain_and_gzip() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("plain.idx");
        write_index_bytes(&plain, b"index-bytes", false).unwrap();
        assert_eq!(fs::read(&plain).unwrap(), b"index-bytes");

        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        enc.write_all(b"gzipped-index").unwrap();
        let gz = enc.finish().unwrap();
        let unz = dir.path().join("gz.idx");
        write_index_bytes(&unz, &gz, true).unwrap();
        assert_eq!(fs::read(&unz).unwrap(), b"gzipped-index");
        assert!(write_index_bytes(&dir.path().join("bad.idx"), b"not gzip", true).is_err());
    }

    #[cfg(unix)]
    #[test]
    fn test_run_cancellable_output_and_cancel() {
        let _lock = native_run_lock();
        let mut cmd = Command::new("sh");
        cmd.args(["-c", "echo out; echo err 1>&2"]);
        let out = run_cancellable(cmd).unwrap();
        assert!(out.status.success());
        assert_eq!(String::from_utf8_lossy(&out.stdout).trim(), "out");
        assert_eq!(String::from_utf8_lossy(&out.stderr).trim(), "err");

        let canceller = std::thread::spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(400));
            commands::cancel_native_quant()
        });
        let mut slow = Command::new("sleep");
        slow.arg("30");
        let started = std::time::Instant::now();
        let res = run_cancellable(slow);
        assert!(canceller.join().unwrap(), "cancel found the running process");
        assert!(res.unwrap_err().starts_with("Cancelled"));
        assert!(started.elapsed() < std::time::Duration::from_secs(10));
    }

    #[test]
    fn test_find_kallisto() {
        let Some(path) = kallisto_for_test() else { return };
        assert!(path.is_file(), "kallisto path must be an existing file");
        let version = get_kallisto_version(&path).expect("should get version");
        assert!(version.to_lowercase().contains("kallisto"), "version string should contain 'kallisto'");
    }

    #[test]
    fn conda_candidates_cover_root_bins_and_named_envs() {
        let home = tempfile::tempdir().unwrap();
        let h = home.path();
        for env in ["zeta", "kallisto"] {
            fs::create_dir_all(h.join("micromamba/envs").join(env).join("bin")).unwrap();
        }
        fs::create_dir_all(h.join("miniforge3/envs/rnaseq/bin")).unwrap();
        let c = conda_kallisto_candidates(&h.to_string_lossy());
        let pos = |rel: &str| c.iter().position(|p| p == &h.join(rel)).unwrap_or_else(|| panic!("{rel} missing from {c:?}"));
        // each root's bin first, then its environments sorted by name; roots in the documented order
        assert!(pos("micromamba/bin/kallisto") < pos("micromamba/envs/kallisto/bin/kallisto"));
        assert!(pos("micromamba/envs/kallisto/bin/kallisto") < pos("micromamba/envs/zeta/bin/kallisto"));
        assert!(pos("micromamba/envs/zeta/bin/kallisto") < pos("miniforge3/bin/kallisto"));
        assert!(pos("miniforge3/bin/kallisto") < pos("miniforge3/envs/rnaseq/bin/kallisto"));
        assert!(c.contains(&PathBuf::from("/opt/conda/bin/kallisto")));
    }

    #[test]
    fn test_check_capabilities() {
        let _lock = native_run_lock();
        let caps = commands::check_native_capabilities();
        assert!(caps.num_cpus >= 1);
        if kallisto_for_test().is_none() {
            // Without kallisto the app says so and offers the page's own kallisto.wasm.
            assert!(!caps.available && caps.kallisto_path.is_none());
            assert!(caps.warning.as_deref().unwrap_or("").contains("QRNAMAP_KALLISTO"), "{:?}", caps.warning);
            return;
        }
        assert!(caps.available);
        assert!(caps.kallisto_path.is_some());
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

    /// A stand-in kallisto (a shell script): `index` writes the -i file, `quant` prints what `quant_stderr` says and exits
    /// with `quant_exit` (a negative value aborts with SIGABRT, like the conda build on a kallisto.wasm index).
    #[cfg(unix)]
    fn fake_kallisto(dir: &Path, quant_stderr: &str, quant_exit: i32, index_ok: bool) -> PathBuf {
        use std::os::unix::fs::PermissionsExt;
        let bin = dir.join("kallisto");
        let quant_end = if quant_exit < 0 { "kill -ABRT $$".to_string() } else { format!("exit {quant_exit}") };
        let index_body = if index_ok {
            "while [ $# -gt 0 ]; do if [ \"$1\" = -i ]; then shift; printf idx > \"$1\"; fi; shift; done; exit 0"
        } else {
            "echo '[build] reading transcripts' >&2; echo 'Error: bad FASTA' >&2; exit 1"
        };
        fs::write(&bin, format!("#!/bin/sh\ncase \"$1\" in\n  index) {index_body} ;;\n  quant) printf '%b\\n' '{quant_stderr}' >&2; {quant_end} ;;\nesac\n")).unwrap();
        fs::set_permissions(&bin, fs::Permissions::from_mode(0o755)).unwrap();
        bin
    }

    #[cfg(unix)]
    #[test]
    fn index_loads_reads_the_processed_line_not_the_exit_status() {
        let _lock = native_run_lock();
        let dir = tempfile::tempdir().unwrap();
        let idx = dir.path().join("t.idx");
        fs::write(&idx, b"idx").unwrap();
        // A readable index: one made-up read, zero pseudoaligned, exit 1. Still a pass.
        let ok = fake_kallisto(dir.path(), "[quant] processed 1 reads, 0 reads pseudoaligned", 1, true);
        index_loads(&ok, &idx, &dir.path().join("w1")).expect("the processed line means the index loaded");

        let sub = dir.path().join("b");
        fs::create_dir_all(&sub).unwrap();
        let crash = fake_kallisto(&sub, "terminate called after throwing an instance of std::bad_array_new_length", -1, true);
        let e = index_loads(&crash, &idx, &dir.path().join("w2")).unwrap_err();
        assert!(e.contains("cannot read this index") && e.contains("bad_array_new_length"), "{e}");

        let sub = dir.path().join("c");
        fs::create_dir_all(&sub).unwrap();
        let corrupt = fake_kallisto(&sub, "Error: Corrupted index; unitig not found: ACGT\n  --verbose   Print out progress information every 1M proccessed reads", 1, true);
        let e = index_loads(&corrupt, &idx, &dir.path().join("w3")).unwrap_err();
        assert!(e.contains("exit 1") && e.contains("Corrupted index") && !e.contains("proccessed"), "{e}");
    }

    #[test]
    fn native_index_path_is_per_genome_and_kallisto_version() {
        let d = Path::new("/x");
        assert_eq!(native_index_path(d, "GCF_000146045.2", "kallisto, version 0.52.0"), d.join("GCF_000146045.2.kallisto-0.52.0.idx"));
        assert_eq!(native_index_path(d, "../../etc/passwd", "0.52.0"), d.join("etc_passwd.kallisto-0.52.0.idx"));
        assert_eq!(native_index_path(d, "", "0.51.1"), d.join("genome.kallisto-0.51.1.idx"));
    }

    #[cfg(unix)]
    #[test]
    fn build_native_index_file_renames_only_a_finished_checked_index() {
        let _lock = native_run_lock();
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("t.fa");
        fs::write(&fa, ">t\nACGT\n").unwrap();
        let good = fake_kallisto(dir.path(), "[quant] processed 1 reads", 1, true);
        let target = dir.path().join("g.kallisto-0.52.0.idx");
        build_native_index_file(&good, &fa, &target, 2).expect("built and checked");
        assert!(target.is_file());
        assert!(!target.with_extension("idx.partial").exists());

        // kallisto index fails: no file is left, and the error has kallisto's last lines.
        let sub = dir.path().join("bad");
        fs::create_dir_all(&sub).unwrap();
        let bad = fake_kallisto(&sub, "[quant] processed 1 reads", 1, false);
        let t2 = dir.path().join("h.idx");
        let e = build_native_index_file(&bad, &fa, &t2, 1).unwrap_err();
        assert!(e.contains("bad FASTA"), "{e}");
        assert!(!t2.exists() && !t2.with_extension("idx.partial").exists());

        // Built, but it does not load in the same kallisto: removed.
        let sub = dir.path().join("noload");
        fs::create_dir_all(&sub).unwrap();
        let noload = fake_kallisto(&sub, "Error: Corrupted index", 1, true);
        let t3 = dir.path().join("i.idx");
        assert!(build_native_index_file(&noload, &fa, &t3, 1).is_err());
        assert!(!t3.exists());
    }

    /// With a real kallisto on this machine: an index it builds passes the check, and a file that is not an index fails it.
    #[test]
    fn real_kallisto_checks_its_own_index_and_refuses_junk() {
        let _lock = native_run_lock();
        let Some(kallisto) = kallisto_for_test() else { return };
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("t.fa");
        let tx1 = "ACGTTGCAAGCTTGACCGTAGGCTAACGTTAGCCATGGTACCGATCGATTGCAAGGCTAGCTAGGCATCGATCGGATCCTAGCTAAGCTTGCA";
        let tx2 = "TTGACCAGTCGGATCATGCAGTCCAGTACGATCAGGCTTAGCATGCAGGTCAAGTCGATCGGCTAGCAGTCAGCTAGGCAATCGACTGACTGA";
        fs::write(&fa, format!(">tx1\n{tx1}\n>tx2\n{tx2}\n")).unwrap();
        let target = dir.path().join("real.idx");
        build_native_index_file(&kallisto, &fa, &target, 1).expect("kallisto builds an index it can load");
        let junk = dir.path().join("junk.idx");
        fs::write(&junk, vec![7u8; 4096]).unwrap();
        assert!(index_loads(&kallisto, &junk, &dir.path().join("w")).is_err());
    }

    /// With R and DESeq2 on this computer (set QRNAMAP_RSCRIPT, e.g. to the repo's pixi `de` environment): DESeq2 runs on a
    /// small simulated count matrix and the results CSV has a row per gene. Skipped without R and DESeq2.
    #[test]
    fn local_deseq2_runs_on_page_style_inputs() {
        let Ok((rscript, version)) = qrnamap_cli::de::find_deseq2() else { return };
        // Counts with real-looking spread (a fixed LCG, so the matrix is the same every run): mean 20 to 920, a gamma-like
        // multiplier per sample (DESeq2 refuses data whose dispersions are all near the minimum), 1 gene in 10 up 4-fold.
        let mut state: u64 = 42;
        let mut unif = move || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((state >> 11) as f64) / ((1u64 << 53) as f64)
        };
        let mut counts = String::from("gene\ta1\ta2\ta3\tb1\tb2\tb3\n");
        for g in 0..300u64 {
            let base = (20 + (g * 37) % 900) as f64;
            let fold = if g % 10 == 0 { 4.0 } else { 1.0 };
            let row: Vec<String> = (0..6u64)
                .map(|j| {
                    let mult = (0..4).map(|_| unif()).sum::<f64>() / 2.0 + 0.2 * unif(); // about 0.2 to 2.2
                    let mean = if j < 3 { base } else { base * fold };
                    ((mean * mult).round() as u64).to_string()
                })
                .collect();
            counts.push_str(&format!("gene{g}\t{}\n", row.join("\t")));
        }
        let samples = "sample\tfastq1\tfastq2\tcondition\na1\ta1.fq\t\tctl\na2\ta2.fq\t\tctl\na3\ta3.fq\t\tctl\nb1\tb1.fq\t\ttrt\nb2\tb2.fq\t\ttrt\nb3\tb3.fq\t\ttrt\n";
        let dir = tempfile::tempdir().unwrap();
        let req = LocalDeRequest { counts_tsv: counts, samples_tsv: samples.to_string(), design: "~ condition".into(), contrast: "condition,trt,ctl".into() };
        let r = run_local_deseq2(&req, &rscript, &dir.path().join("de")).expect("DESeq2 ran");
        assert_eq!(r.deseq2_version, version);
        let lines: Vec<&str> = r.results_csv.lines().collect();
        assert_eq!(lines.len(), 301, "header and one row per gene");
        assert!(lines[0].contains("log2FoldChange") && lines[0].contains("padj"), "{}", lines[0]);
        assert!(Path::new(&r.out_dir).join("de.normalized.csv").is_file());
    }

    #[test]
    fn test_native_quant_flow() {
        let _lock = native_run_lock();
        if kallisto_for_test().is_none() {
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
        let _lock = native_run_lock();
        if kallisto_for_test().is_none() {
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
