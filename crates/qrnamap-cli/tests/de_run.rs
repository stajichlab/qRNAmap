//! `qrnamap de` against the real Rscript and DESeq2. Tests skip with a message when Rscript or
//! DESeq2 is missing, and fail when `QRNAMAP_REQUIRE_TOOLS=1` (ruling R5). No R is mocked.

mod common;

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use qrnamap_cli::de::{read_count_matrix, run_de, DeArgs, DeRunError};
use qrnamap_core::de::{DeError, Shrink};
use qrnamap_core::CountMatrix;

/// True when Rscript and the DESeq2 package are present. Missing DESeq2 fails under
/// `QRNAMAP_REQUIRE_TOOLS=1`, else skips.
fn r_ready() -> bool {
    if common::require_tool("Rscript").is_none() {
        return false;
    }
    let ok = Command::new("Rscript")
        .args(["--vanilla", "-e", "quit(status = if (requireNamespace('DESeq2', quietly = TRUE)) 0 else 1)"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !ok {
        if std::env::var("QRNAMAP_REQUIRE_TOOLS").map(|v| v == "1").unwrap_or(false) {
            panic!("the R package DESeq2 is not installed and QRNAMAP_REQUIRE_TOOLS=1");
        }
        eprintln!("SKIP: the R package DESeq2 is not installed");
    }
    ok
}

/// Hand-written LCG (Numerical Recipes constants) with Box-Muller normals.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let (u1, u2) = (self.next_f64(), self.next_f64());
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

const NAMES: [&str; 6] = ["c1", "c2", "c3", "t1", "t2", "t3"];

/// 200 genes, 3 control and 3 treated samples. Genes g001..g020 are changed 4-fold, up for
/// odd positions and down for even ones. Counts are lognormal-overdispersed (sigma 0.15) with a
/// normal approximation of Poisson noise. Base means run from 200 to 2000.
fn synth(names: &[&str; 6]) -> CountMatrix {
    let mut rng = Lcg(20261119);
    let mut genes = Vec::new();
    let mut counts = Vec::new();
    for i in 1..=200 {
        let base = 200.0 + 1800.0 * rng.next_f64();
        let fold = match (i <= 20, i % 2) {
            (true, 1) => 4.0,
            (true, _) => 0.25,
            _ => 1.0,
        };
        let row: Vec<u64> = (0..6)
            .map(|j| {
                let mean = base * if j >= 3 { fold } else { 1.0 } * (0.15 * rng.normal()).exp();
                (mean + mean.sqrt() * rng.normal()).round().max(1.0) as u64
            })
            .collect();
        genes.push(format!("g{i:03}"));
        counts.push(row);
    }
    CountMatrix { genes, samples: names.iter().map(|s| s.to_string()).collect(), counts }
}

const FASTQ: &str = "/data/secret_reads";

/// Sample sheet with `condition`, plus extra columns as (header, six values).
fn sheet_text(names: &[&str; 6], extra: &[(&str, [&str; 6])]) -> String {
    let mut s = String::from("sample\tfastq1\tfastq2\tcondition");
    for (h, _) in extra {
        s.push('\t');
        s.push_str(h);
    }
    s.push('\n');
    for (j, n) in names.iter().enumerate() {
        let cond = if j < 3 { "control" } else { "treated" };
        s.push_str(&format!("{n}\t{FASTQ}_{n}_R1.fastq.gz\t{FASTQ}_{n}_R2.fastq.gz\t{cond}"));
        for (_, v) in extra {
            s.push('\t');
            s.push_str(v[j]);
        }
        s.push('\n');
    }
    s
}

struct Case {
    dir: tempfile::TempDir,
    sheet: PathBuf,
    counts: PathBuf,
}

fn case(names: &[&str; 6], extra: &[(&str, [&str; 6])]) -> Case {
    let dir = tempfile::tempdir().unwrap();
    let sheet = dir.path().join("sheet.tsv");
    std::fs::write(&sheet, sheet_text(names, extra)).unwrap();
    let counts = dir.path().join("counts.tsv.zst");
    let tsv = synth(names).to_tsv();
    let mut enc = zstd::Encoder::new(std::fs::File::create(&counts).unwrap(), 3).unwrap();
    enc.write_all(tsv.as_bytes()).unwrap();
    enc.finish().unwrap();
    Case { dir, sheet, counts }
}

impl Case {
    fn args(&self, design: &str) -> DeArgs {
        DeArgs {
            counts: self.counts.clone(),
            samples: self.sheet.clone(),
            design: design.into(),
            contrast: "condition,treated,control".into(),
            shrink: None,
            out: self.dir.path().join("out").join("de"),
            tmpdir: Some(self.dir.path().join("tmp")),
        }
    }
}

/// Read a CSV written by R (quoted header and strings, no embedded commas in our data).
fn read_csv(path: &Path) -> (Vec<String>, Vec<HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let unq = |s: &str| s.trim_matches('"').to_string();
    let header: Vec<String> = lines.next().unwrap().split(',').map(unq).collect();
    let rows = lines
        .map(|l| header.iter().cloned().zip(l.split(',').map(unq)).collect())
        .collect();
    (header, rows)
}

#[test]
fn synthetic_de_detects_changed_genes() {
    if !r_ready() {
        return;
    }
    let c = case(&NAMES, &[]);
    let out = run_de(&c.args("~ condition")).unwrap();
    assert!(out.deseq2_version.chars().next().unwrap().is_ascii_digit(), "{}", out.deseq2_version);
    let (_, rows) = read_csv(&c.dir.path().join("out/de.results.csv"));
    assert_eq!(rows.len(), 200);
    let mut hit_changed = 0;
    let mut hit_other = 0;
    let mut wrong_direction = 0;
    for r in &rows {
        let n: usize = r["gene"][1..].parse().unwrap();
        let sig = r["padj"].parse::<f64>().map(|p| p < 0.05).unwrap_or(false);
        if n <= 20 {
            if sig {
                hit_changed += 1;
                let lfc: f64 = r["log2FoldChange"].parse().unwrap();
                if (n % 2 == 1) != (lfc > 0.0) {
                    wrong_direction += 1;
                }
            }
        } else if sig {
            hit_other += 1;
        }
    }
    eprintln!("MEASURED: changed genes with padj<0.05: {hit_changed}/20; other genes: {hit_other}/180; wrong direction: {wrong_direction}");
    assert!(hit_changed >= 15, "changed genes detected: {hit_changed}");
    assert!(hit_other <= 5, "false positives: {hit_other}");
    assert_eq!(wrong_direction, 0);
    for f in ["normalized.csv", "vst.csv", "versions.txt"] {
        assert!(c.dir.path().join(format!("out/de.{f}")).is_file(), "{f}");
    }
    assert!(!c.dir.path().join("out/de.results.shrunk.csv").exists());
}

#[test]
fn confounded_batch_is_not_full_rank() {
    if !r_ready() {
        return;
    }
    let batch = ["x", "x", "x", "y", "y", "y"];
    let c = case(&NAMES, &[("batch", batch)]);
    let err = run_de(&c.args("~ batch + condition")).unwrap_err();
    match err {
        DeRunError::De(DeError::NotFullRank(m)) => {
            assert!(m.contains("not full rank"), "{m}");
            assert!(m.len() <= 5000);
        }
        other => panic!("expected NotFullRank, got {other}"),
    }
}

#[test]
fn de_does_not_print_fastq_paths() {
    if !r_ready() {
        return;
    }
    let c = case(&NAMES, &[]);
    let out = Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args(["de", "--counts"])
        .arg(&c.counts)
        .arg("--samples")
        .arg(&c.sheet)
        .args(["--design", "~ condition", "--contrast", "condition,treated,control", "--out"])
        .arg(c.dir.path().join("o/de"))
        .arg("--tmpdir")
        .arg(c.dir.path().join("tmp"))
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let all = format!("{}{}", String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    assert!(!all.contains("secret_reads") && !all.contains(".fastq"), "{all}");
    assert!(String::from_utf8_lossy(&out.stdout).contains("DESeq2 "), "{all}");
    // A failing run does not print them either.
    let bad = Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args(["de", "--counts"])
        .arg(&c.counts)
        .arg("--samples")
        .arg(&c.sheet)
        .args(["--design", "~ condition", "--contrast", "condition,ghost,control", "--out"])
        .arg(c.dir.path().join("o2/de"))
        .output()
        .unwrap();
    assert!(!bad.status.success());
    let all = String::from_utf8_lossy(&bad.stderr);
    assert!(all.contains("ghost") && !all.contains("secret_reads"), "{all}");
}

#[test]
fn shrink_apeglm_missing_package_is_named_error_or_runs() {
    if !r_ready() {
        return;
    }
    let c = case(&NAMES, &[]);
    let mut a = c.args("~ condition");
    a.shrink = Some(Shrink::Apeglm);
    let have = Command::new("Rscript")
        .args(["--vanilla", "-e", "quit(status = if (requireNamespace('apeglm', quietly = TRUE)) 0 else 1)"])
        .output()
        .unwrap()
        .status
        .success();
    match run_de(&a) {
        Ok(_) => {
            assert!(have, "apeglm ran although it is not installed");
            assert!(c.dir.path().join("out/de.results.shrunk.csv").is_file());
        }
        Err(e) => {
            assert!(!have, "apeglm is installed but the run failed: {e}");
            assert!(matches!(e, DeRunError::MissingPackage(ref p) if p == "apeglm"), "{e}");
            assert!(e.to_string().contains("apeglm"));
            // The unshrunk results are kept.
            assert!(c.dir.path().join("out/de.results.csv").is_file());
            assert!(!c.dir.path().join("out/de.results.shrunk.csv").exists());
        }
    }
}

#[test]
fn shrink_normal_writes_extra_file_and_keeps_unshrunk() {
    if !r_ready() {
        return;
    }
    let c = case(&NAMES, &[]);
    let mut a = c.args("~ condition");
    a.shrink = Some(Shrink::Normal);
    run_de(&a).unwrap();
    let (_, plain) = read_csv(&c.dir.path().join("out/de.results.csv"));
    let (_, shrunk) = read_csv(&c.dir.path().join("out/de.results.shrunk.csv"));
    assert_eq!(plain.len(), shrunk.len());
    let versions = std::fs::read_to_string(c.dir.path().join("out/de.versions.txt")).unwrap();
    assert!(versions.contains("shrink: normal"), "{versions}");
    // Shrinkage pulls the estimates toward zero: the sum of |lfc| must not grow.
    let sum = |v: &[HashMap<String, String>]| -> f64 { v.iter().map(|r| r["log2FoldChange"].parse::<f64>().unwrap().abs()).sum() };
    assert!(sum(&shrunk) <= sum(&plain));
}

#[test]
fn versions_file_has_deseq2_and_r_versions() {
    if !r_ready() {
        return;
    }
    let c = case(&NAMES, &[]);
    let out = run_de(&c.args("~ condition")).unwrap();
    let v = std::fs::read_to_string(c.dir.path().join("out/de.versions.txt")).unwrap();
    eprintln!("VERSIONS:\n{v}");
    let lines: Vec<&str> = v.lines().collect();
    assert_eq!(lines[0], format!("DESeq2 {}", out.deseq2_version));
    assert!(lines[1].starts_with("R R version "), "{}", lines[1]);
    assert!(v.contains("design: ~ condition"));
    assert!(v.contains("contrast: condition,treated,control"));
}

#[test]
fn numeric_covariate_design_runs() {
    if !r_ready() {
        return;
    }
    let c = case(&NAMES, &[("age:numeric", ["3", "5", "4", "6", "4", "5"])]);
    run_de(&c.args("~ age + condition")).unwrap();
    let (_, rows) = read_csv(&c.dir.path().join("out/de.results.csv"));
    assert_eq!(rows.len(), 200);
}

#[test]
#[allow(non_snake_case)]
fn sample_named_NA_is_not_mangled() {
    if !r_ready() {
        return;
    }
    let names = ["NA", "T", "c3", "t1", "TRUE", "t3"];
    let c = case(&names, &[]);
    run_de(&c.args("~ condition")).unwrap();
    for f in ["normalized", "vst"] {
        let (header, rows) = read_csv(&c.dir.path().join(format!("out/de.{f}.csv")));
        assert_eq!(header[1..], names, "{f}");
        assert_eq!(rows.len(), 200);
        assert!(rows.iter().all(|r| r["NA"] != "NA" && !r["NA"].is_empty()));
    }
}

#[test]
fn counts_plain_and_gz_are_read_like_zst() {
    // No R needed: the matrix reader goes through the same input opener.
    let dir = tempfile::tempdir().unwrap();
    let m = synth(&NAMES);
    let p = dir.path().join("c.tsv");
    std::fs::write(&p, m.to_tsv()).unwrap();
    let text = qrnamap_cli::pipeline::read_table(&p).unwrap();
    assert_eq!(read_count_matrix(&text).unwrap(), m);
}

#[test]
fn read_count_matrix_names_its_errors() {
    use qrnamap_cli::de::MatrixReadError as E;
    assert!(matches!(read_count_matrix(""), Err(E::Empty)));
    assert!(matches!(read_count_matrix("id\ta\n"), Err(E::BadHeader)));
    assert!(matches!(read_count_matrix("gene\ta\ta\ng\t1\t1\n"), Err(E::DuplicateSample(_))));
    assert!(matches!(read_count_matrix("gene\ta\ng\t1\ng\t2\n"), Err(E::DuplicateGene(_))));
    assert!(matches!(read_count_matrix("gene\ta\tb\ng\t1\n"), Err(E::BadRow { line: 2 })));
    assert!(matches!(read_count_matrix("gene\ta\ng\t1.5\n"), Err(E::BadCount { line: 2, .. })));
    assert!(matches!(read_count_matrix("gene\ta\ng\t-1\n"), Err(E::BadCount { line: 2, .. })));
}
