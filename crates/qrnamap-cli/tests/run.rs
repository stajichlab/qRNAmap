//! `qrnamap run` and the reads-per-sample check.

mod common;

use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use qrnamap_cli::io::open_reader;
use qrnamap_cli::run::{check_reads_per_sample, count_pairs, count_records};
use qrnamap_core::SampleSheet;

fn fastq(n: usize) -> Vec<u8> {
    let mut v = Vec::new();
    for i in 0..n {
        writeln!(v, "@r{i} extra\nACGTACGTAC\n+\nIIIIIIIIII").unwrap();
    }
    v
}

fn gz(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

#[test]
fn count_pairs_plain_and_gzip() {
    let d = tempfile::tempdir().unwrap();
    let plain = d.path().join("a.fq");
    fs::write(&plain, fastq(123)).unwrap();
    // Multi-member gzip: two members of 100 and 23 records.
    let multi = d.path().join("a.fq.gz");
    let mut bytes = gz(&fastq(100));
    bytes.extend(gz(&fastq(23)));
    fs::write(&multi, bytes).unwrap();
    let zst = d.path().join("a.fq.zst");
    fs::write(&zst, zstd::encode_all(&fastq(123)[..], 3).unwrap()).unwrap();
    for p in [&plain, &multi, &zst] {
        assert_eq!(count_pairs(p, None).unwrap(), 123, "{}", p.display());
        assert_eq!(count_pairs(p, Some(1000)).unwrap(), 123);
        assert_eq!(count_pairs(p, Some(123)).unwrap(), 123, "equal to the limit is not over it");
        assert_eq!(count_pairs(p, Some(50)).unwrap(), 51, "stops at limit + 1");
    }
    // A truncated record is an error, not a short count.
    let bad = d.path().join("bad.fq");
    fs::write(&bad, b"@r1\nACGT\n+\n").unwrap();
    assert!(count_pairs(&bad, None).is_err());
    // A truncated gzip stream is a named error.
    let full = gz(&fastq(500));
    let cut = d.path().join("cut.fq.gz");
    fs::write(&cut, &full[..full.len() / 2]).unwrap();
    let e = count_pairs(&cut, None).unwrap_err();
    assert!(e.to_string().contains("truncated"), "{e}");
}

/// A reader that records how many bytes were taken from it.
struct Counting {
    inner: std::io::Cursor<Vec<u8>>,
    taken: Arc<AtomicUsize>,
}

impl Read for Counting {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.taken.fetch_add(n, Ordering::SeqCst);
        Ok(n)
    }
}

#[test]
fn count_pairs_stops_after_limit() {
    let data = fastq(200_000);
    let total = data.len();
    let taken = Arc::new(AtomicUsize::new(0));
    let r = open_reader(Counting { inner: std::io::Cursor::new(data), taken: taken.clone() }, Path::new("mem.fq")).unwrap();
    assert_eq!(count_records(r, Path::new("mem.fq"), Some(10)).unwrap(), 11);
    let used = taken.load(Ordering::SeqCst);
    eprintln!("read {used} of {total} bytes");
    assert!(used < 64 * 1024, "read {used} of {total} bytes; it must stop after the limit");
}

fn sheet_with(d: &Path, rows: &[(&str, usize)]) -> SampleSheet {
    let mut t = String::from("sample\tfastq1\tfastq2\n");
    for (name, n) in rows {
        let r1 = d.join(format!("{name}_1.fq"));
        let r2 = d.join(format!("{name}_2.fq"));
        fs::write(&r1, fastq(*n)).unwrap();
        fs::write(&r2, fastq(*n)).unwrap();
        t.push_str(&format!("{name}\t{}\t{}\n", r1.display(), r2.display()));
    }
    SampleSheet::parse(&t).unwrap()
}

#[test]
fn exceeds_reads_per_sample_gives_warning_not_error() {
    let d = tempfile::tempdir().unwrap();
    let sheet = sheet_with(d.path(), &[("A", 30), ("B", 5), ("C", 6)]);
    let w = check_reads_per_sample(&sheet, 5).unwrap();
    let names: Vec<&str> = w.iter().map(|x| x.sample.as_str()).collect();
    assert_eq!(names, ["A", "C"]);
    let text = w[0].to_string();
    assert!(text.contains("'A'") && text.contains("more than 5") && text.contains("reads_per_sample"), "{text}");
    assert!(check_reads_per_sample(&sheet, 30).unwrap().is_empty());
}

#[test]
fn reads_per_sample_sums_lanes() {
    let d = tempfile::tempdir().unwrap();
    let (a, b) = (d.path().join("l1.fq"), d.path().join("l2.fq"));
    fs::write(&a, fastq(4)).unwrap();
    fs::write(&b, fastq(4)).unwrap();
    let t = format!("sample\tfastq1\nS\t{},{}\n", a.display(), b.display());
    let sheet = SampleSheet::parse(&t).unwrap();
    assert!(check_reads_per_sample(&sheet, 8).unwrap().is_empty());
    assert_eq!(check_reads_per_sample(&sheet, 7).unwrap().len(), 1);
}

// ---------------------------------------------------------------------------------------------
// `qrnamap run` through the binary

use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::Instant;

use common::{records, sim, COUNTS};
use qrnamap_cli::pipeline::read_table;
use qrnamap_core::Strandedness;

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn qrnamap(args: &[&str], path_env: Option<&Path>) -> Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_qrnamap"));
    c.args(args).env_remove("SCRATCH").env_remove("TMPDIR");
    if let Some(p) = path_env {
        c.env("PATH", p);
    }
    c.output().unwrap()
}

fn ok(o: &Output) -> (String, String) {
    let (out, err) = (String::from_utf8_lossy(&o.stdout).into_owned(), String::from_utf8_lossy(&o.stderr).into_owned());
    assert!(o.status.success(), "stdout:\n{out}\nstderr:\n{err}");
    (out, err)
}

fn fail(o: &Output) -> String {
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    assert!(!o.status.success(), "expected failure; stderr:\n{err}");
    err
}

/// Per-transcript scale factors for c1, c2, t1, t2 (variation between replicates and a change
/// of tA1 and tC1 between the groups). Chosen by hand.
const FACTORS: [[f64; 6]; 4] = [
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [1.3, 0.8, 1.1, 0.7, 1.2, 0.9],
    [4.0, 1.0, 0.9, 1.2, 0.25, 1.1],
    [3.2, 1.2, 1.0, 0.8, 0.35, 0.9],
];
const NAMES: [&str; 4] = ["c1", "c2", "t1", "t2"];

/// A 2-condition, 2-replicate sheet of simulated pairs from the mini fixture.
fn mini_sheet(d: &Path) -> PathBuf {
    let recs = records();
    let mut t = String::from("sample\tfastq1\tfastq2\tcondition\n");
    for (i, name) in NAMES.iter().enumerate() {
        let counts: Vec<u64> = COUNTS.iter().zip(FACTORS[i]).map(|(&c, f)| (c as f64 * f).round() as u64).collect();
        let pairs = sim::simulate(&recs, &counts, Strandedness::Unstranded, 100 + i as u64);
        let (r1, r2) = (d.join(format!("{name}_1.fq")), d.join(format!("{name}_2.fq")));
        sim::write_fastq(&pairs, &r1, &r2);
        let cond = if name.starts_with('c') { "control" } else { "treated" };
        t.push_str(&format!("{name}\t{}\t{}\t{cond}\n", r1.display(), r2.display()));
    }
    let p = d.join("sheet.tsv");
    fs::write(&p, t).unwrap();
    p
}

fn mini() -> (PathBuf, PathBuf) {
    (common::mini_genome(), common::mini_gff())
}

/// True when Rscript and DESeq2 are present (fails under QRNAMAP_REQUIRE_TOOLS=1 otherwise).
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

#[test]
fn run_end_to_end_quant_on_mini_fixture() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    let with_de = r_ready();
    let tmp = d.join("tmp");
    let mut args = vec![
        "run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out),
        "--tmpdir", s(&tmp),
    ];
    if with_de {
        args.extend(["--contrast", "condition,treated,control"]);
    }
    let (stdout, stderr) = ok(&qrnamap(&args, None));
    for n in NAMES {
        assert!(stdout.contains(&format!("{n}\treads_processed=")), "{stdout}");
    }
    let m = read_table(&out.join("counts.tsv.zst")).unwrap();
    assert!(m.starts_with("gene\tc1\tc2\tt1\tt2\n"), "{m}");
    // Row set and order from tx2gene: the five fixture genes in transcript order.
    let genes: Vec<&str> = m.lines().skip(1).map(|l| l.split('\t').next().unwrap()).collect();
    assert_eq!(genes, ["geneA", "geneB", "geneR", "geneC", "geneU"], "{m}");
    assert!(out.join(".prepared").is_file());
    assert!(!stderr.contains("read pairs (reads_per_sample"), "no reads warning at the default: {stderr}");
    if with_de {
        for f in ["deseq2.results.csv", "deseq2.normalized.csv", "deseq2.vst.csv", "deseq2.versions.txt"] {
            assert!(out.join(f).is_file(), "{f}\n{stderr}");
        }
        assert!(stdout.contains("DESeq2 "), "{stdout}");
        let res = fs::read_to_string(out.join("deseq2.results.csv")).unwrap();
        assert!(res.contains("geneA"), "{res}");
    }
    eprintln!("WALL run_end_to_end_quant_on_mini_fixture: {:.2?}", t0.elapsed());
}

#[test]
fn run_end_to_end_align_with_task_slices() {
    if !common::require_tools(&["hisat2-build", "hisat2-align-s", "featureCounts"]) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    // Local run.
    let local = d.join("local");
    ok(&qrnamap(
        &["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "align", "--out", s(&local), "--tmpdir", s(&d.join("tmp"))],
        None,
    ));
    // Cluster layout: write the scripts, then run them with bash in order (no SLURM): prepare,
    // each array task (3 samples per task: tasks 0 and 1), final.
    let cl = d.join("cluster");
    let (stdout, _) = ok(&qrnamap(
        &[
            "run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "align", "--out", s(&cl),
            "--slurm", "--samples-per-task", "3", "--qrnamap-bin", env!("CARGO_BIN_EXE_qrnamap"),
        ],
        None,
    ));
    let sdir = cl.join("slurm");
    assert!(stdout.contains(s(&sdir.join("submit.sh"))), "{stdout}");
    let samples = fs::read_to_string(sdir.join("samples.sbatch")).unwrap();
    assert!(samples.contains("#SBATCH --array=0-1\n"), "{samples}");
    let scratch = d.join("node_scratch");
    fs::create_dir_all(&scratch).unwrap();
    let job = |script: &str, task: Option<&str>| {
        let mut c = Command::new("bash");
        c.arg(sdir.join(script)).env("SCRATCH", &scratch);
        if let Some(t) = task {
            c.env("SLURM_ARRAY_TASK_ID", t);
        }
        ok(&c.output().unwrap());
    };
    job("prepare.sbatch", None);
    assert!(cl.join(".prepared").is_file());
    let index_file = cl.join("index").join("genome.1.ht2");
    let built = fs::metadata(&index_file).unwrap().modified().unwrap();
    job("samples.sbatch", Some("0"));
    let made: Vec<bool> = NAMES.iter().map(|n| cl.join(format!("{n}.genecounts.tsv.zst")).is_file()).collect();
    assert_eq!(made, [true, true, true, false], "task 0 holds samples 0..3");
    job("samples.sbatch", Some("1"));
    assert_eq!(fs::metadata(&index_file).unwrap().modified().unwrap(), built, "array tasks reuse the index");
    job("final.sbatch", None);
    let a = read_table(&local.join("counts.tsv.zst")).unwrap();
    let b = read_table(&cl.join("counts.tsv.zst")).unwrap();
    assert_eq!(a, b);
    assert!(a.starts_with("gene\tc1\tc2\tt1\tt2\n"), "{a}");
    // A task id past the last sample is a named error.
    let o = qrnamap(
        &[
            "align", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--out", s(&cl), "--task-id", "2",
            "--samples-per-task", "3", "--tmpdir", s(&scratch),
        ],
        None,
    );
    let err = fail(&o);
    assert!(err.contains("task 2") && err.contains("4 samples"), "{err}");
    eprintln!("WALL run_end_to_end_align_with_task_slices: {:.2?}", t0.elapsed());
}

#[test]
fn task_id_needs_prepared_run_dir() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = d.join("sheet.tsv");
    fs::write(&sh, "sample\tfastq1\tfastq2\nA\t/x/a_1.fq\t/x/a_2.fq\n").unwrap();
    let out = d.join("out");
    let o = qrnamap(
        &["quant", "--samples", s(&sh), "--transcripts", "/x/t.fa", "--tx2gene", "/x/t2g", "--out", s(&out), "--task-id", "0", "--tmpdir", s(d)],
        None,
    );
    let err = fail(&o);
    assert!(err.contains(".prepared") && err.contains("--step prepare"), "{err}");
}

/// Empty directory to use as PATH, so no tool can be found.
fn no_tools(d: &Path) -> PathBuf {
    let p = d.join("empty_path");
    fs::create_dir_all(&p).unwrap();
    p
}

#[test]
fn run_removes_old_outputs_of_all_samples_first() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    // c1 (the first sample) holds reads that match no transcript, so kallisto fails on it.
    let sh = mini_sheet(d);
    let mut rng = sim::Lcg::new(7);
    let junk: Vec<sim::Pair> = (0..200)
        .map(|i| {
            let r: Vec<u8> = (0..75).map(|_| b"ACGT"[rng.below(4)]).collect();
            sim::Pair { name: format!("junk{i}"), r1: r.clone(), r2: sim::revcomp(&r), tx: 0 }
        })
        .collect();
    sim::write_fastq(&junk, &d.join("c1_1.fq"), &d.join("c1_2.fq"));
    let (g, gff) = mini();
    let out = d.join("run");
    fs::create_dir_all(&out).unwrap();
    let stale = out.join("t1.genecounts.tsv.zst");
    qrnamap_cli::pipeline::write_gene_counts(&stale, "t1", &[("geneA".to_string(), 99u64)].into_iter().collect()).unwrap();
    let stale_matrix = out.join("counts.tsv.zst");
    fs::write(&stale_matrix, b"old").unwrap();
    let stale_de = out.join("deseq2.results.csv");
    fs::write(&stale_de, b"old").unwrap();
    let o = qrnamap(
        &["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out), "--tmpdir", s(&d.join("tmp"))],
        None,
    );
    let err = fail(&o);
    assert!(err.contains("sample 'c1'") && err.contains("pseudoaligned 0 reads"), "{err}");
    assert!(!stale.exists(), "the stale count file of sample t1 must be gone");
    assert!(!stale_matrix.exists() && !stale_de.exists());
}

#[test]
fn run_removes_old_outputs_even_when_the_index_cannot_be_built() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    fs::create_dir_all(&out).unwrap();
    let stale = out.join("t2.genecounts.tsv.zst");
    fs::write(&stale, b"old").unwrap();
    fs::write(out.join(".prepared"), b"old").unwrap();
    // reads_per_sample 100 is below every sample: a warning each, then the missing kallisto fails.
    let o = qrnamap(
        &[
            "run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out),
            "--tmpdir", s(&d.join("tmp")), "--reads-per-sample", "100",
        ],
        Some(&no_tools(d)),
    );
    let err = fail(&o);
    assert!(err.contains("kallisto"), "{err}");
    for n in NAMES {
        assert!(err.contains(&format!("warning: sample '{n}' has more than 100 read pairs")), "{err}");
    }
    assert!(!stale.exists());
    assert!(!out.join(".prepared").exists(), "a failed prepare leaves no completion marker");
}

#[test]
fn run_bad_design_fails_before_any_tool() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let empty = no_tools(d);
    for (contrast, want) in [("condition,treated,nosuch", "nosuch"), ("batch,a,b", "batch"), ("condition,treated,treated", "the same")] {
        let out = d.join("run");
        let o = qrnamap(
            &["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out), "--contrast", contrast],
            Some(&empty),
        );
        let err = fail(&o);
        assert!(err.contains(want), "{contrast}: {err}");
        assert!(!err.contains("kallisto"), "{err}");
        assert!(!out.exists(), "no output for a bad design");
    }
    // One replicate in a tested group.
    let one = d.join("one.tsv");
    let text = fs::read_to_string(&sh).unwrap();
    let text: String = text.lines().filter(|l| !l.starts_with("t2\t")).map(|l| format!("{l}\n")).collect();
    fs::write(&one, text).unwrap();
    let out = d.join("run1");
    let o = qrnamap(
        &[
            "run", "--samples", s(&one), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out),
            "--contrast", "condition,treated,control",
        ],
        Some(&empty),
    );
    let err = fail(&o);
    assert!(err.contains("treated"), "{err}");
    assert!(!out.exists());
}

#[test]
fn run_single_end_without_fragment_length_fails_before_any_tool() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = d.join("se.tsv");
    fs::write(&sh, "sample\tfastq1\nSE1\t/x/r.fq\n").unwrap();
    let (g, gff) = mini();
    let out = d.join("run");
    for extra in [&[][..], &["--slurm"][..]] {
        let mut a = vec!["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out)];
        a.extend_from_slice(extra);
        let err = fail(&qrnamap(&a, Some(&no_tools(d))));
        assert!(err.contains("SE1") && err.contains("fragment length"), "{err}");
        assert!(!out.exists());
    }
}
