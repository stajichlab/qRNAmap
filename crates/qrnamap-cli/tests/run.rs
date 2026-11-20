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

#[test]
fn count_pairs_ignores_trailing_blank_lines() {
    // R29-12: blank lines after the last record (plain and gzip) are not an error.
    let d = tempfile::tempdir().unwrap();
    let mut data = fastq(7);
    data.extend_from_slice(b"\n\n");
    let plain = d.path().join("t.fq");
    fs::write(&plain, &data).unwrap();
    let gzp = d.path().join("t.fq.gz");
    fs::write(&gzp, gz(&data)).unwrap();
    for p in [&plain, &gzp] {
        assert_eq!(count_pairs(p, None).unwrap(), 7, "{}", p.display());
    }
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
    // The marker written by the local run matches a task run with the same inputs (R27).
    let (fa, t2g) = (out.join("transcripts.fa.zst"), out.join("tx2gene.tsv.zst"));
    ok(&qrnamap(
        &[
            "quant", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&out), "--reference", s(&g),
            "--gff", s(&gff), "--task-id", "1", "--samples-per-task", "2", "--tmpdir", s(&tmp),
        ],
        None,
    ));
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

/// The reviewer's two C1 inputs for one sample S1 of 10,000 simulated pairs: (a) R2 `.gz` cut
/// at half its bytes; (b) plain R1 of 10,000 records with R2 of 2,000 records. Returns the two
/// sheet paths.
fn c1_sheets(d: &Path) -> (PathBuf, PathBuf) {
    let pairs = sim::simulate(&records(), &COUNTS, Strandedness::Unstranded, 41);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let (g1, g2) = (d.join("s_1.fq.gz"), d.join("s_2.fq.gz"));
    fs::write(&g1, gz(&fs::read(&r1).unwrap())).unwrap();
    let full = gz(&fs::read(&r2).unwrap());
    fs::write(&g2, &full[..full.len() / 2]).unwrap();
    let short: String = fs::read_to_string(&r2).unwrap().lines().take(8000).map(|l| format!("{l}\n")).collect();
    let r2s = d.join("short_2.fq");
    fs::write(&r2s, short).unwrap();
    let (a, b) = (d.join("trunc.tsv"), d.join("short.tsv"));
    fs::write(&a, format!("sample\tfastq1\tfastq2\nS1\t{}\t{}\n", g1.display(), g2.display())).unwrap();
    fs::write(&b, format!("sample\tfastq1\tfastq2\nS1\t{}\t{}\n", r1.display(), r2s.display())).unwrap();
    (a, b)
}

#[test]
fn quant_and_run_quant_refuse_truncated_gz_and_short_mate() {
    // C1: both used to exit 0 with partial counts.
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (trunc, short) = c1_sheets(d);
    let (g, gff) = mini();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    let tmp = d.join("tmp");
    for (sheet, want) in [(&trunc, "truncated"), (&short, "10000")] {
        let out = d.join("q");
        let err = fail(&qrnamap(
            &["quant", "--samples", s(sheet), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&out), "--tmpdir", s(&tmp)],
            None,
        ));
        assert!(err.contains("sample 'S1'") && err.contains(want), "quant {}: {err}", sheet.display());
        assert!(!out.join("S1.genecounts.tsv.zst").exists());
        let out = d.join("r");
        let err = fail(&qrnamap(
            &["run", "--samples", s(sheet), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out), "--tmpdir", s(&tmp)],
            None,
        ));
        assert!(err.contains("sample 'S1'") && err.contains(want), "run {}: {err}", sheet.display());
        assert!(!out.join("S1.genecounts.tsv.zst").exists());
        assert!(!out.join("counts.tsv.zst").exists());
    }
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
    let (g, gff) = mini();
    // R27: with --task-id, quant needs --reference and --gff (checked against the marker).
    let o = qrnamap(
        &["quant", "--samples", s(&sh), "--transcripts", "/x/t.fa", "--tx2gene", "/x/t2g", "--out", s(&out), "--task-id", "0", "--tmpdir", s(d)],
        None,
    );
    let err = fail(&o);
    assert!(err.contains("--task-id") && err.contains("--reference") && err.contains("--gff"), "{err}");
    let o = qrnamap(
        &[
            "quant", "--samples", s(&sh), "--transcripts", "/x/t.fa", "--tx2gene", "/x/t2g", "--out", s(&out), "--task-id", "0",
            "--tmpdir", s(d), "--reference", s(&g), "--gff", s(&gff),
        ],
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

// ---------------------------------------------------------------------------------------------
// Ruling R27: the prepared marker records the inputs; a stale prepared state is refused.

use qrnamap_cli::run::{check_prepared, write_prepared_marker, RunError, PREPARED_MARKER};
use qrnamap_core::Mode;

/// Copies of the mini reference and GFF (so tests can change them) and an empty run dir.
struct Prepared {
    _dir: tempfile::TempDir,
    d: PathBuf,
    reference: PathBuf,
    gff: PathBuf,
    out: PathBuf,
    sheet: PathBuf,
}

fn prepared_align_dir() -> Prepared {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path().to_path_buf();
    let (g, gff) = mini();
    let (reference, gff_copy) = (d.join("genome.fa"), d.join("genes.gff3"));
    fs::copy(g, &reference).unwrap();
    fs::copy(gff, &gff_copy).unwrap();
    let out = d.join("run");
    fs::create_dir_all(&out).unwrap();
    write_prepared_marker(&out, Mode::Align, &[("reference", &reference), ("gff", &gff_copy)]).unwrap();
    let sheet = mini_sheet(&d);
    Prepared { _dir: dir, d, reference, gff: gff_copy, out, sheet }
}

fn not_prepared_reason(r: Result<(), RunError>) -> String {
    match r {
        Err(RunError::NotPrepared { reason, .. }) => reason,
        other => panic!("expected NotPrepared, got {other:?}"),
    }
}

fn align_task(p: &Prepared, gff: &Path) -> Output {
    qrnamap(
        &[
            "align", "--samples", s(&p.sheet), "--reference", s(&p.reference), "--gff", s(gff), "--out", s(&p.out),
            "--task-id", "0", "--tmpdir", s(&p.d),
        ],
        Some(&no_tools(&p.d)),
    )
}

#[test]
fn marker_records_inputs_and_version() {
    let p = prepared_align_dir();
    let text = fs::read_to_string(p.out.join(PREPARED_MARKER)).unwrap();
    eprintln!("{text}");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "qrnamap-prepared\t1");
    assert_eq!(lines[1], format!("version\t{}", env!("CARGO_PKG_VERSION")));
    assert_eq!(lines[2], "mode\talign");
    for (i, (role, path)) in [("reference", &p.reference), ("gff", &p.gff)].iter().enumerate() {
        let m = fs::metadata(path).unwrap();
        let t = m.modified().unwrap().duration_since(std::time::UNIX_EPOCH).unwrap();
        let canon = fs::canonicalize(path).unwrap();
        assert_eq!(
            lines[3 + i],
            format!("input\t{role}\t{}\t{}\t{}\t{}", m.len(), t.as_secs(), t.subsec_nanos(), canon.display())
        );
    }
    assert_eq!(lines.len(), 5);
    assert!(check_prepared(&p.out, Mode::Align, &[("reference", &p.reference), ("gff", &p.gff)]).is_ok());
    // Mode and version are checked too.
    let r = not_prepared_reason(check_prepared(&p.out, Mode::Quant, &[("reference", &p.reference), ("gff", &p.gff)]));
    assert!(r.contains("mode") && r.contains("align") && r.contains("quant"), "{r}");
    let text2 = text.replace(&format!("version\t{}", env!("CARGO_PKG_VERSION")), "version\t0.0.0-old");
    fs::write(p.out.join(PREPARED_MARKER), text2).unwrap();
    let r = not_prepared_reason(check_prepared(&p.out, Mode::Align, &[("reference", &p.reference), ("gff", &p.gff)]));
    assert!(r.contains("version") && r.contains("0.0.0-old"), "{r}");
}

#[test]
fn task_run_with_different_gff_is_not_prepared() {
    let p = prepared_align_dir();
    let other = p.d.join("other.gff3");
    fs::copy(&p.gff, &other).unwrap();
    let r = not_prepared_reason(check_prepared(&p.out, Mode::Align, &[("reference", &p.reference), ("gff", &other)]));
    assert!(r.starts_with("gff: path differs"), "{r}");
    assert!(r.contains("genes.gff3") && r.contains("other.gff3"), "{r}");
    let err = fail(&align_task(&p, &other));
    assert!(err.contains("not prepared") && err.contains("gff: path differs"), "{err}");
    assert!(!err.contains("hisat2"), "refused before any tool: {err}");
}

#[test]
fn task_run_with_modified_reference_mtime_is_not_prepared() {
    let p = prepared_align_dir();
    let f = fs::OpenOptions::new().write(true).open(&p.reference).unwrap();
    f.set_modified(std::time::SystemTime::now() + std::time::Duration::from_secs(3600)).unwrap();
    drop(f);
    let r = not_prepared_reason(check_prepared(&p.out, Mode::Align, &[("reference", &p.reference), ("gff", &p.gff)]));
    assert!(r.starts_with("reference: mtime changed"), "{r}");
    let err = fail(&align_task(&p, &p.gff));
    assert!(err.contains("reference: mtime changed"), "{err}");
}

#[test]
fn task_run_with_modified_reference_size_is_not_prepared() {
    let p = prepared_align_dir();
    let before = fs::metadata(&p.reference).unwrap().modified().unwrap();
    let mut f = fs::OpenOptions::new().append(true).open(&p.reference).unwrap();
    f.write_all(b"ACGT\n").unwrap();
    f.set_modified(before).unwrap();
    drop(f);
    let r = not_prepared_reason(check_prepared(&p.out, Mode::Align, &[("reference", &p.reference), ("gff", &p.gff)]));
    assert!(r.starts_with("reference: size changed"), "{r}");
    let err = fail(&align_task(&p, &p.gff));
    assert!(err.contains("reference: size changed"), "{err}");
}

fn final_step(p: &Prepared, gff: &Path) -> Output {
    qrnamap(
        &[
            "run", "--samples", s(&p.sheet), "--reference", s(&p.reference), "--gff", s(gff), "--mode", "align", "--out",
            s(&p.out), "--step", "final", "--tmpdir", s(&p.d),
        ],
        Some(&no_tools(&p.d)),
    )
}

#[test]
fn final_step_with_different_gff_is_not_prepared() {
    let p = prepared_align_dir();
    let other = p.d.join("other.gff3");
    fs::copy(&p.gff, &other).unwrap();
    let err = fail(&final_step(&p, &other));
    assert!(err.contains("not prepared") && err.contains("gff: path differs"), "{err}");
    assert!(!p.out.join("counts.tsv.zst").exists());
}

#[test]
fn final_step_without_marker_is_not_prepared() {
    let p = prepared_align_dir();
    fs::remove_file(p.out.join(PREPARED_MARKER)).unwrap();
    let r = not_prepared_reason(check_prepared(&p.out, Mode::Align, &[("reference", &p.reference), ("gff", &p.gff)]));
    assert!(r.contains("is missing"), "{r}");
    let err = fail(&final_step(&p, &p.gff));
    assert!(err.contains("not prepared") && err.contains("is missing"), "{err}");
}

#[test]
fn plain_quant_removes_stale_marker() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    fs::create_dir_all(&out).unwrap();
    let (fa, t2g) = (out.join("transcripts.fa.zst"), out.join("tx2gene.tsv.zst"));
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    write_prepared_marker(&out, Mode::Quant, &[("reference", &g), ("gff", &gff), ("transcripts", &fa), ("tx2gene", &t2g)]).unwrap();
    // No kallisto on PATH: the index build fails, but the old marker must already be gone.
    let o = qrnamap(
        &["quant", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&out), "--tmpdir", s(d)],
        Some(&no_tools(d)),
    );
    let err = fail(&o);
    assert!(err.contains("kallisto"), "{err}");
    assert!(!out.join(PREPARED_MARKER).exists(), "plain quant must remove the stale marker");
    // The same for align.
    write_prepared_marker(&out, Mode::Align, &[("reference", &g), ("gff", &gff)]).unwrap();
    let o = qrnamap(
        &["align", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--out", s(&out), "--tmpdir", s(d)],
        Some(&no_tools(d)),
    );
    fail(&o);
    assert!(!out.join(PREPARED_MARKER).exists(), "plain align must remove the stale marker");
}

// ---------------------------------------------------------------------------------------------
// Final review I1, I2, I3: stale state.

#[test]
fn strand_check_into_prepared_dir_removes_marker() {
    // I1: strand-check rebuilds the index, genes.saf and splice_sites.txt in --out from its own
    // GFF; the old marker must not survive to vouch for them.
    if !common::require_tools(&["hisat2-build", "hisat2-align-s", "featureCounts"]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let gff_b = d.join("genes_b.gff3");
    fs::write(&gff_b, fs::read_to_string(&gff).unwrap().replace("Name=geneA", "Name=geneA_renamed")).unwrap();
    let out = d.join("run");
    let tmp = d.join("tmp");
    ok(&qrnamap(
        &["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "align", "--out", s(&out), "--tmpdir", s(&tmp), "--step", "prepare"],
        None,
    ));
    assert!(out.join(PREPARED_MARKER).is_file());
    ok(&qrnamap(
        &["strand-check", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff_b), "--out", s(&out), "--pairs", "500", "--tmpdir", s(&tmp)],
        None,
    ));
    assert!(!out.join(PREPARED_MARKER).exists(), "strand-check must remove the marker");
    let err = fail(&qrnamap(
        &["align", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--out", s(&out), "--task-id", "0", "--tmpdir", s(&tmp)],
        None,
    ));
    assert!(err.contains("not prepared"), "{err}");
    let err = fail(&qrnamap(
        &["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "align", "--out", s(&out), "--tmpdir", s(&tmp), "--step", "final"],
        None,
    ));
    assert!(err.contains("not prepared"), "{err}");
}

/// A sheet S1, S2, S3 where S1's mate 2 is a gzip file cut in half (it fails), and S2, S3 are
/// valid plain pairs.
fn sheet_first_sample_fails(d: &Path) -> PathBuf {
    let recs = records();
    let mut t = String::from("sample\tfastq1\tfastq2\n");
    for (i, name) in ["S1", "S2", "S3"].iter().enumerate() {
        let pairs = sim::simulate(&recs, &COUNTS, Strandedness::Unstranded, 200 + i as u64);
        let (r1, mut r2) = (d.join(format!("{name}_1.fq")), d.join(format!("{name}_2.fq")));
        sim::write_fastq(&pairs, &r1, &r2);
        if i == 0 {
            let full = gz(&fs::read(&r2).unwrap());
            r2 = d.join("S1_2.fq.gz");
            fs::write(&r2, &full[..full.len() / 2]).unwrap();
        }
        t.push_str(&format!("{name}\t{}\t{}\n", r1.display(), r2.display()));
    }
    let p = d.join("sheet3.tsv");
    fs::write(&p, t).unwrap();
    p
}

#[test]
fn plain_quant_and_align_remove_old_outputs_of_all_samples_first() {
    // I2: the first sample fails; a later sample must not keep a count file of an earlier run.
    let have_quant = common::require_tool("kallisto").is_some();
    let have_align = common::require_tools(&["hisat2-build", "hisat2-align-s", "featureCounts"]);
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet_first_sample_fails(d);
    let (g, gff) = mini();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    let tmp = d.join("tmp");
    let cases: Vec<(&str, Vec<&str>)> = vec![
        ("quant", vec!["--transcripts", s(&fa), "--tx2gene", s(&t2g)]),
        ("align", vec!["--reference", s(&g), "--gff", s(&gff)]),
    ];
    for (cmd, extra) in cases {
        if (cmd == "quant" && !have_quant) || (cmd == "align" && !have_align) {
            continue;
        }
        let out = d.join(format!("out_{cmd}"));
        fs::create_dir_all(&out).unwrap();
        let stale = out.join("S3.genecounts.tsv.zst");
        qrnamap_cli::pipeline::write_gene_counts(&stale, "S3", &[("geneA".to_string(), 99u64)].into_iter().collect()).unwrap();
        let mut a = vec![cmd, "--samples", s(&sh), "--out", s(&out), "--tmpdir", s(&tmp)];
        a.extend(extra);
        let err = fail(&qrnamap(&a, None));
        assert!(err.contains("sample 'S1'") && err.contains("truncated"), "{cmd}: {err}");
        assert!(!stale.exists(), "{cmd}: the stale count file of S3 must be gone");
    }
}

#[test]
fn missing_or_unreadable_fastq_is_named_before_any_tool() {
    // R29-11: quant, align and run check every sheet FASTQ before any index build. No tool is on
    // PATH, so an error that names a tool means the check came too late.
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    let empty = no_tools(d);
    let missing = d.join("t1_2.fq");
    let unreadable = d.join("c2_1.fq");
    for (case, file, sample) in [("missing", &missing, "t1"), ("unreadable", &unreadable, "c2")] {
        if case == "missing" {
            fs::remove_file(&missing).unwrap();
        } else {
            fs::set_permissions(&unreadable, fs::Permissions::from_mode(0o000)).unwrap();
            if fs::File::open(&unreadable).is_ok() {
                eprintln!("SKIP unreadable case: this user can read a mode-000 file (root?)");
                continue;
            }
        }
        let out = d.join(format!("out_{case}"));
        let runs: [Vec<&str>; 3] = [
            vec!["quant", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&out), "--tmpdir", s(d)],
            vec!["align", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--out", s(&out), "--tmpdir", s(d)],
            vec!["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "align", "--out", s(&out), "--tmpdir", s(d)],
        ];
        for a in runs {
            let err = fail(&qrnamap(&a, Some(&empty)));
            assert!(err.contains(&format!("sample '{sample}'")) && err.contains(s(file)), "{case} {}: {err}", a[0]);
            for tool in ["kallisto", "hisat2", "featureCounts"] {
                assert!(!err.contains(tool), "{case} {}: checked after a tool start: {err}", a[0]);
            }
            assert!(!out.join("index").exists() && !out.join("transcripts.kidx").exists());
        }
    }
    fs::set_permissions(&unreadable, fs::Permissions::from_mode(0o644)).unwrap();
}

#[test]
fn out_under_scratch_gives_shared_storage_warning() {
    // R29-13: --out on node-local scratch is a warning (index files must be on shared storage).
    use qrnamap_cli::run::scratch_out_warning;
    let w = scratch_out_warning(Path::new("/scratch/u/123/run"), None).expect("a /scratch/ path warns");
    assert!(w.contains("/scratch/u/123/run") && w.contains("shared storage"), "{w}");
    assert!(scratch_out_warning(Path::new("/node/tmp/x/run"), Some(std::ffi::OsStr::new("/node/tmp/x"))).is_some());
    assert!(scratch_out_warning(Path::new("/node/tmp/x"), Some(std::ffi::OsStr::new("/node/tmp/x"))).is_some());
    assert!(scratch_out_warning(Path::new("/bigdata/lab/run"), Some(std::ffi::OsStr::new("/node/tmp/x"))).is_none());
    assert!(scratch_out_warning(Path::new("/bigdata/lab/run"), None).is_none());
    assert!(scratch_out_warning(Path::new("/scratchy/run"), None).is_none(), "prefix must be the directory /scratch/");
    assert!(scratch_out_warning(Path::new("/node/tmp/xy"), Some(std::ffi::OsStr::new("/node/tmp/x"))).is_none());
    // Through the binary: SCRATCH set, --out under it. It is a warning, not an error: the run
    // goes on and here fails later for a reason that has nothing to do with the warning.
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = mini_sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    let o = Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args(["run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out)])
        .args(["--contrast", "condition,treated,nosuch"])
        .env("SCRATCH", d)
        .env("PATH", no_tools(d))
        .output()
        .unwrap();
    let err = fail(&o);
    assert!(err.contains("warning:") && err.contains("shared storage"), "{err}");
    assert!(err.contains("nosuch"), "the run continued past the warning: {err}");
}

/// Write a gene count file and a summary file for every sheet sample of `p`.
fn write_counts_for_sheet(p: &Prepared) {
    for n in NAMES {
        let m = [("geneA".to_string(), 5u64), ("geneB".to_string(), 7u64)].into_iter().collect();
        qrnamap_cli::pipeline::write_gene_counts(&p.out.join(format!("{n}.genecounts.tsv.zst")), n, &m).unwrap();
        let summary = qrnamap_cli::pipeline::SampleSummary { reads_processed: 1, fraction_mapped: 1.0, fraction_assigned: 1.0 };
        qrnamap_cli::pipeline::write_sample_summary(&p.out, n, qrnamap_core::Mode::Quant, &summary, Strandedness::Unstranded).unwrap();
    }
}

fn set_mtime(path: &Path, t: std::time::SystemTime) {
    let f = fs::OpenOptions::new().write(true).open(path).unwrap();
    f.set_modified(t).unwrap();
}

#[test]
fn final_step_refuses_count_files_older_than_marker() {
    // I3: count files from before the current prepare must not reach the matrix.
    let p = prepared_align_dir();
    write_counts_for_sheet(&p);
    let now = std::time::SystemTime::now();
    set_mtime(&p.out.join(PREPARED_MARKER), now + std::time::Duration::from_secs(3600));
    let err = fail(&final_step(&p, &p.gff));
    assert!(err.contains("sample 'c1'") && err.contains("older than"), "{err}");
    assert!(!p.out.join("counts.tsv.zst").exists());
    // Count files newer than the marker pass.
    set_mtime(&p.out.join(PREPARED_MARKER), now - std::time::Duration::from_secs(3600));
    ok(&final_step(&p, &p.gff));
    assert!(p.out.join("counts.tsv.zst").is_file());
}
