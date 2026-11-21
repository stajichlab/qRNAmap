//! End-to-end runs of the `qrnamap` binary with the real tools on simulated reads.

mod common;

use std::path::Path;
use std::process::{Command, Output};
use std::time::Instant;

use common::{records, sim, COUNTS};
use qrnamap_cli::pipeline::{read_gene_counts, read_table};
use qrnamap_core::Strandedness;

fn qrnamap(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_qrnamap")).args(args).output().unwrap()
}

fn ok(o: &Output) -> (String, String) {
    let (out, err) = (String::from_utf8_lossy(&o.stdout).into_owned(), String::from_utf8_lossy(&o.stderr).into_owned());
    assert!(o.status.success(), "stdout:\n{out}\nstderr:\n{err}");
    (out, err)
}

fn sheet(d: &Path, lib: Strandedness) -> std::path::PathBuf {
    let pairs = sim::simulate(&records(), &COUNTS, lib, 41);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let p = d.join("sheet.tsv");
    std::fs::write(&p, format!("sample\tfastq1\tfastq2\tcondition\nS1\t{}\t{}\tctl\n", r1.display(), r2.display())).unwrap();
    p
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn cli_quant_align_matrix_end_to_end() {
    if !common::require_tools(&["kallisto", "hisat2-align-s", "hisat2-build", "featureCounts"]) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Unstranded);
    let fx = common::mini_fixture();
    let (genome, gff) = (fx.join("genome.fa"), fx.join("genes.gff3"));
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    ok(&qrnamap(&["transcripts", "--reference", s(&genome), "--gff", s(&gff), "--out", s(&fa), "--tx2gene", s(&t2g)]));

    let q = d.join("q");
    let (out, _) = ok(&qrnamap(&[
        "quant", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&q),
        "--tmpdir", s(&d.join("tmp")), "--threads", "1",
    ]));
    assert!(out.starts_with("S1\treads_processed=10000\tpseudoaligned=100.00%"), "{out}");
    ok(&qrnamap(&["matrix", "--out", s(&q), "--samples", s(&sh)]));

    let a = d.join("a");
    let (out, _) = ok(&qrnamap(&[
        "align", "--samples", s(&sh), "--reference", s(&genome), "--gff", s(&gff), "--out", s(&a),
        "--tmpdir", s(&d.join("tmp")),
    ]));
    assert!(out.starts_with("S1\treads_processed=10000\tmapped=100.00%\tassigned_fraction=1.0000"), "{out}");
    ok(&qrnamap(&["matrix", "--out", s(&a), "--samples", s(&sh)]));

    let qm = read_table(&q.join("counts.tsv.zst")).unwrap();
    let am = read_table(&a.join("counts.tsv.zst")).unwrap();
    assert_eq!(qm, "gene\tS1\ngeneA\t3000\ngeneB\t2500\ngeneC\t1500\ngeneR\t2000\ngeneU\t1000\n");
    assert_eq!(am, qm);
    assert_eq!(read_gene_counts(&a.join("S1.genecounts.tsv.zst"), "S1").unwrap()["geneB"], 2500);
    // The run directory in scratch is removed.
    assert_eq!(std::fs::read_dir(d.join("tmp")).unwrap().count(), 0);
    eprintln!("WALL cli_quant_align_matrix_end_to_end: {:.2?}", t0.elapsed());
}

#[test]
fn cli_strand_check_reports_reverse_and_align_warns_on_wrong_setting() {
    if !common::require_tools(&["hisat2-align-s", "hisat2-build", "featureCounts"]) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Reverse);
    let fx = common::mini_fixture();
    let (genome, gff) = (fx.join("genome.fa"), fx.join("genes.gff3"));
    let sc = d.join("sc");
    let (out, _) = ok(&qrnamap(&[
        "strand-check", "--samples", s(&sh), "--reference", s(&genome), "--gff", s(&gff), "--out", s(&sc),
        "--pairs", "500", "--tmpdir", s(&d.join("tmp")),
    ]));
    assert_eq!(out, "S1\t0\t500\t0\t0.0000\t1.0000\t0.0000\treverse\tfalse\n");
    assert!(read_table(&sc.join("strand_check.tsv")).unwrap().ends_with(&out));

    let (_, err) = ok(&qrnamap(&[
        "align", "--samples", s(&sh), "--reference", s(&genome), "--gff", s(&gff), "--out", s(&d.join("a")),
        "--strandedness", "forward", "--tmpdir", s(&d.join("tmp")),
    ]));
    assert!(err.contains("assigned fraction 0.0000 is below"), "{err}");
    assert!(err.contains("all-zero sample"), "{err}");
    eprintln!("WALL cli_strand_check_reports_reverse_and_align_warns_on_wrong_setting: {:.2?}", t0.elapsed());
}

#[test]
fn cli_names_bad_sample_before_tools_run() {
    // Single-end sample without --fragment-length: named error, no tool runs, no output dir.
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = d.join("sheet.tsv");
    std::fs::write(&sh, "sample\tfastq1\nSE1\t/x/r.fq\n").unwrap();
    let out_dir = d.join("q");
    let o = qrnamap(&["quant", "--samples", s(&sh), "--transcripts", "/x/tx.fa", "--tx2gene", "/x/t2g", "--out", s(&out_dir)]);
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("SE1") && err.contains("fragment length"), "{err}");
    assert!(!out_dir.exists());
}

#[test]
fn cli_tmpdir_order_flag_then_scratch_then_tmpdir() {
    // R23: --tmpdir, then $SCRATCH, then $TMPDIR. A compressed transcript FASTA needs scratch.
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Unstranded);
    let fx = common::mini_fixture();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv"));
    ok(&qrnamap(&["transcripts", "--reference", s(&fx.join("genome.fa")), "--gff", s(&fx.join("genes.gff3")), "--out", s(&fa), "--tx2gene", s(&t2g)]));
    let quant = |out: &str, extra: &[&str], env: &[(&str, &Path)]| {
        let mut c = Command::new(env!("CARGO_BIN_EXE_qrnamap"));
        c.args(["quant", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&d.join(out))])
            .args(extra)
            .env_remove("SCRATCH")
            .env_remove("TMPDIR");
        for (k, v) in env {
            c.env(k, v);
        }
        c.output().unwrap()
    };
    // None of the three: named error.
    let o = quant("q0", &[], &[]);
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--tmpdir") && err.contains("$SCRATCH") && err.contains("$TMPDIR"), "{err}");
    // Only $TMPDIR: used (the run creates the directory).
    let only_tmp = d.join("only_tmpdir");
    ok(&quant("q1", &[], &[("TMPDIR", &only_tmp)]));
    assert!(only_tmp.is_dir());
    // $SCRATCH beats $TMPDIR.
    let (sc, tm) = (d.join("scratch2"), d.join("tmpdir2"));
    ok(&quant("q2", &[], &[("SCRATCH", &sc), ("TMPDIR", &tm)]));
    assert!(sc.is_dir() && !tm.exists());
    // --tmpdir beats both.
    let (flag, sc, tm) = (d.join("flag3"), d.join("scratch3"), d.join("tmpdir3"));
    ok(&quant("q3", &["--tmpdir", s(&flag)], &[("SCRATCH", &sc), ("TMPDIR", &tm)]));
    assert!(flag.is_dir() && !sc.exists() && !tm.exists());
}

#[test]
fn quant_salmon_stdout_says_mapped() {
    if !common::require_tools(&["salmon"]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Unstranded);
    let fx = common::mini_fixture();
    let (genome, gff) = (fx.join("genome.fa"), fx.join("genes.gff3"));
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    ok(&qrnamap(&["transcripts", "--reference", s(&genome), "--gff", s(&gff), "--out", s(&fa), "--tx2gene", s(&t2g)]));
    let q = d.join("q");
    let (out, _) = ok(&qrnamap(&[
        "quant", "--quantifier", "salmon", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g),
        "--out", s(&q), "--tmpdir", s(&d.join("tmp")), "--threads", "1",
    ]));
    assert!(out.starts_with("S1\treads_processed=10000\tmapped="), "{out}");
    assert!(!out.contains("pseudoaligned"), "{out}");
}

#[test]
fn quantifier_with_align_mode_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = d.join("sheet.tsv");
    std::fs::write(&sh, "sample\tfastq1\tfastq2\nA\t/x/a_1.fq\t/x/a_2.fq\n").unwrap();
    let out = d.join("out");
    let base = [
        "run", "--samples", s(&sh), "--reference", "/x/g.fa", "--gff", "/x/g.gff3",
        "--mode", "align", "--out", s(&out), "--tmpdir", s(d), "--step", "prepare",
    ];
    let o = qrnamap(&[&base[..], &["--quantifier", "salmon"]].concat());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success(), "{err}");
    assert!(err.contains("--quantifier") && err.contains("--mode align"), "{err}");
    assert!(!out.exists(), "refused before any file is written");
    // The default quantifier does not trigger the error: the run fails later, on the missing
    // FASTQ files.
    let o = qrnamap(&base);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(!err.contains("--quantifier"), "{err}");
    // `qrnamap align` has no --quantifier flag.
    let o2 = d.join("o2");
    let o = qrnamap(&[
        "align", "--samples", s(&sh), "--reference", "/x/g.fa", "--gff", "/x/g.gff3",
        "--out", s(&o2),
        "--quantifier", "salmon",
    ]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success() && err.contains("--quantifier"), "{err}");
}
