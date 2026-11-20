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
    assert_eq!(read_gene_counts(&a.join("S1.genecounts.tsv.zst")).unwrap()["geneB"], 2500);
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
fn cli_needs_tmpdir_when_decompression_is_needed() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Unstranded);
    let fx = common::mini_fixture();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv"));
    ok(&qrnamap(&["transcripts", "--reference", s(&fx.join("genome.fa")), "--gff", s(&fx.join("genes.gff3")), "--out", s(&fa), "--tx2gene", s(&t2g)]));
    let o = Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args(["quant", "--samples", s(&sh), "--transcripts", s(&fa), "--tx2gene", s(&t2g), "--out", s(&d.join("q"))])
        .env_remove("SCRATCH")
        .output()
        .unwrap();
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--tmpdir") && err.contains("$SCRATCH"), "{err}");
}
