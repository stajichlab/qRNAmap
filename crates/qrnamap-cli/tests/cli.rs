//! End-to-end runs of the `qrnamap` binary with the real tools on simulated reads.

mod common;

use std::path::Path;
use std::process::{Command, Output};
use std::time::Instant;

use common::{records, sim, COUNTS};
use qrnamap_cli::pipeline::{read_gene_counts, read_table};
use qrnamap_core::Strandedness;

fn qrnamap(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args(args)
        .output()
        .unwrap()
}

fn ok(o: &Output) -> (String, String) {
    let (out, err) = (
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    );
    assert!(o.status.success(), "stdout:\n{out}\nstderr:\n{err}");
    (out, err)
}

fn sheet(d: &Path, lib: Strandedness) -> std::path::PathBuf {
    let pairs = sim::simulate(&records(), &COUNTS, lib, 41);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let p = d.join("sheet.tsv");
    std::fs::write(
        &p,
        format!(
            "sample\tfastq1\tfastq2\tcondition\nS1\t{}\t{}\tctl\n",
            r1.display(),
            r2.display()
        ),
    )
    .unwrap();
    p
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn cli_quant_align_matrix_end_to_end() {
    if !common::require_tools(&[
        "kallisto",
        "hisat2-align-s",
        "hisat2-build",
        "featureCounts",
    ]) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Unstranded);
    let fx = common::mini_fixture();
    let (genome, gff) = (fx.join("genome.fa"), fx.join("genes.gff3"));
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    ok(&qrnamap(&[
        "transcripts",
        "--reference",
        s(&genome),
        "--gff",
        s(&gff),
        "--out",
        s(&fa),
        "--tx2gene",
        s(&t2g),
    ]));

    let q = d.join("q");
    let (out, _) = ok(&qrnamap(&[
        "quant",
        "--samples",
        s(&sh),
        "--transcripts",
        s(&fa),
        "--tx2gene",
        s(&t2g),
        "--out",
        s(&q),
        "--tmpdir",
        s(&d.join("tmp")),
        "--threads",
        "1",
    ]));
    assert!(
        out.starts_with("S1\treads_processed=10000\tpseudoaligned=100.00%"),
        "{out}"
    );
    ok(&qrnamap(&["matrix", "--out", s(&q), "--samples", s(&sh)]));

    let a = d.join("a");
    let (out, _) = ok(&qrnamap(&[
        "align",
        "--samples",
        s(&sh),
        "--reference",
        s(&genome),
        "--gff",
        s(&gff),
        "--out",
        s(&a),
        "--tmpdir",
        s(&d.join("tmp")),
    ]));
    assert!(
        out.starts_with("S1\treads_processed=10000\tmapped=100.00%\tassigned_fraction=1.0000"),
        "{out}"
    );
    ok(&qrnamap(&["matrix", "--out", s(&a), "--samples", s(&sh)]));

    let qm = read_table(&q.join("counts.tsv.zst")).unwrap();
    let am = read_table(&a.join("counts.tsv.zst")).unwrap();
    assert_eq!(
        qm,
        "gene\tS1\ngeneA\t3000\ngeneB\t2500\ngeneC\t1500\ngeneR\t2000\ngeneU\t1000\n"
    );
    assert_eq!(am, qm);
    assert_eq!(
        read_gene_counts(&a.join("S1.genecounts.tsv.zst"), "S1").unwrap()["geneB"],
        2500
    );
    // The run directory in scratch is removed.
    assert_eq!(std::fs::read_dir(d.join("tmp")).unwrap().count(), 0);
    eprintln!(
        "WALL cli_quant_align_matrix_end_to_end: {:.2?}",
        t0.elapsed()
    );
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
        "strand-check",
        "--samples",
        s(&sh),
        "--reference",
        s(&genome),
        "--gff",
        s(&gff),
        "--out",
        s(&sc),
        "--pairs",
        "500",
        "--tmpdir",
        s(&d.join("tmp")),
    ]));
    assert_eq!(
        out,
        "S1\t0\t500\t0\t0.0000\t1.0000\t0.0000\treverse\tfalse\n"
    );
    assert!(read_table(&sc.join("strand_check.tsv"))
        .unwrap()
        .ends_with(&out));

    let (_, err) = ok(&qrnamap(&[
        "align",
        "--samples",
        s(&sh),
        "--reference",
        s(&genome),
        "--gff",
        s(&gff),
        "--out",
        s(&d.join("a")),
        "--strandedness",
        "forward",
        "--tmpdir",
        s(&d.join("tmp")),
    ]));
    assert!(err.contains("assigned fraction 0.0000 is below"), "{err}");
    assert!(err.contains("all-zero sample"), "{err}");
    eprintln!(
        "WALL cli_strand_check_reports_reverse_and_align_warns_on_wrong_setting: {:.2?}",
        t0.elapsed()
    );
}

#[test]
fn cli_names_bad_sample_before_tools_run() {
    // Single-end sample without --fragment-length: named error, no tool runs, no output dir.
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = d.join("sheet.tsv");
    std::fs::write(&sh, "sample\tfastq1\nSE1\t/x/r.fq\n").unwrap();
    let out_dir = d.join("q");
    let o = qrnamap(&[
        "quant",
        "--samples",
        s(&sh),
        "--transcripts",
        "/x/tx.fa",
        "--tx2gene",
        "/x/t2g",
        "--out",
        s(&out_dir),
    ]);
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("SE1") && err.contains("fragment length"),
        "{err}"
    );
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
    ok(&qrnamap(&[
        "transcripts",
        "--reference",
        s(&fx.join("genome.fa")),
        "--gff",
        s(&fx.join("genes.gff3")),
        "--out",
        s(&fa),
        "--tx2gene",
        s(&t2g),
    ]));
    let quant = |out: &str, extra: &[&str], env: &[(&str, &Path)]| {
        let mut c = Command::new(env!("CARGO_BIN_EXE_qrnamap"));
        c.args([
            "quant",
            "--samples",
            s(&sh),
            "--transcripts",
            s(&fa),
            "--tx2gene",
            s(&t2g),
            "--out",
            s(&d.join(out)),
        ])
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
    assert!(
        err.contains("--tmpdir") && err.contains("$SCRATCH") && err.contains("$TMPDIR"),
        "{err}"
    );
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
    ok(&quant(
        "q3",
        &["--tmpdir", s(&flag)],
        &[("SCRATCH", &sc), ("TMPDIR", &tm)],
    ));
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
    ok(&qrnamap(&[
        "transcripts",
        "--reference",
        s(&genome),
        "--gff",
        s(&gff),
        "--out",
        s(&fa),
        "--tx2gene",
        s(&t2g),
    ]));
    let q = d.join("q");
    let (out, _) = ok(&qrnamap(&[
        "quant",
        "--quantifier",
        "salmon",
        "--samples",
        s(&sh),
        "--transcripts",
        s(&fa),
        "--tx2gene",
        s(&t2g),
        "--out",
        s(&q),
        "--tmpdir",
        s(&d.join("tmp")),
        "--threads",
        "1",
    ]));
    assert!(
        out.starts_with("S1\treads_processed=10000\tmapped="),
        "{out}"
    );
    assert!(!out.contains("pseudoaligned"), "{out}");
}

#[test]
fn quant_salmon_with_overhang_flags_runs() {
    if !common::require_tools(&["salmon"]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d, Strandedness::Unstranded);
    let fx = common::mini_fixture();
    let (genome, gff) = (fx.join("genome.fa"), fx.join("genes.gff3"));
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    ok(&qrnamap(&[
        "transcripts",
        "--reference",
        s(&genome),
        "--gff",
        s(&gff),
        "--out",
        s(&fa),
        "--tx2gene",
        s(&t2g),
    ]));
    let q = d.join("q");
    let (out, _) = ok(&qrnamap(&[
        "quant",
        "--quantifier",
        "salmon",
        "--samples",
        s(&sh),
        "--transcripts",
        s(&fa),
        "--tx2gene",
        s(&t2g),
        "--out",
        s(&q),
        "--tmpdir",
        s(&d.join("tmp")),
        "--threads",
        "1",
        "--salmon-softclip-overhangs",
        "--salmon-min-score-fraction",
        "0.5",
    ]));
    assert!(
        out.starts_with("S1\treads_processed=10000\tmapped="),
        "{out}"
    );
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
        "run",
        "--samples",
        s(&sh),
        "--reference",
        "/x/g.fa",
        "--gff",
        "/x/g.gff3",
        "--mode",
        "align",
        "--out",
        s(&out),
        "--tmpdir",
        s(d),
        "--step",
        "prepare",
    ];
    let o = qrnamap(&[&base[..], &["--quantifier", "salmon"]].concat());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success(), "{err}");
    assert!(
        err.contains("--quantifier") && err.contains("--mode align"),
        "{err}"
    );
    assert!(!out.exists(), "refused before any file is written");
    // The default quantifier does not trigger the error: the run fails later, on the missing
    // FASTQ files.
    let o = qrnamap(&base);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(!err.contains("--quantifier"), "{err}");
    assert!(
        err.contains("error: sample 'A': input error: /x/a_1.fq"),
        "{err}"
    );
    // `qrnamap align` has no --quantifier flag.
    let o2 = d.join("o2");
    let o = qrnamap(&[
        "align",
        "--samples",
        s(&sh),
        "--reference",
        "/x/g.fa",
        "--gff",
        "/x/g.gff3",
        "--out",
        s(&o2),
        "--quantifier",
        "salmon",
    ]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success() && err.contains("--quantifier"), "{err}");
}

// ---------------------------------------------------------------------------------------------
// Salmon overhang flags

/// Arguments of `quant` for a sheet and an output directory that need no tool: the checks under
/// test run before any file is read or written.
fn quant_args<'a>(sh: &'a str, out: &'a str) -> Vec<&'a str> {
    vec![
        "quant",
        "--samples",
        sh,
        "--transcripts",
        "/x/tx.fa.zst",
        "--tx2gene",
        "/x/t2g.tsv.zst",
        "--out",
        out,
    ]
}

fn write_small_sheet(d: &Path) -> std::path::PathBuf {
    let sh = d.join("sheet.tsv");
    std::fs::write(&sh, "sample\tfastq1\tfastq2\nA\t/x/a_1.fq\t/x/a_2.fq\n").unwrap();
    sh
}

#[test]
fn salmon_flags_rejected_with_kallisto() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = write_small_sheet(d);
    let out = d.join("out");
    let base = quant_args(s(&sh), s(&out));
    for flag in [
        &["--salmon-softclip-overhangs"][..],
        &["--salmon-min-score-fraction", "0.5"][..],
    ] {
        let o = qrnamap(&[&base[..], &["--quantifier", "kallisto"], flag].concat());
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{flag:?}: {err}");
        assert!(
            err.contains(flag[0]) && err.contains("--quantifier salmon"),
            "{flag:?}: {err}"
        );
        assert!(!out.exists(), "refused before any file is written");
    }
    // `run` as well.
    let o = qrnamap(&[
        "run",
        "--samples",
        s(&sh),
        "--reference",
        "/x/g.fa",
        "--gff",
        "/x/g.gff3",
        "--mode",
        "quant",
        "--quantifier",
        "kallisto",
        "--salmon-softclip-overhangs",
        "--out",
        s(&out),
    ]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success(), "{err}");
    assert!(
        err.contains("--salmon-softclip-overhangs") && err.contains("--quantifier salmon"),
        "{err}"
    );
    assert!(!out.exists());
}

#[test]
fn salmon_flags_rejected_with_default_quantifier() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = write_small_sheet(d);
    let out = d.join("out");
    let base = quant_args(s(&sh), s(&out));
    let o = qrnamap(&[&base[..], &["--salmon-min-score-fraction", "0.5"]].concat());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success(), "{err}");
    assert!(
        err.contains("--salmon-min-score-fraction") && err.contains("--quantifier salmon"),
        "{err}"
    );
    assert!(!out.exists());
    // Both flags: both are named.
    let o = qrnamap(
        &[
            &base[..],
            &[
                "--salmon-softclip-overhangs",
                "--salmon-min-score-fraction",
                "0.5",
            ],
        ]
        .concat(),
    );
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(
        err.contains("--salmon-softclip-overhangs") && err.contains("--salmon-min-score-fraction"),
        "{err}"
    );
}

#[test]
fn salmon_flags_rejected_with_align_mode() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = write_small_sheet(d);
    let out = d.join("out");
    let base = [
        "run",
        "--samples",
        s(&sh),
        "--reference",
        "/x/g.fa",
        "--gff",
        "/x/g.gff3",
        "--mode",
        "align",
        "--out",
        s(&out),
        "--tmpdir",
        s(d),
        "--step",
        "prepare",
    ];
    for flag in [
        &["--salmon-softclip-overhangs"][..],
        &["--salmon-min-score-fraction", "0.5"][..],
    ] {
        for q in [&[][..], &["--quantifier", "salmon"][..]] {
            let o = qrnamap(&[&base[..], flag, q].concat());
            let err = String::from_utf8_lossy(&o.stderr);
            assert!(!o.status.success(), "{err}");
            assert!(err.contains("--mode align"), "{flag:?} {q:?}: {err}");
            // With --quantifier salmon, the older --quantifier check fires first and names
            // that option. Without it, the Salmon flag error names the flag.
            let named = if q.is_empty() {
                flag[0]
            } else {
                "--quantifier"
            };
            assert!(err.contains(named), "{flag:?} {q:?}: {err}");
            assert!(!out.exists(), "refused before any file is written");
        }
    }
    // Without the flags, align mode passes the check: the run fails later on the missing FASTQ.
    let o = qrnamap(&base);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(
        !err.contains("--salmon-") && !err.contains("--mode align"),
        "{err}"
    );
    assert!(
        err.contains("error: sample 'A': input error: /x/a_1.fq"),
        "{err}"
    );
}

#[test]
fn salmon_min_score_fraction_invalid_value_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = write_small_sheet(d);
    let out = d.join("out");
    for v in ["0", "1.5", "abc", "-0.1", "NaN", "inf"] {
        let mut a = quant_args(s(&sh), s(&out));
        a.extend(["--quantifier", "salmon", "--salmon-min-score-fraction", v]);
        let o = qrnamap(&a);
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{v}: {err}");
        assert!(err.contains("--salmon-min-score-fraction"), "{v}: {err}");
        assert!(!out.exists(), "{v}: a file was written");
        // `run --slurm` plans nothing for an invalid value.
        let o = qrnamap(&[
            "run",
            "--samples",
            s(&sh),
            "--reference",
            "/x/g.fa",
            "--gff",
            "/x/g.gff3",
            "--mode",
            "quant",
            "--quantifier",
            "salmon",
            "--salmon-min-score-fraction",
            v,
            "--slurm",
            "--out",
            s(&out),
        ]);
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{v}: {err}");
        assert!(err.contains("--salmon-min-score-fraction"), "{v}: {err}");
        assert!(!out.exists(), "{v}: --slurm wrote files");
    }
    // 1 is the upper bound and is valid: the run fails later, on the missing files.
    let mut a = quant_args(s(&sh), s(&out));
    a.extend(["--quantifier", "salmon", "--salmon-min-score-fraction", "1"]);
    let o = qrnamap(&a);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(!err.contains("--salmon-min-score-fraction"), "{err}");
}

#[test]
fn sheet_command_pairs_mates_from_a_folder_and_names_problems() {
    let d = tempfile::tempdir().unwrap();
    for n in [
        "WT1_S1_L001_R1_001.fastq.gz",
        "WT1_S1_L001_R2_001.fastq.gz",
        "SRR9_1.fq.gz",
        "SRR9_2.fq.gz",
        "x_R2.fq",
        "README.txt",
    ] {
        std::fs::write(d.path().join(n), b"").unwrap();
    }
    let dir = d.path().to_str().unwrap();
    let (out, err) = ok(&qrnamap(&["sheet", dir, "--column", "condition"]));
    let p = |n: &str| d.path().join(n).display().to_string();
    assert_eq!(
        out,
        format!(
            "sample\tfastq1\tfastq2\tcondition\nSRR9\t{}\t{}\t\nWT1\t{}\t{}\t\nx\t\t{}\t\n",
            p("SRR9_1.fq.gz"),
            p("SRR9_2.fq.gz"),
            p("WT1_S1_L001_R1_001.fastq.gz"),
            p("WT1_S1_L001_R2_001.fastq.gz"),
            p("x_R2.fq"),
        )
    );
    assert!(err.contains("README.txt: not a FASTQ"), "{err}");
    assert!(err.contains("problem: sample 'x': only mate 2"), "{err}");
    assert!(
        err.contains("3 samples (2 paired-end, 0 single-end, 1 with a problem), 1 files skipped"),
        "{err}"
    );

    // --out does not replace a file; --no-numeric-mates keeps _1/_2 apart.
    let target = d.path().join("samples.tsv");
    ok(&qrnamap(&[
        "sheet",
        &p("SRR9_1.fq.gz"),
        &p("SRR9_2.fq.gz"),
        "--no-numeric-mates",
        "--out",
        target.to_str().unwrap(),
    ]));
    let text = std::fs::read_to_string(&target).unwrap();
    assert!(
        text.contains("SRR9_1\t") && text.contains("SRR9_2\t"),
        "{text}"
    );
    let again = qrnamap(&["sheet", dir, "--out", target.to_str().unwrap()]);
    assert!(!again.status.success());
    let bad = qrnamap(&["sheet", dir, "--column", "1st"]);
    assert!(String::from_utf8_lossy(&bad.stderr).contains("--column \"1st\""));
}
