//! kallisto pipeline on simulated reads with known truth. Real tools, no mocks.
//! Skips when kallisto is missing, fails under `QRNAMAP_REQUIRE_TOOLS=1` (ruling R5).

mod common;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

use common::sim;
use qrnamap_cli::pipeline::{
    check_samples, prepare_quant, quant_sample, read_gene_counts, run_transcripts, PipelineError,
    Scratch,
};
use qrnamap_core::{Mode, RunConfig, SampleSheet, Strandedness};

use common::{records, sample, COUNTS};

/// Transcripts and tx2gene from the fixture through the `transcripts` step.
fn transcripts(dir: &Path) -> (PathBuf, PathBuf) {
    let fa = dir.join("tx.fa.zst");
    let t2g = dir.join("tx2gene.tsv.zst");
    run_transcripts(&common::mini_genome(), &common::mini_gff(), &fa, &t2g, None, None).unwrap();
    (fa, t2g)
}

fn max_rel_error(got: &BTreeMap<String, u64>, truth: &BTreeMap<String, u64>) -> f64 {
    let mut worst: f64 = 0.0;
    for (g, &t) in truth {
        let c = got.get(g).copied().unwrap_or(0);
        let e = (c as f64 - t as f64).abs() / t as f64;
        eprintln!("  {g}: truth {t}, got {c}, rel error {e:.5}");
        if t >= 200 {
            worst = worst.max(e);
        }
    }
    worst
}

#[test]
fn quant_recovers_simulated_counts_within_10_percent() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let recs = records();
    let pairs = sim::simulate(&recs, &COUNTS, Strandedness::Unstranded, 11);
    let (r1, r2) = (dir.path().join("s_1.fq"), dir.path().join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let (fa, t2g) = transcripts(dir.path());
    let out = dir.path().join("quant");
    let ctx = prepare_quant(&fa, &t2g, &out, Scratch::new(dir.path().join("tmp")), RunConfig::default()).unwrap();
    let res = quant_sample(&sample("S1", &[&r1], &[&r2]), &ctx).unwrap();
    let truth = sim::truth_per_gene(&recs, &COUNTS);
    let worst = max_rel_error(&res.counts, &truth);
    eprintln!("MEASURED quant max relative error (genes >= 200 pairs): {worst:.5}");
    assert!(worst <= 0.10, "max relative error {worst}");
    assert_eq!(res.summary.reads_processed, 10_000);
    assert!(res.summary.fraction_mapped > 0.9, "{:?}", res.summary);
    // Outputs where the CLI promises them.
    assert!(out.join("S1").join("abundance.tsv").is_file());
    let back = read_gene_counts(&out.join("S1.genecounts.tsv.zst"), "S1").unwrap();
    assert_eq!(back, res.counts);
    eprintln!("WALL quant_recovers_simulated_counts_within_10_percent: {:.2?}", t0.elapsed());
}

#[test]
fn quant_with_lanes_equals_concatenated() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let recs = records();
    let pairs = sim::simulate(&recs, &COUNTS, Strandedness::Unstranded, 23);
    let (a1, a2, b1, b2) = (d.join("a_1.fq"), d.join("a_2.fq"), d.join("b_1.fq"), d.join("b_2.fq"));
    sim::write_fastq(&pairs[..4000], &a1, &a2);
    sim::write_fastq(&pairs[4000..], &b1, &b2);
    // Lane 2 is zstd-compressed: kallisto cannot read it, so the pipeline must stage it.
    let (z1, z2) = (d.join("b_1.fq.zst"), d.join("b_2.fq.zst"));
    for (src, dst) in [(&b1, &z1), (&b2, &z2)] {
        let data = std::fs::read(src).unwrap();
        std::fs::write(dst, zstd::encode_all(&data[..], 3).unwrap()).unwrap();
    }
    // Lane 1 is gzip: kallisto reads it directly (no staging).
    let (g1, g2) = (d.join("a_1.fq.gz"), d.join("a_2.fq.gz"));
    for (src, dst) in [(&a1, &g1), (&a2, &g2)] {
        let data = std::fs::read(src).unwrap();
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut e, &data).unwrap();
        std::fs::write(dst, e.finish().unwrap()).unwrap();
    }
    let (c1, c2) = (d.join("all_1.fq"), d.join("all_2.fq"));
    sim::write_fastq(&pairs, &c1, &c2);

    let (fa, t2g) = transcripts(d);
    let out = d.join("quant");
    let ctx = prepare_quant(&fa, &t2g, &out, Scratch::new(d.join("tmp")), RunConfig::default()).unwrap();
    let lanes = quant_sample(&sample("Lanes", &[&g1, &z1], &[&g2, &z2]), &ctx).unwrap();
    let whole = quant_sample(&sample("Whole", &[&c1], &[&c2]), &ctx).unwrap();
    assert_eq!(lanes.counts, whole.counts);
    assert_eq!(lanes.summary.reads_processed, 10_000);
    let est = |s: &str| {
        let t = std::fs::read_to_string(out.join(s).join("abundance.tsv")).unwrap();
        t.lines().map(|l| l.split('\t').nth(3).unwrap().to_string()).collect::<Vec<_>>()
    };
    assert_eq!(est("Lanes"), est("Whole"), "est_counts differ between lanes and concatenated");
    // The staged decompressed copies are removed after the sample.
    let leftovers: Vec<_> = walk(&d.join("tmp"));
    assert!(leftovers.iter().all(|p| !p.to_string_lossy().ends_with(".fq")), "{leftovers:?}");
    eprintln!("WALL quant_with_lanes_equals_concatenated: {:.2?}", t0.elapsed());
}

fn walk(p: &Path) -> Vec<PathBuf> {
    let mut v = Vec::new();
    if let Ok(rd) = std::fs::read_dir(p) {
        for e in rd.flatten() {
            let path = e.path();
            if path.is_dir() {
                v.extend(walk(&path));
            } else {
                v.push(path);
            }
        }
    }
    v
}

#[test]
fn quant_wrong_strandedness_is_reported() {
    // Review focus 6. Two measured cases on simulated data.
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let recs = records();
    let (fa, t2g) = transcripts(d);
    let cfg = RunConfig { strandedness: Strandedness::Forward, ..RunConfig::default() };
    let ctx = prepare_quant(&fa, &t2g, &d.join("q"), Scratch::new(d.join("tmp")), cfg).unwrap();

    // (a) dUTP-like library run as forward: kallisto pseudoaligns nothing and exits 1.
    let pairs = sim::simulate(&recs, &COUNTS, Strandedness::Reverse, 5);
    let (r1, r2) = (d.join("rev_1.fq"), d.join("rev_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let err = quant_sample(&sample("Rev", &[&r1], &[&r2]), &ctx).unwrap_err();
    assert!(matches!(err, PipelineError::NoReadsPseudoaligned { .. }), "{err:?}");
    let text = err.to_string();
    assert!(text.contains("'Rev'") && text.contains("strand"), "{text}");

    // (b) unstranded library run as forward: about half the pairs pseudoalign.
    let pairs = sim::simulate(&recs, &COUNTS, Strandedness::Unstranded, 6);
    let (u1, u2) = (d.join("un_1.fq"), d.join("un_2.fq"));
    sim::write_fastq(&pairs, &u1, &u2);
    let res = quant_sample(&sample("Un", &[&u1], &[&u2]), &ctx).unwrap();
    eprintln!("MEASURED unstranded-run-as-forward quant summary: {:?}", res.summary);
    assert!(res.summary.fraction_assigned < 0.6, "{:?}", res.summary);
    assert!(res.summary.fraction_assigned < qrnamap_cli::pipeline::LOW_ASSIGNED_FRACTION + 0.1);
    eprintln!("WALL quant_wrong_strandedness_is_reported: {:.2?}", t0.elapsed());
}

#[test]
fn config_error_names_sample_before_any_tool_runs() {
    // Single-end, no fragment length: R20 / review focus 3. No tool is needed.
    let s = SampleSheet::parse("sample\tfastq1\nSE1\t/nonexistent/r.fq\n").unwrap().samples;
    let err = check_samples(&s, &RunConfig::default(), Mode::Quant).unwrap_err();
    assert!(matches!(err, PipelineError::Config { .. }), "{err:?}");
    assert!(err.to_string().contains("SE1"), "{err}");
    // Align mode does not need a fragment length.
    check_samples(&s, &RunConfig::default(), Mode::Align).unwrap();
}

#[test]
fn missing_scratch_is_named_error() {
    let err = Scratch::none().dir("x").unwrap_err();
    assert!(matches!(err, PipelineError::NoScratch), "{err:?}");
    assert!(err.to_string().contains("--tmpdir"), "{err}");
}

#[test]
fn rerun_with_failing_sample_does_not_leave_old_counts() {
    // R24: a rerun that fails must not leave the previous run's gene counts for `matrix`.
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let pairs = sim::simulate(&records(), &COUNTS, Strandedness::Unstranded, 13);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let (fa, t2g) = transcripts(d);
    let out = d.join("q");
    let ctx = prepare_quant(&fa, &t2g, &out, Scratch::new(d.join("tmp")), RunConfig::default()).unwrap();
    let good = sample("S1", &[&r1], &[&r2]);
    quant_sample(&good, &ctx).unwrap();
    let counts = out.join("S1.genecounts.tsv.zst");
    assert!(counts.is_file());
    // Rerun with reads that kallisto pseudoaligns none of.
    let (b1, b2) = (d.join("bad_1.fq"), d.join("bad_2.fq"));
    let junk = "@x/1\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    std::fs::write(&b1, junk).unwrap();
    std::fs::write(&b2, junk.replace("/1", "/2")).unwrap();
    let err = quant_sample(&sample("S1", &[&b1], &[&b2]), &ctx).unwrap_err();
    eprintln!("rerun error: {err}");
    assert!(!counts.exists(), "old gene counts survived a failed rerun");
    // kallisto 0.52.0 writes abundance.tsv before it exits with status 1, so the file can
    // exist here. If it does, it is from the failed run (all est_counts 0), not the old one.
    if let Ok(t) = std::fs::read_to_string(out.join("S1").join("abundance.tsv")) {
        assert!(t.lines().skip(1).all(|l| l.split('\t').nth(3) == Some("0")), "{t}");
    }
    let err = qrnamap_cli::pipeline::run_matrix(&out, &[good], None).unwrap_err();
    assert!(err.to_string().contains("S1.genecounts.tsv.zst"), "{err}");
}

// ---------------------------------------------------------------------------------------------
// C1 (final review): quant validates every mate file to its end before kallisto runs.

fn gz_file(src: &Path, dst: &Path) {
    let data = std::fs::read(src).unwrap();
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(&mut e, &data).unwrap();
    std::fs::write(dst, e.finish().unwrap()).unwrap();
}

/// A quant context on the mini fixture and 10,000 simulated pairs in plain files.
fn quant_setup(d: &Path, seed: u64) -> (qrnamap_cli::pipeline::QuantContext, PathBuf, PathBuf) {
    let pairs = sim::simulate(&records(), &COUNTS, Strandedness::Unstranded, seed);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let (fa, t2g) = transcripts(d);
    let ctx = prepare_quant(&fa, &t2g, &d.join("q"), Scratch::new(d.join("tmp")), RunConfig::default()).unwrap();
    (ctx, r1, r2)
}

#[test]
fn quant_truncated_r2_gz_is_named_error() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (ctx, r1, r2) = quant_setup(d, 31);
    let (g1, g2) = (d.join("s_1.fq.gz"), d.join("s_2.fq.gz"));
    gz_file(&r1, &g1);
    gz_file(&r2, &g2);
    let full = std::fs::read(&g2).unwrap();
    std::fs::write(&g2, &full[..full.len() / 2]).unwrap();
    let err = quant_sample(&sample("S1", &[&g1], &[&g2]), &ctx).unwrap_err();
    let text = err.to_string();
    assert!(
        matches!(&err, PipelineError::Input { sample: Some(s), source: qrnamap_cli::io::InputError::Truncated { .. } } if s == "S1"),
        "{err:?}"
    );
    assert!(text.contains("'S1'") && text.contains("truncated") && text.contains("s_2.fq.gz"), "{text}");
    assert!(!ctx.out_dir.join("S1.genecounts.tsv.zst").exists(), "no gene counts for a truncated input");
}

#[test]
fn quant_mate_record_count_mismatch_is_named_error() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (ctx, r1, r2) = quant_setup(d, 32);
    // Keep the first 2,000 records of mate 2 (8,000 lines).
    let text = std::fs::read_to_string(&r2).unwrap();
    let short: String = text.lines().take(8000).map(|l| format!("{l}\n")).collect();
    let r2s = d.join("short_2.fq");
    std::fs::write(&r2s, short).unwrap();
    let err = quant_sample(&sample("S1", &[&r1], &[&r2s]), &ctx).unwrap_err();
    match &err {
        PipelineError::MateCountMismatch { sample, lane, r1, r2 } => {
            assert_eq!((sample.as_str(), *lane, *r1, *r2), ("S1", 1, 10_000, 2_000));
        }
        other => panic!("expected MateCountMismatch, got {other:?}"),
    }
    let t = err.to_string();
    assert!(t.contains("'S1'") && t.contains("10000") && t.contains("2000"), "{t}");
    assert!(!ctx.out_dir.join("S1.genecounts.tsv.zst").exists());
}

#[test]
fn quant_zero_record_input_is_named_error() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (ctx, _r1, _r2) = quant_setup(d, 33);
    let (e1, e2) = (d.join("empty_1.fq"), d.join("empty_2.fq"));
    std::fs::write(&e1, b"").unwrap();
    std::fs::write(&e2, b"").unwrap();
    let err = quant_sample(&sample("E1", &[&e1], &[&e2]), &ctx).unwrap_err();
    match &err {
        PipelineError::NoReads { sample, file } => {
            assert_eq!(sample, "E1");
            assert!(file.ends_with("empty_1.fq"), "{}", file.display());
        }
        other => panic!("expected NoReads, got {other:?}"),
    }
    assert!(err.to_string().contains("'E1'"), "{err}");
}

#[test]
fn quant_valid_gzip_pair_passes_validation() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (ctx, r1, r2) = quant_setup(d, 34);
    let (g1, g2) = (d.join("s_1.fq.gz"), d.join("s_2.fq.gz"));
    gz_file(&r1, &g1);
    gz_file(&r2, &g2);
    let res = quant_sample(&sample("G1", &[&g1], &[&g2]), &ctx).unwrap();
    assert_eq!(res.summary.reads_processed, 10_000);
    assert!(ctx.out_dir.join("G1.genecounts.tsv.zst").is_file());
}
