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
    let back = read_gene_counts(&out.join("S1.genecounts.tsv.zst")).unwrap();
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
