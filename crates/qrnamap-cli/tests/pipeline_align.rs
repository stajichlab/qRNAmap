//! HISAT2 + featureCounts pipeline and strand check on simulated reads with known truth.
//! Real tools, no mocks. Skips when a tool is missing, fails under `QRNAMAP_REQUIRE_TOOLS=1`.

mod common;

use std::path::Path;
use std::time::Instant;

use common::{records, sample, sim, COUNTS};
use qrnamap_cli::pipeline::{
    align_sample, prepare_align, prepare_quant, quant_sample, read_gene_counts, read_table,
    run_matrix, run_transcripts, strand_check_sample, GenesFrom, Scratch,
};
use qrnamap_core::{RunConfig, SampleSheet, Strandedness};

const ALIGN_TOOLS: [&str; 3] = ["hisat2-align-s", "hisat2-build", "featureCounts"];

fn simulate(dir: &Path, lib: Strandedness, seed: u64) -> (std::path::PathBuf, std::path::PathBuf) {
    let pairs = sim::simulate(&records(), &COUNTS, lib, seed);
    let (r1, r2) = (dir.join("s_1.fq"), dir.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    (r1, r2)
}

#[test]
fn align_counts_match_truth_within_5_percent() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = simulate(d, Strandedness::Unstranded, 11);
    let out = d.join("align");
    let ctx = prepare_align(&common::mini_genome(), &common::mini_gff(), &out, Scratch::new(d.join("tmp")), RunConfig::default(), false).unwrap();
    let res = align_sample(&sample("S1", &[&r1], &[&r2]), &ctx).unwrap();
    let truth = sim::truth_per_gene(&records(), &COUNTS);
    let mut worst: f64 = 0.0;
    for (g, &t) in &truth {
        let c = res.counts.get(g).copied().unwrap_or(0);
        let e = (c as f64 - t as f64).abs() / t as f64;
        eprintln!("  {g}: truth {t}, got {c}, rel error {e:.5}");
        if t >= 200 {
            worst = worst.max(e);
        }
    }
    eprintln!("MEASURED align max relative error (genes >= 200 pairs): {worst:.5}; summary {:?}", res.summary);
    assert!(worst <= 0.05, "max relative error {worst}");
    assert_eq!(res.summary.reads_processed, 10_000);
    assert_eq!(read_gene_counts(&out.join("S1.genecounts.tsv.zst")).unwrap(), res.counts);
    // The SAM lived in scratch and is gone (no --keep-sam).
    assert!(!out.join("S1").join("S1.sam.zst").exists());
    eprintln!("WALL align_counts_match_truth_within_5_percent: {:.2?}", t0.elapsed());
}

#[test]
fn align_wrong_strandedness_reports_low_fraction() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = simulate(d, Strandedness::Reverse, 3);
    let cfg = RunConfig { strandedness: Strandedness::Forward, ..RunConfig::default() };
    let ctx = prepare_align(&common::mini_genome(), &common::mini_gff(), &d.join("a"), Scratch::new(d.join("tmp")), cfg, true).unwrap();
    let res = align_sample(&sample("S1", &[&r1], &[&r2]), &ctx).unwrap();
    eprintln!("MEASURED wrong-strand align summary: {:?}", res.summary);
    assert!(res.summary.fraction_assigned < 0.5, "{:?}", res.summary);
    // --keep-sam keeps a compressed copy in the sample directory.
    assert!(d.join("a").join("S1").join("S1.sam.zst").is_file());
    eprintln!("WALL align_wrong_strandedness_reports_low_fraction: {:.2?}", t0.elapsed());
}

#[test]
fn both_modes_same_gene_rows() {
    if !common::require_tools(&ALIGN_TOOLS) || common::require_tool("kallisto").is_none() {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = simulate(d, Strandedness::Unstranded, 7);
    let s = sample("S1", &[&r1], &[&r2]);
    let samples = vec![s.clone()];
    let (fa, t2g, saf) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"), d.join("genes.saf.zst"));
    run_transcripts(&common::mini_genome(), &common::mini_gff(), &fa, &t2g, Some(&saf), None).unwrap();

    let qdir = d.join("quant");
    let qctx = prepare_quant(&fa, &t2g, &qdir, Scratch::new(d.join("tmp")), RunConfig::default()).unwrap();
    quant_sample(&s, &qctx).unwrap();
    let adir = d.join("align");
    let actx = prepare_align(&common::mini_genome(), &common::mini_gff(), &adir, Scratch::new(d.join("tmp")), RunConfig::default(), false).unwrap();
    align_sample(&s, &actx).unwrap();

    let rows = |p: &Path| -> Vec<String> {
        let t = read_table(p).unwrap();
        t.lines().skip(1).map(|l| l.split('\t').next().unwrap().to_string()).collect()
    };
    // Default: union of the sample genes.
    run_matrix(&qdir, &samples, None).unwrap();
    run_matrix(&adir, &samples, None).unwrap();
    let (q, a) = (rows(&qdir.join("counts.tsv.zst")), rows(&adir.join("counts.tsv.zst")));
    assert_eq!(q, a);
    assert_eq!(q, ["geneA", "geneB", "geneC", "geneR", "geneU"]);
    // Gene lists from the two annotation tables give the same rows too.
    run_matrix(&qdir, &samples, Some((GenesFrom::Tx2Gene, t2g.as_path()))).unwrap();
    run_matrix(&adir, &samples, Some((GenesFrom::Saf, saf.as_path()))).unwrap();
    let (q, a) = (rows(&qdir.join("counts.tsv.zst")), rows(&adir.join("counts.tsv.zst")));
    assert_eq!(q, a);
    assert_eq!(q.len(), 5);
    eprintln!("WALL both_modes_same_gene_rows: {:.2?}", t0.elapsed());
}

fn strand_check(lib: Strandedness, seed: u64) -> qrnamap_core::StrandReport {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = simulate(d, lib, seed);
    let ctx = prepare_align(&common::mini_genome(), &common::mini_gff(), &d.join("sc"), Scratch::new(d.join("tmp")), RunConfig::default(), false).unwrap();
    // Fewer pairs than the files hold: only the first 2,000 are used.
    let r = strand_check_sample(&sample("S1", &[&r1], &[&r2]), &ctx, 2000).unwrap();
    eprintln!("MEASURED strand check {lib:?}: {r:?}");
    assert!(r.sense + r.antisense <= 2000, "{r:?}");
    assert!(r.sense + r.antisense > 1800, "{r:?}");
    r
}

#[test]
fn strand_check_detects_reverse_library() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let r = strand_check(Strandedness::Reverse, 19);
    assert_eq!(r.suggested, Strandedness::Reverse);
    assert!(!r.low_fraction_warning);
    eprintln!("WALL strand_check_detects_reverse_library: {:.2?}", t0.elapsed());
}

#[test]
fn strand_check_detects_forward_and_unstranded_libraries() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    assert_eq!(strand_check(Strandedness::Forward, 29).suggested, Strandedness::Forward);
    let u = strand_check(Strandedness::Unstranded, 31);
    assert_eq!(u.suggested, Strandedness::Unstranded);
    assert!(u.low_fraction_warning);
    eprintln!("WALL strand_check_detects_forward_and_unstranded_libraries: {:.2?}", t0.elapsed());
}

#[test]
fn sheet_with_space_in_path_reaches_tool_as_one_argument() {
    // R18: a path with a space and a `$` is one argv element. No shell sees it.
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path().join("odd dir $HOME");
    std::fs::create_dir(&d).unwrap();
    let (r1, r2) = simulate(&d, Strandedness::Unstranded, 2);
    let s = SampleSheet::parse(&format!(
        "sample\tfastq1\tfastq2\nS1\t{}\t{}\n",
        r1.display(),
        r2.display()
    ))
    .unwrap()
    .samples
    .remove(0);
    let ctx = prepare_align(&common::mini_genome(), &common::mini_gff(), &d.join("out"), Scratch::new(d.join("tmp")), RunConfig::default(), false).unwrap();
    let res = align_sample(&s, &ctx).unwrap();
    assert_eq!(res.summary.reads_processed, 10_000);
}
