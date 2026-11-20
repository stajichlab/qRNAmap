//! Recovery tests on the larger synthetic fixture (`common/large.rs`, ruling R24 (3)).
//! Reads come from the generator's own true transcripts, not from `build_transcripts`.
//! Real tools, no mocks. R5: skip when a tool is missing, fail under QRNAMAP_REQUIRE_TOOLS=1.

mod common;

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::Instant;

use common::large::{generate, Large};
use common::sim::{self, Pair, SimParams};
use qrnamap_cli::pipeline::{
    align_sample, prepare_align, prepare_quant, quant_sample, run_transcripts, PipelineError,
    SampleResult, Scratch,
};
use qrnamap_core::{build_transcripts, RunConfig, Strandedness};
use qseq_io::{parse_gff3, read_fasta};

const PARAMS: SimParams = SimParams { read_len: 100, frag_mean: 250, frag_spread: 50 };
const QUANT_TOL: f64 = 0.10;
const ALIGN_TOL: f64 = 0.05;

fn large() -> &'static Large {
    static L: OnceLock<Large> = OnceLock::new();
    L.get_or_init(generate)
}

/// Genome FASTA and GFF3 of the large fixture, written into `d`.
fn write_reference(d: &Path) -> (PathBuf, PathBuf) {
    let (fa, gff) = (d.join("large.fa"), d.join("large.gff3"));
    std::fs::write(&fa, large().fasta()).unwrap();
    std::fs::write(&gff, large().gff3()).unwrap();
    (fa, gff)
}

fn pairs(lib: Strandedness, seed: u64) -> Vec<Pair> {
    sim::simulate_with(&large().truth, &large().counts, lib, seed, PARAMS)
}

/// Gzip with two members (first half, second half): spec 11f multi-member input.
fn gzip_two_members(data: &[u8]) -> Vec<u8> {
    let cut = data.len() / 2;
    let mut out = Vec::new();
    for part in [&data[..cut], &data[cut..]] {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(part).unwrap();
        out.extend(e.finish().unwrap());
    }
    out
}

/// Two lanes: lane 1 multi-member gzip, lane 2 zstd. Returns (R1 lanes, R2 lanes).
fn write_compressed_lanes(d: &Path, p: &[Pair]) -> (Vec<PathBuf>, Vec<PathBuf>) {
    let cut = p.len() * 2 / 5;
    let (a1, a2, b1, b2) = (d.join("a_1.fq"), d.join("a_2.fq"), d.join("b_1.fq"), d.join("b_2.fq"));
    sim::write_fastq(&p[..cut], &a1, &a2);
    sim::write_fastq(&p[cut..], &b1, &b2);
    let conv = |src: &Path, ext: &str| {
        let data = std::fs::read(src).unwrap();
        let dst = PathBuf::from(format!("{}.{ext}", src.display()));
        let bytes = if ext == "gz" { gzip_two_members(&data) } else { zstd::encode_all(&data[..], 3).unwrap() };
        std::fs::write(&dst, bytes).unwrap();
        std::fs::remove_file(src).unwrap();
        dst
    };
    (vec![conv(&a1, "gz"), conv(&b1, "zst")], vec![conv(&a2, "gz"), conv(&b2, "zst")])
}

fn sample_of(name: &str, r1: &[PathBuf], r2: &[PathBuf]) -> qrnamap_core::Sample {
    let a: Vec<&Path> = r1.iter().map(PathBuf::as_path).collect();
    let b: Vec<&Path> = r2.iter().map(PathBuf::as_path).collect();
    common::sample(name, &a, &b)
}

/// Max relative error over genes with at least 200 true pairs; prints every gene.
fn max_error(label: &str, got: &BTreeMap<String, u64>) -> f64 {
    let truth = large().truth_per_gene();
    let mut worst: f64 = 0.0;
    for (g, &t) in &truth {
        let c = got.get(g).copied().unwrap_or(0);
        let e = (c as f64 - t as f64).abs() / t as f64;
        eprintln!("  {label} {g}: truth {t}, got {c}, rel error {e:.5}");
        if t >= 200 {
            worst = worst.max(e);
        }
    }
    eprintln!("MEASURED {label} max relative error (genes >= 200 pairs): {worst:.5}");
    worst
}

fn run_quant(d: &Path, lib: Strandedness, flag: Strandedness, r1: &[PathBuf], r2: &[PathBuf]) -> Result<SampleResult, PipelineError> {
    let (fa, gff) = write_reference(d);
    let (tx, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    run_transcripts(&fa, &gff, &tx, &t2g, None, None).unwrap();
    let cfg = RunConfig { strandedness: flag, ..RunConfig::default() };
    let ctx = prepare_quant(&tx, &t2g, &d.join(format!("q_{lib:?}_{flag:?}")), Scratch::new(d.join("tmp")), cfg).unwrap();
    quant_sample(&sample_of("S1", r1, r2), &ctx)
}

fn run_align(d: &Path, flag: Strandedness, samples: &[(&str, &[PathBuf], &[PathBuf])]) -> Vec<SampleResult> {
    let (fa, gff) = write_reference(d);
    let cfg = RunConfig { strandedness: flag, ..RunConfig::default() };
    let ctx = prepare_align(&fa, &gff, &d.join(format!("a_{flag:?}")), Scratch::new(d.join("tmp")), cfg, false).unwrap();
    samples.iter().map(|(n, r1, r2)| align_sample(&sample_of(n, r1, r2), &ctx).unwrap()).collect()
}

const QUANT_TOOLS: [&str; 1] = ["kallisto"];
const ALIGN_TOOLS: [&str; 3] = ["hisat2-align-s", "hisat2-build", "featureCounts"];

#[test]
fn large_builder_matches_generator_transcripts() {
    // Fixture-independent check: the builder, fed the generated FASTA and GFF3, must give the
    // generator's own transcript sequences.
    let l = large();
    let genome = read_fasta(l.fasta().as_bytes()).unwrap();
    let model = parse_gff3(l.gff3().as_bytes()).unwrap();
    let built = build_transcripts(&genome, &model).unwrap();
    assert_eq!(built.len(), l.truth.len());
    for (b, t) in built.iter().zip(&l.truth) {
        assert_eq!((&b.id, &b.gene, &b.biotype), (&t.id, &t.gene, &t.biotype));
        assert!(b.seq == t.seq, "sequence of {} differs from the generator's truth", t.id);
    }
    let total: u64 = l.counts.iter().sum();
    assert_eq!(total, 40_000);
    eprintln!(
        "large fixture: {} genes, {} transcripts, lengths {}..{}",
        l.genes.len(),
        l.truth.len(),
        l.truth.iter().map(|t| t.seq.len()).min().unwrap(),
        l.truth.iter().map(|t| t.seq.len()).max().unwrap()
    );
}

#[test]
fn large_quant_recovers_truth_within_10_percent() {
    if !common::require_tools(&QUANT_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let p = pairs(Strandedness::Unstranded, 101);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&p, &r1, &r2);
    let res = run_quant(d, Strandedness::Unstranded, Strandedness::Unstranded, &[r1], &[r2]).unwrap();
    let worst = max_error("quant", &res.counts);
    assert!(worst <= QUANT_TOL, "max relative error {worst}");
    assert_eq!(res.summary.reads_processed, 40_000);
    eprintln!("WALL large_quant_recovers_truth_within_10_percent: {:.2?}", t0.elapsed());
}

#[test]
fn large_align_recovers_truth_within_5_percent() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = write_compressed_lanes(d, &pairs(Strandedness::Unstranded, 102));
    let res = run_align(d, Strandedness::Unstranded, &[("S1", &r1, &r2)]).remove(0);
    eprintln!("summary {:?}", res.summary);
    let fc = std::fs::read_to_string(d.join("a_Unstranded/S1/featurecounts.summary")).unwrap();
    eprintln!("featureCounts summary (non-zero rows):");
    for l in fc.lines().filter(|l| !l.ends_with("\t0")) {
        eprintln!("  {l}");
    }
    let worst = max_error("align", &res.counts);
    assert!(worst <= ALIGN_TOL, "max relative error {worst}");
    assert_eq!(res.summary.reads_processed, 40_000);
    eprintln!("WALL large_align_recovers_truth_within_5_percent: {:.2?}", t0.elapsed());
}

fn stranded_full_counts(lib: Strandedness, seed: u64) {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = write_compressed_lanes(d, &pairs(lib, seed));
    let q = run_quant(d, lib, lib, &r1, &r2).unwrap();
    let wq = max_error(&format!("quant {lib:?} as {lib:?}"), &q.counts);
    assert!(wq <= QUANT_TOL, "quant max relative error {wq}");
    let a = run_align(d, lib, &[("S1", &r1, &r2)]).remove(0);
    eprintln!("align {lib:?} summary {:?}", a.summary);
    let wa = max_error(&format!("align {lib:?} as {lib:?}"), &a.counts);
    assert!(wa <= ALIGN_TOL, "align max relative error {wa}");
}

#[test]
fn large_reverse_library_with_reverse_flag_gives_full_counts_in_quant_and_align() {
    if !common::require_tools(&QUANT_TOOLS) || !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    stranded_full_counts(Strandedness::Reverse, 103);
    eprintln!("WALL large_reverse_library_with_reverse_flag_gives_full_counts_in_quant_and_align: {:.2?}", t0.elapsed());
}

#[test]
fn large_forward_library_with_forward_flag_gives_full_counts_in_quant_and_align() {
    if !common::require_tools(&QUANT_TOOLS) || !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    stranded_full_counts(Strandedness::Forward, 104);
    eprintln!("WALL large_forward_library_with_forward_flag_gives_full_counts_in_quant_and_align: {:.2?}", t0.elapsed());
}

#[test]
fn large_reverse_library_run_as_forward_loses_counts() {
    if !common::require_tools(&QUANT_TOOLS) || !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = write_compressed_lanes(d, &pairs(Strandedness::Reverse, 105));
    match run_quant(d, Strandedness::Reverse, Strandedness::Forward, &r1, &r2) {
        Err(PipelineError::NoReadsPseudoaligned { .. }) => {
            eprintln!("MEASURED quant reverse-as-forward: NoReadsPseudoaligned (kallisto exit 1)")
        }
        Ok(r) => {
            eprintln!("MEASURED quant reverse-as-forward summary {:?}", r.summary);
            assert!(r.summary.fraction_assigned < 0.2, "{:?}", r.summary);
        }
        Err(e) => panic!("unexpected error: {e}"),
    }
    let a = run_align(d, Strandedness::Forward, &[("S1", &r1, &r2)]).remove(0);
    eprintln!("MEASURED align reverse-as-forward summary {:?}", a.summary);
    assert!(a.summary.fraction_assigned < 0.2, "{:?}", a.summary);
    eprintln!("WALL large_reverse_library_run_as_forward_loses_counts: {:.2?}", t0.elapsed());
}

#[test]
fn large_align_lanes_equal_concatenated() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let t0 = Instant::now();
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let p = pairs(Strandedness::Unstranded, 106);
    let (l1, l2) = write_compressed_lanes(d, &p);
    let (c1, c2) = (d.join("all_1.fq"), d.join("all_2.fq"));
    sim::write_fastq(&p, &c1, &c2);
    let res = run_align(d, Strandedness::Unstranded, &[("Lanes", &l1, &l2), ("Whole", &[c1], &[c2])]);
    assert_eq!(res[0].counts, res[1].counts);
    assert_eq!(res[0].summary, res[1].summary);
    eprintln!("WALL large_align_lanes_equal_concatenated: {:.2?}", t0.elapsed());
}
