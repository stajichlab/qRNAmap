//! Per-sample summary files and `mapping_summary.tsv.zst`. The pure tests need no tool; the
//! two end-to-end tests skip when a tool is missing and fail under `QRNAMAP_REQUIRE_TOOLS=1`.

mod common;

use std::collections::BTreeMap;
use std::path::Path;

use common::{records, sample, sim, COUNTS};
use qrnamap_cli::pipeline::{
    align_sample, build_mapping_summary, prepare_align, prepare_quant, quant_sample, read_table,
    run_matrix, run_transcripts, write_gene_counts, write_output, write_sample_summary, PipelineError,
    QuantContext, SampleSummary, Scratch, MAPPING_SUMMARY,
};
use qrnamap_core::{Mode, Quantifier, RunConfig, SampleSheet, Strandedness};

const ALIGN_TOOLS: [&str; 3] = ["hisat2-align-s", "hisat2-build", "featureCounts"];
const HEADER: &str = "sample\tmode\tquantifier\treads_processed\tfraction_mapped\tfraction_assigned\tstrandedness\tqrnamap_version";

fn summ(n: u64, mapped: f64, assigned: f64) -> SampleSummary {
    SampleSummary { reads_processed: n, fraction_mapped: mapped, fraction_assigned: assigned }
}

fn names(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

fn samples_of(v: &[&str]) -> Vec<qrnamap_core::Sample> {
    let mut sheet = String::from("sample\tfastq1\n");
    for n in v {
        sheet.push_str(&format!("{n}\t/x/{n}.fq\n"));
    }
    SampleSheet::parse(&sheet).unwrap().samples
}

/// Write count and summary files for each sample with the real writers.
fn write_both(d: &Path, v: &[&str]) {
    for (i, n) in v.iter().enumerate() {
        let m: BTreeMap<String, u64> = [("geneA".to_string(), 5u64 + i as u64)].into_iter().collect();
        write_gene_counts(&d.join(format!("{n}.genecounts.tsv.zst")), n, &m).unwrap();
        write_sample_summary(d, n, Mode::Quant, Some(Quantifier::Kallisto), &summ(1000 + i as u64, 0.9, 0.85), Strandedness::Reverse).unwrap();
    }
}

fn summary_rows(d: &Path) -> Vec<Vec<String>> {
    read_table(&d.join(format!("{}.summary.tsv.zst", "S1"))).unwrap().lines().map(|l| l.split('\t').map(String::from).collect()).collect()
}

#[test]
fn summary_header_exact() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_sample_summary(d, "S1", Mode::Align, None, &summ(10_000, 0.9876, 0.5), Strandedness::Forward).unwrap();
    let rows = summary_rows(d);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].join("\t"), HEADER);
    assert_eq!(&rows[1][..7], ["S1", "align", "none", "10000", "0.9876", "0.5", "forward"]);
    assert_eq!(rows[1][7], env!("CARGO_PKG_VERSION"));
    // zstd magic bytes: the file is compressed by its suffix.
    assert_eq!(&std::fs::read(d.join("S1.summary.tsv.zst")).unwrap()[..4], &[0x28, 0xb5, 0x2f, 0xfd]);
    assert!(!d.join("S1.summary.tsv.zst.partial").exists());
}

#[test]
fn summary_removed_before_a_failing_rerun() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_sample_summary(d, "S1", Mode::Quant, Some(Quantifier::Kallisto), &summ(5, 0.5, 0.5), Strandedness::Unstranded).unwrap();
    let old = d.join("S1.summary.tsv.zst");
    assert!(old.is_file());
    // A quant rerun whose FASTQ file does not exist fails before kallisto runs.
    let ctx = QuantContext {
        index: d.join("none.idx"),
        tx2gene: Default::default(),
        out_dir: d.to_path_buf(),
        scratch: Scratch::new(d.join("tmp")),
        config: RunConfig::default(),
    };
    let s = samples_of(&["S1"]).remove(0);
    let err = quant_sample(&s, &ctx).unwrap_err();
    eprintln!("rerun error: {err}");
    assert!(!old.exists(), "old summary survived a failed rerun");
    assert!(!d.join("S1.summary.tsv.zst.partial").exists());
}

#[test]
fn matrix_writes_mapping_summary_in_sheet_order() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S2", "S1", "S3"]);
    let samples = samples_of(&["S3", "S1", "S2"]);
    run_matrix(d, &samples, None).unwrap();
    let text = read_table(&d.join(MAPPING_SUMMARY)).unwrap();
    let lines: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(lines[0].join("\t"), format!("{HEADER}\tstrand_check_suggestion"));
    let order: Vec<&str> = lines[1..].iter().map(|r| r[0]).collect();
    assert_eq!(order, ["S3", "S1", "S2"]);
    // S1 was written second in write_both (i = 1): 1001 reads.
    assert_eq!(&lines[2][1..7], ["quant", "kallisto", "1001", "0.9", "0.85", "reverse"]);
    assert_eq!(lines[2].len(), 9);
    assert_eq!(lines[2][8], "");
    assert!(!d.join(format!("{MAPPING_SUMMARY}.partial")).exists());
}

#[test]
fn matrix_refuses_missing_summary_naming_the_sample() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    assert!(!d.join("counts.tsv.zst").exists());
    std::fs::remove_file(d.join("S2.summary.tsv.zst")).unwrap();
    let err = run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap_err();
    assert!(matches!(&err, PipelineError::MissingSummaryFile { sample, .. } if sample == "S2"), "{err:?}");
    let msg = err.to_string();
    assert!(msg.contains("S2") && msg.contains("S2.summary.tsv.zst") && msg.contains("rerun"), "{msg}");
    // The directory started without a counts file, so this proves the matrix was not written
    // before the summary check failed.
    assert!(!d.join("counts.tsv.zst").exists());
    assert!(!d.join(MAPPING_SUMMARY).exists());
}

#[test]
fn matrix_refuses_summary_for_another_sample() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    // The file S1.summary.tsv.zst holds a row for S9.
    write_sample_summary(d, "S9", Mode::Quant, Some(Quantifier::Kallisto), &summ(1, 0.1, 0.1), Strandedness::Unstranded).unwrap();
    std::fs::rename(d.join("S9.summary.tsv.zst"), d.join("S1.summary.tsv.zst")).unwrap();
    let err = run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap_err();
    assert!(
        matches!(&err, PipelineError::SummaryFileSample { expected, found, .. } if expected == "S1" && found == "S9"),
        "{err:?}"
    );
    let msg = err.to_string();
    assert!(msg.contains("S1") && msg.contains("S9") && msg.contains("rerun"), "{msg}");
    // The directory started without a counts file, so this proves the matrix was not written
    // before the summary check failed.
    assert!(!d.join("counts.tsv.zst").exists());
    assert!(!d.join(MAPPING_SUMMARY).exists());
}

#[test]
fn old_seven_column_summary_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1"]);
    let old = "sample\tmode\treads_processed\tfraction_mapped\tfraction_assigned\tstrandedness\tqrnamap_version\nS1\tquant\t5\t0.5\t0.5\tunstranded\t0.1.0\n";
    write_output(&d.join("S1.summary.tsv.zst"), old.as_bytes()).unwrap();
    let err = build_mapping_summary(d, &names(&["S1"]), None).unwrap_err();
    assert!(
        matches!(&err, PipelineError::BadTable { detail, .. } if detail.contains("header")),
        "{err:?}"
    );
    let t = err.to_string();
    assert!(t.contains("written by an older qrnamap; run that sample again"), "{t}");
}

#[test]
fn mixed_quantifier_rows_keep_their_tool() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    let s = summ(7, 0.5, 0.5);
    write_sample_summary(d, "S2", Mode::Quant, Some(Quantifier::Salmon), &s, Strandedness::Reverse)
        .unwrap();
    run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap();
    let m = read_table(&d.join(MAPPING_SUMMARY)).unwrap();
    let cells: Vec<Vec<&str>> = m.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(cells[0][2], "quantifier");
    assert_eq!(cells[1][2], "kallisto");
    assert_eq!(cells[2][2], "salmon");
}

fn strand_table(d: &Path) -> std::path::PathBuf {
    let p = d.join("strand_check.tsv");
    std::fs::write(
        &p,
        "sample\tsense\tantisense\tambiguous\tsense_fraction\tantisense_fraction\tambiguous_fraction\tsuggested\tlow_fraction_warning\n\
         S2\t90\t5\t1\t0.9474\t0.0526\t0.0104\tforward\tfalse\n",
    )
    .unwrap();
    p
}

#[test]
fn strand_check_suggestion_is_joined_when_present() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    let t = strand_table(d);
    let text = build_mapping_summary(d, &names(&["S1", "S2"]), Some(&t)).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows[0].last().copied(), Some("strand_check_suggestion"));
    assert_eq!(rows[1].last().copied(), Some(""), "S1 has no strand-check row");
    assert_eq!(rows[2].last().copied(), Some("forward"));
    assert!(rows.iter().all(|r| r.len() == 9));
    // Through matrix: <out>/strand_check.tsv is used when it exists.
    run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap();
    let m = read_table(&d.join(MAPPING_SUMMARY)).unwrap();
    assert_eq!(m, text);
}

#[test]
fn strand_check_suggestion_is_empty_when_absent() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    let text = build_mapping_summary(d, &names(&["S1", "S2"]), None).unwrap();
    for (i, l) in text.lines().enumerate() {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f.len(), 9, "{l}");
        if i > 0 {
            assert_eq!(f[8], "", "{l}");
        }
    }
    assert!(text.ends_with('\n'));
}

#[test]
fn stale_mapping_summary_removed_by_matrix_before_writing() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1"]);
    std::fs::write(d.join(MAPPING_SUMMARY), b"old").unwrap();
    std::fs::write(d.join(format!("{MAPPING_SUMMARY}.partial")), b"old").unwrap();
    // S2 has no gene counts: the matrix fails.
    let err = run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap_err();
    assert!(matches!(err, PipelineError::MissingCountFile { .. }), "{err:?}");
    assert!(!d.join(MAPPING_SUMMARY).exists(), "stale mapping summary survived a failing matrix");
    assert!(!d.join(format!("{MAPPING_SUMMARY}.partial")).exists());
}

#[test]
fn matrix_output_other_files_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    let (path, warnings) = run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap();
    assert_eq!(path, d.join("counts.tsv.zst"));
    assert!(warnings.is_empty(), "{warnings:?}");
    let counts = read_table(&path).unwrap();
    assert_eq!(counts.lines().next().unwrap(), "gene\tS1\tS2");
    assert_eq!(counts, "gene\tS1\tS2\ngeneA\t5\t6\n");
}

// ---------------------------------------------------------------------------------------------
// End to end (tools)

fn simulate(d: &Path, lib: Strandedness, seed: u64) -> (std::path::PathBuf, std::path::PathBuf) {
    let pairs = sim::simulate(&records(), &COUNTS, lib, seed);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    (r1, r2)
}

fn check_summary_file(path: &Path, mode: &str, strand: &str, s: &SampleSummary) {
    let t = read_table(path).unwrap();
    let lines: Vec<&str> = t.lines().collect();
    assert_eq!(lines.len(), 2, "{t}");
    assert_eq!(lines[0], HEADER);
    let f: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(f.len(), 8, "{t}");
    assert_eq!(f[0], "S1");
    assert_eq!(f[1], mode);
    assert_eq!(f[2], if mode == "quant" { "kallisto" } else { "none" });
    assert_eq!(f[3].parse::<u64>().unwrap(), s.reads_processed);
    assert_eq!(f[4].parse::<f64>().unwrap(), s.fraction_mapped);
    assert_eq!(f[5].parse::<f64>().unwrap(), s.fraction_assigned);
    assert_eq!(f[6], strand);
    assert_eq!(f[7], env!("CARGO_PKG_VERSION"));
}

#[test]
fn quant_sample_writes_summary_equal_to_printed_numbers() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = simulate(d, Strandedness::Reverse, 11);
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    run_transcripts(&common::mini_genome(), &common::mini_gff(), &fa, &t2g, None, None).unwrap();
    let out = d.join("quant");
    let cfg = RunConfig { strandedness: Strandedness::Reverse, ..RunConfig::default() };
    let ctx = prepare_quant(&fa, &t2g, &out, Scratch::new(d.join("tmp")), cfg, None).unwrap();
    let res = quant_sample(&sample("S1", &[&r1], &[&r2]), &ctx).unwrap();
    eprintln!("quant summary {:?}", res.summary);
    check_summary_file(&out.join("S1.summary.tsv.zst"), "quant", "reverse", &res.summary);
    assert!(!out.join("S1.summary.tsv.zst.partial").exists());
}

#[test]
fn align_sample_writes_summary_equal_to_printed_numbers() {
    if !common::require_tools(&ALIGN_TOOLS) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (r1, r2) = simulate(d, Strandedness::Forward, 11);
    let out = d.join("align");
    let cfg = RunConfig { strandedness: Strandedness::Forward, ..RunConfig::default() };
    let ctx = prepare_align(&common::mini_genome(), &common::mini_gff(), &out, Scratch::new(d.join("tmp")), cfg, false, &[]).unwrap();
    let res = align_sample(&sample("S1", &[&r1], &[&r2]), &ctx).unwrap();
    eprintln!("align summary {:?}", res.summary);
    check_summary_file(&out.join("S1.summary.tsv.zst"), "align", "forward", &res.summary);
    assert!(!out.join("S1.summary.tsv.zst.partial").exists());
}

// ---------------------------------------------------------------------------------------------
// Final review fix wave: I3 (R35-2) and the error paths of the summary readers (R35-8)

#[test]
fn failed_matrix_removes_old_counts_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    // A run directory with an old matrix; S2 has lost its summary (for example a run made
    // before summary files existed).
    std::fs::write(d.join("counts.tsv.zst"), b"old matrix").unwrap();
    std::fs::write(d.join("counts.tsv.zst.partial"), b"old").unwrap();
    std::fs::remove_file(d.join("S2.summary.tsv.zst")).unwrap();
    let err = run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap_err();
    assert!(matches!(&err, PipelineError::MissingSummaryFile { sample, .. } if sample == "S2"), "{err:?}");
    assert!(!d.join("counts.tsv.zst").exists(), "the old matrix survived a failed matrix");
    assert!(!d.join("counts.tsv.zst.partial").exists());
    assert!(!d.join(MAPPING_SUMMARY).exists());
    // A failure before any summary is read (a missing count file) removes it too.
    std::fs::write(d.join("counts.tsv.zst"), b"old matrix").unwrap();
    let err = run_matrix(d, &samples_of(&["S1", "S9"]), None).unwrap_err();
    assert!(matches!(err, PipelineError::MissingCountFile { .. }), "{err:?}");
    assert!(!d.join("counts.tsv.zst").exists());
}

#[test]
fn successful_matrix_still_writes_both() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2"]);
    std::fs::write(d.join("counts.tsv.zst"), b"old matrix").unwrap();
    std::fs::write(d.join(MAPPING_SUMMARY), b"old").unwrap();
    run_matrix(d, &samples_of(&["S1", "S2"]), None).unwrap();
    assert_eq!(read_table(&d.join("counts.tsv.zst")).unwrap(), "gene\tS1\tS2\ngeneA\t5\t6\n");
    let m = read_table(&d.join(MAPPING_SUMMARY)).unwrap();
    assert_eq!(m.lines().count(), 3, "{m}");
    assert!(m.starts_with(HEADER), "{m}");
}

/// Write `text` as `<d>/S1.summary.tsv.zst` and return the error of the summary reader.
fn summary_error(text: &str) -> PipelineError {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_output(&d.join("S1.summary.tsv.zst"), text.as_bytes()).unwrap();
    build_mapping_summary(d, &names(&["S1"]), None).unwrap_err()
}

fn bad_table(e: &PipelineError) -> (usize, String) {
    match e {
        PipelineError::BadTable { line, detail, .. } => (*line, detail.clone()),
        _ => panic!("expected BadTable, got {e:?}"),
    }
}

#[test]
fn sample_summary_reader_error_paths() {
    let row = format!("S1\tquant\tkallisto\t10\t0.9\t0.8\treverse\t{}", env!("CARGO_PKG_VERSION"));
    // Bad header.
    let e = summary_error(&format!("sample\tmode\n{row}\n"));
    assert_eq!(bad_table(&e), (1, "the summary header is not the expected one".into()));
    assert!(e.to_string().contains("S1.summary.tsv.zst line 1:"), "{e}");
    // Missing row.
    let e = summary_error(&format!("{HEADER}\n"));
    assert_eq!(bad_table(&e), (2, "the summary row is missing".into()));
    // Extra row.
    let e = summary_error(&format!("{HEADER}\n{row}\n{row}\n"));
    assert_eq!(bad_table(&e), (2, "expected one row of 8 tab-separated fields".into()));
    // Wrong field count (7 and 9).
    let e = summary_error(&format!("{HEADER}\nS1\tquant\tkallisto\t10\t0.9\t0.8\treverse\n"));
    assert_eq!(bad_table(&e), (2, "expected one row of 8 tab-separated fields".into()));
    let e = summary_error(&format!("{HEADER}\n{row}\textra\n"));
    assert_eq!(bad_table(&e), (2, "expected one row of 8 tab-separated fields".into()));
}

/// Return the error of the strand-check table reader for `text` (S1 has a valid summary).
fn strand_error(text: &str) -> PipelineError {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1"]);
    let t = d.join("strand_check.tsv");
    std::fs::write(&t, text).unwrap();
    build_mapping_summary(d, &names(&["S1"]), Some(&t)).unwrap_err()
}

#[test]
fn strand_suggestion_reader_error_paths() {
    // Column not found.
    let e = strand_error("sample\tsense\nS1\t5\n");
    assert_eq!(bad_table(&e), (1, "column 'suggested' not found".into()));
    let e = strand_error("name\tsuggested\nS1\tforward\n");
    assert_eq!(bad_table(&e), (1, "column 'sample' not found".into()));
    // Too few fields on line 3 (line numbers count from 1 in the file).
    let e = strand_error("sample\tsense\tsuggested\nS1\t5\tforward\nS2\t7\n");
    assert_eq!(bad_table(&e), (3, "row has too few fields".into()));
    assert!(e.to_string().contains("strand_check.tsv line 3:"), "{e}");
    // An empty table: line 0 means the whole file, and the message names no line.
    let e = strand_error("");
    assert_eq!(bad_table(&e), (0, "the table is empty".into()));
    assert!(e.to_string().ends_with("strand_check.tsv: the table is empty"), "{e}");
}

#[test]
fn matrix_warns_when_quantifiers_differ() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    write_both(d, &["S1", "S2", "S3"]);
    let mixed = |w: &str| w.contains("quantifier");
    // One tool: no quantifier warning.
    let (_, warnings) = run_matrix(d, &samples_of(&["S1", "S2", "S3"]), None).unwrap();
    assert!(!warnings.iter().any(|w| mixed(w)), "{warnings:?}");
    // S2 counted with salmon: a warning (not an error) that names both tools and the samples.
    let s = summ(7, 0.5, 0.5);
    write_sample_summary(d, "S2", Mode::Quant, Some(Quantifier::Salmon), &s, Strandedness::Reverse)
        .unwrap();
    let (path, warnings) = run_matrix(d, &samples_of(&["S1", "S2", "S3"]), None).unwrap();
    assert!(path.is_file());
    let w: Vec<&String> = warnings.iter().filter(|w| mixed(w)).collect();
    assert_eq!(w.len(), 1, "{warnings:?}");
    assert!(w[0].contains("kallisto (S1, S3)") && w[0].contains("salmon (S2)"), "{}", w[0]);
}
