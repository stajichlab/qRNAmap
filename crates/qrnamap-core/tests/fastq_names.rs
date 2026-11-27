//! Tests of the file-name rules of `fastq_names` (sample sheet rows from FASTQ file names).

use qrnamap_core::fastq_names::{
    group_fastq_names, grouping_sheet_tsv, natural_cmp, safe_sample_name, GroupOptions, Rule,
};
use qrnamap_core::SampleSheet;

fn group(names: &[&str]) -> qrnamap_core::fastq_names::Grouping {
    group_fastq_names(names, &GroupOptions::default())
}

/// (name, fastq1, fastq2, rule) of each row.
fn rows(names: &[&str]) -> Vec<(String, String, String, Rule)> {
    group(names)
        .samples
        .into_iter()
        .map(|r| (r.name, r.fastq1.join(","), r.fastq2.join(","), r.rule))
        .collect()
}

fn row(name: &str, f1: &str, f2: &str, rule: Rule) -> (String, String, String, Rule) {
    (name.into(), f1.into(), f2.into(), rule)
}

#[test]
fn illumina_names_pair_drop_the_s_number_and_sort_naturally() {
    let got = rows(&[
        "WT10_S10_L001_R2_001.fastq.gz",
        "WT10_S10_L001_R1_001.fastq.gz",
        "WT2_S2_L001_R1_001.fastq.gz",
        "WT2_S2_L001_R2_001.fastq.gz",
    ]);
    assert_eq!(
        got,
        vec![
            row(
                "WT2",
                "WT2_S2_L001_R1_001.fastq.gz",
                "WT2_S2_L001_R2_001.fastq.gz",
                Rule::Illumina
            ),
            row(
                "WT10",
                "WT10_S10_L001_R1_001.fastq.gz",
                "WT10_S10_L001_R2_001.fastq.gz",
                Rule::Illumina
            ),
        ]
    );
}

#[test]
fn illumina_lanes_become_comma_lists_in_lane_order() {
    let g = group(&[
        "A_S1_L002_R1_001.fastq.gz",
        "A_S1_L001_R2_001.fastq.gz",
        "A_S1_L001_R1_001.fastq.gz",
        "A_S1_L002_R2_001.fastq.gz",
    ]);
    assert_eq!(g.samples.len(), 1);
    let r = &g.samples[0];
    assert_eq!(
        r.fastq1,
        ["A_S1_L001_R1_001.fastq.gz", "A_S1_L002_R1_001.fastq.gz"]
    );
    assert_eq!(
        r.fastq2,
        ["A_S1_L001_R2_001.fastq.gz", "A_S1_L002_R2_001.fastq.gz"]
    );
    assert!(r.problems.is_empty());
    assert!(
        r.notes.iter().any(|n| n.contains("2 lanes")),
        "{:?}",
        r.notes
    );
}

#[test]
fn illumina_without_lane_and_casava_with_index_sequence() {
    assert_eq!(
        rows(&["B_S3_R1_001.fastq.gz", "B_S3_R2_001.fastq.gz"]),
        vec![row(
            "B",
            "B_S3_R1_001.fastq.gz",
            "B_S3_R2_001.fastq.gz",
            Rule::Illumina
        )]
    );
    assert_eq!(
        rows(&[
            "ctl_ATCACG_L001_R1_001.fastq.gz",
            "ctl_ATCACG_L001_R2_001.fastq.gz"
        ]),
        vec![row(
            "ctl_ATCACG",
            "ctl_ATCACG_L001_R1_001.fastq.gz",
            "ctl_ATCACG_L001_R2_001.fastq.gz",
            Rule::Illumina
        )]
    );
}

#[test]
fn index_reads_undetermined_and_other_files_are_skipped_with_reasons() {
    let g = group(&[
        "A_S1_L001_I1_001.fastq.gz",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        "md5sums.txt",
        "A_S1_L001_R1_001.fastq.gz",
        "A_S1_L001_R1_001.fastq.gz",
    ]);
    let reasons: Vec<(&str, &str)> = g
        .skipped
        .iter()
        .map(|s| (s.path.as_str(), s.reason.as_str()))
        .collect();
    assert_eq!(reasons.len(), 4, "{reasons:?}");
    assert!(reasons[0].1.starts_with("index read"));
    assert!(reasons[1].1.starts_with("Undetermined"));
    assert!(reasons[2].1.starts_with("not a FASTQ"));
    assert_eq!(reasons[3].1, "listed twice");
    assert_eq!(g.samples.len(), 1);
    assert!(g.samples[0].notes.iter().any(|n| n.contains("no mate 2")));
}

#[test]
fn s_number_is_kept_when_dropping_it_would_merge_two_samples() {
    let got = rows(&[
        "rep_S1_R1_001.fastq.gz",
        "rep_S1_R2_001.fastq.gz",
        "rep_S2_R1_001.fastq.gz",
        "rep_S2_R2_001.fastq.gz",
    ]);
    let names: Vec<&str> = got.iter().map(|r| r.0.as_str()).collect();
    assert_eq!(names, ["rep_S1", "rep_S2"]);
}

#[test]
fn read_tags_in_their_common_spellings() {
    for (m1, m2) in [
        ("x_R1.fq.gz", "x_R2.fq.gz"),
        ("x.R1.fastq", "x.R2.fastq"),
        ("x_r1.fq", "x_r2.fq"),
        ("x-R1.fastq.gz", "x-R2.fastq.gz"),
        ("x_read1.fastq.gz", "x_read2.fastq.gz"),
        ("x_L003_R1.fastq.gz", "x_L003_R2.fastq.gz"),
        ("X_R1.FASTQ.GZ", "X_R2.FASTQ.GZ"),
        ("x_R1.fq.zst", "x_R2.fq.zst"),
    ] {
        let got = rows(&[m2, m1]);
        assert_eq!(got.len(), 1, "{m1}: {got:?}");
        assert_eq!(
            (got[0].1.as_str(), got[0].2.as_str(), got[0].3),
            (m1, m2, Rule::ReadTag),
            "{m1}"
        );
        assert!(got[0].0.eq_ignore_ascii_case("x"), "{m1}: {}", got[0].0);
    }
}

#[test]
fn sra_numbers_pair_only_when_both_mates_exist() {
    assert_eq!(
        rows(&[
            "SRR1234567_1.fastq.gz",
            "SRR1234567_2.fastq.gz",
            "SRR7654321_1.fastq.gz"
        ]),
        vec![
            row(
                "SRR1234567",
                "SRR1234567_1.fastq.gz",
                "SRR1234567_2.fastq.gz",
                Rule::Number
            ),
            row("SRR7654321_1", "SRR7654321_1.fastq.gz", "", Rule::Single),
        ]
    );
    // Without numeric mates, ctrl_1 and ctrl_2 are two single-end replicates.
    let g = group_fastq_names(
        &["ctrl_1.fq.gz", "ctrl_2.fq.gz"],
        &GroupOptions {
            numeric_mates: false,
        },
    );
    let names: Vec<&str> = g.samples.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(names, ["ctrl_1", "ctrl_2"]);
    assert!(g
        .samples
        .iter()
        .all(|r| r.rule == Rule::Single && r.fastq2.is_empty()));
}

#[test]
fn missing_mates_and_mismatched_lanes_are_problems() {
    let g = group(&["y_R2.fastq.gz"]);
    assert!(g.samples[0].problems[0].contains("only mate 2"));
    let g = group(&["z_L001_R1.fq.gz", "z_L002_R1.fq.gz", "z_L001_R2.fq.gz"]);
    assert!(
        g.samples[0].problems[0].contains("2 mate 1 files but 1 mate 2"),
        "{:?}",
        g.samples[0].problems
    );
    let g = group(&["w_L001_R1.fq.gz", "w_L002_R2.fq.gz"]);
    assert!(g.samples[0].problems[0].contains("different lanes"));
}

#[test]
fn folders_are_kept_in_paths_and_ignored_in_names() {
    let got = rows(&[
        "/data/run1/a_R1.fq.gz",
        "/data/run1/a_R2.fq.gz",
        "C:\\reads\\b.fastq",
    ]);
    assert_eq!(
        got,
        vec![
            row(
                "a",
                "/data/run1/a_R1.fq.gz",
                "/data/run1/a_R2.fq.gz",
                Rule::ReadTag
            ),
            row("b", "C:\\reads\\b.fastq", "", Rule::Single),
        ]
    );
}

#[test]
fn paths_the_sheet_cannot_hold_are_skipped() {
    let g = group(&["a,b_R1.fq.gz", "tab\there.fq"]);
    assert!(g.samples.is_empty());
    assert!(g.skipped[0].reason.contains("comma"));
    assert!(g.skipped[1].reason.contains("control character"));
}

#[test]
fn unsafe_names_are_made_safe_unique_and_reported() {
    assert_eq!(safe_sample_name("WT 1 (rep)"), "WT_1_rep");
    assert_eq!(safe_sample_name("01-ctrl"), "X01-ctrl");
    assert_eq!(safe_sample_name("__"), "sample");
    assert_eq!(safe_sample_name("matrix.tsv"), "matrix.tsv"); // not a reserved output name
    let g = group(&["WT 1.fq.gz", "WT_1.fastq.gz", "deseq2.vst.csv.fq"]);
    let names: Vec<&str> = g.samples.iter().map(|r| r.name.as_str()).collect();
    assert_eq!(names, ["deseq2.vst.csv_sample", "WT_1", "WT_1-2"]);
    assert!(g.samples[0]
        .notes
        .iter()
        .chain(&g.samples[1].notes)
        .any(|n| n.contains("changed from \"WT 1\"")));
}

#[test]
fn natural_order_compares_digit_runs_as_numbers() {
    let mut v = vec!["S10", "s3", "S2", "S1", "S02", "a"];
    v.sort_by(|a, b| natural_cmp(a, b));
    assert_eq!(v, ["a", "S1", "S02", "S2", "s3", "S10"]);
}

#[test]
fn the_sheet_text_parses_with_the_sheet_parser() {
    let g = group(&[
        "A_S1_L001_R1_001.fastq.gz",
        "A_S1_L002_R1_001.fastq.gz",
        "A_S1_L001_R2_001.fastq.gz",
        "A_S1_L002_R2_001.fastq.gz",
        "SRR1_1.fq.gz",
        "SRR1_2.fq.gz",
        "single.fq",
    ]);
    let text = grouping_sheet_tsv(&g.samples, &[]);
    let sheet = SampleSheet::parse(&text).expect("valid sheet");
    assert_eq!(sheet.samples.len(), 3);
    assert_eq!(sheet.samples[0].name, "A");
    assert_eq!(sheet.samples[0].fastq1.len(), 2);
    assert!(sheet
        .samples
        .iter()
        .any(|s| s.name == "single" && s.fastq2.is_empty()));
    // Empty covariate columns are a template to fill: the header is there, the parser wants values.
    let text = grouping_sheet_tsv(&g.samples, &["condition".into()]);
    assert!(text.starts_with("sample\tfastq1\tfastq2\tcondition\n"));
    assert!(SampleSheet::parse(&text).is_err());
}
