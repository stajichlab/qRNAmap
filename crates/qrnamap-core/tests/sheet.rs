//! Sample sheet tests. Expected values are written by hand from the sheet format.

use qrnamap_core::{ColType, SampleSheet, SheetError};

const MINIMAL: &str = "sample\tfastq1\tfastq2\tcondition\tbatch\n\
s1\ta_1.fq.gz\ta_2.fq.gz\tctrl\tb1\n\
s2\tb_1.fq.gz\tb_2.fq.gz\ttreat\tb2\n";

fn errs(tsv: &str) -> Vec<SheetError> {
    SampleSheet::parse(tsv).unwrap_err()
}

#[test]
fn parses_minimal_sheet() {
    let sh = SampleSheet::parse(MINIMAL).unwrap();
    assert_eq!(sh.samples.len(), 2);
    assert_eq!(
        sh.columns,
        vec![
            ("condition".to_string(), ColType::Factor),
            ("batch".to_string(), ColType::Factor)
        ]
    );
    let s = &sh.samples[0];
    assert_eq!(s.name, "s1");
    assert_eq!(s.fastq1, vec!["a_1.fq.gz"]);
    assert_eq!(s.fastq2, vec!["a_2.fq.gz"]);
    assert_eq!(s.covariates.get("condition").map(String::as_str), Some("ctrl"));
    assert_eq!(s.covariates.get("batch").map(String::as_str), Some("b1"));
    assert_eq!(sh.samples[1].name, "s2");
}

#[test]
fn lanes_split_on_comma() {
    let t = "sample\tfastq1\tfastq2\ns1\tl1_1.fq,l2_1.fq\tl1_2.fq,l2_2.fq\n";
    let sh = SampleSheet::parse(t).unwrap();
    assert_eq!(sh.samples[0].fastq1, vec!["l1_1.fq", "l2_1.fq"]);
    assert_eq!(sh.samples[0].fastq2, vec!["l1_2.fq", "l2_2.fq"]);
    assert!(sh.columns.is_empty());
}

#[test]
fn single_end_has_empty_fastq2() {
    let t = "sample\tfastq1\tfastq2\tcondition\ns1\ta.fq\t\tc\ns2\tb_1.fq\tb_2.fq\tc\n";
    let sh = SampleSheet::parse(t).unwrap();
    assert!(sh.samples[0].fastq2.is_empty());
    assert_eq!(sh.samples[1].fastq2, vec!["b_2.fq"]);
    // A sheet with no fastq2 column at all is single-end too.
    let sh = SampleSheet::parse("sample\tfastq1\ns1\ta.fq\n").unwrap();
    assert!(sh.samples[0].fastq2.is_empty());
}

#[test]
fn duplicate_sample_is_error() {
    let t = "sample\tfastq1\ns1\ta.fq\ns1\tb.fq\n";
    assert_eq!(errs(t), vec![SheetError::DuplicateSample("s1".into())]);
}

#[test]
fn name_with_space_rejected() {
    let t = "sample\tfastq1\ns 1\ta.fq\n";
    assert_eq!(errs(t), vec![SheetError::BadName("s 1".into())]);
}

#[test]
fn name_with_quote_rejected() {
    let t = "sample\tfastq1\ns\"1\ta.fq\n";
    assert_eq!(errs(t), vec![SheetError::BadName("s\"1".into())]);
}

#[test]
fn name_with_dollar_rejected() {
    let t = "sample\tfastq1\ns$1\ta.fq\n";
    assert_eq!(errs(t), vec![SheetError::BadName("s$1".into())]);
}

#[test]
fn sample_name_must_start_with_letter() {
    let t = "sample\tfastq1\n1s\ta.fq\n";
    assert_eq!(errs(t), vec![SheetError::BadName("1s".into())]);
    assert!(SampleSheet::parse("sample\tfastq1\nS_1.a-b\ta.fq\n").is_ok());
}

#[test]
fn factor_value_with_unsafe_character_rejected() {
    // space, double quote, single quote, dollar, backtick, backslash, control characters
    // (a raw newline cannot sit inside a cell: it ends the row, so it is not testable here)
    for bad in ["a b", "a\"b", "a'b", "a$b", "a`b", "a\\b", "a\u{b}b", "a\rb", ""] {
        let t = format!("sample\tfastq1\tcondition\ns1\ta.fq\t{bad}\n");
        let e = errs(&t);
        assert_eq!(e.len(), 1, "value {bad:?}");
        match (&e[0], bad.is_empty()) {
            (SheetError::EmptyValue { sample, column }, true) => {
                assert_eq!((sample.as_str(), column.as_str()), ("s1", "condition"));
            }
            (SheetError::BadName(v), false) => assert_eq!(v, bad),
            (other, _) => panic!("value {bad:?} gave {other:?}"),
        }
    }
}

#[test]
fn bad_covariate_column_name_rejected() {
    for bad in ["my cond", "1cond", "a-b", "a$b", "a\"b"] {
        let t = format!("sample\tfastq1\t{bad}\ns1\ta.fq\tx\n");
        assert_eq!(errs(&t), vec![SheetError::BadName(bad.into())], "{bad:?}");
    }
}

#[test]
fn mate_lists_must_match_length() {
    let t = "sample\tfastq1\tfastq2\ns1\ta1.fq,a2.fq\ta1_2.fq\n";
    assert_eq!(
        errs(t),
        vec![SheetError::MateCountMismatch { sample: "s1".into() }]
    );
}

#[test]
fn numeric_header_parses_numbers() {
    let t = "sample\tfastq1\tdose:numeric\tcondition\ns1\ta.fq\t0.5\tc\ns2\tb.fq\t1e2\tc\n";
    let sh = SampleSheet::parse(t).unwrap();
    assert_eq!(
        sh.columns,
        vec![
            ("dose".to_string(), ColType::Numeric),
            ("condition".to_string(), ColType::Factor)
        ]
    );
    assert_eq!(sh.samples[0].covariates.get("dose").map(String::as_str), Some("0.5"));
    assert_eq!(sh.samples[1].covariates.get("dose").map(String::as_str), Some("1e2"));
}

#[test]
fn bad_numeric_value_is_error() {
    for bad in ["abc", "", "inf", "NaN", "1,5"] {
        let t = format!("sample\tfastq1\tdose:numeric\ns1\ta.fq\t{bad}\n");
        assert_eq!(
            errs(&t),
            vec![SheetError::BadNumeric { column: "dose".into(), value: bad.into() }],
            "{bad:?}"
        );
    }
}

#[test]
fn missing_required_columns() {
    assert_eq!(
        errs("sample\tcondition\ns1\tx\n"),
        vec![SheetError::MissingColumn("fastq1".into())]
    );
    assert_eq!(
        errs("fastq1\tcondition\na.fq\tx\n"),
        vec![SheetError::MissingColumn("sample".into())]
    );
}

#[test]
fn all_errors_reported_together() {
    let t = "sample\tfastq1\tfastq2\tdose:numeric\n\
s1\ta1.fq,a2.fq\tb.fq\t1\n\
s1\tc.fq\t\t2\n\
s 3\td.fq\t\tx\n";
    assert_eq!(
        errs(t),
        vec![
            SheetError::MateCountMismatch { sample: "s1".into() },
            SheetError::DuplicateSample("s1".into()),
            SheetError::BadName("s 3".into()),
            SheetError::BadNumeric { column: "dose".into(), value: "x".into() },
        ]
    );
}

#[test]
fn crlf_sheet_parses() {
    let t = MINIMAL.replace('\n', "\r\n");
    let sh = SampleSheet::parse(&t).unwrap();
    assert_eq!(sh, SampleSheet::parse(MINIMAL).unwrap());
    assert_eq!(sh.samples[1].covariates.get("batch").map(String::as_str), Some("b2"));
}

#[test]
fn comments_and_blank_lines_ignored() {
    let t = "# a comment\n\nsample\tfastq1\n# another\ns1\ta.fq\n\n   \ns2\tb.fq\n";
    let sh = SampleSheet::parse(t).unwrap();
    let names: Vec<&str> = sh.samples.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, vec!["s1", "s2"]);
}

#[test]
fn empty_sheet_is_error() {
    assert_eq!(errs(""), vec![SheetError::Empty]);
    assert_eq!(errs("# only a comment\n\n"), vec![SheetError::Empty]);
    assert_eq!(errs("sample\tfastq1\n"), vec![SheetError::Empty]);
}

#[test]
fn duplicate_header_is_error() {
    let t = "sample\tfastq1\tcondition\tcondition\ns1\ta.fq\tx\ty\n";
    assert_eq!(errs(t), vec![SheetError::DuplicateColumn("condition".into())]);
    // The numeric marker does not hide a duplicate.
    let t = "sample\tfastq1\tdose\tdose:numeric\ns1\ta.fq\tx\t1\n";
    assert_eq!(errs(t), vec![SheetError::DuplicateColumn("dose".into())]);
}

#[test]
fn ragged_row_is_error() {
    let t = "sample\tfastq1\tcondition\ns1\ta.fq\n";
    assert_eq!(errs(t), vec![SheetError::BadRow("s1".into())]);
}

#[test]
fn errors_display() {
    let m = SheetError::DuplicateSample("s1".into()).to_string();
    assert!(m.contains("s1"));
    let _: &dyn std::error::Error = &SheetError::Empty;
}
