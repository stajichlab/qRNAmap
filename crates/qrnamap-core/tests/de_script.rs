use qrnamap_core::de::{check_design, coldata_tsv, deseq2_script, DeError, Shrink};
use qrnamap_core::{ColType, Contrast, CountMatrix, Design, SampleSheet};

const SHEET: &str = "sample\tfastq1\tfastq2\tcondition\tbatch\tage:numeric\n\
a1\t/data/secret_reads_a1_R1.fastq.gz\t/data/secret_reads_a1_R2.fastq.gz\tcontrol\tb1\t3\n\
a2\t/data/secret_reads_a2_R1.fastq.gz\t/data/secret_reads_a2_R2.fastq.gz\tcontrol\tb2\t4\n\
b1\t/data/secret_reads_b1_R1.fastq.gz\t/data/secret_reads_b1_R2.fastq.gz\ttreated\tb1\t5\n\
b2\t/data/secret_reads_b2_R1.fastq.gz\t/data/secret_reads_b2_R2.fastq.gz\ttreated\tb2\t6\n";

fn sheet() -> SampleSheet {
    SampleSheet::parse(SHEET).unwrap()
}

fn contrast(sh: &SampleSheet) -> Contrast {
    Contrast::parse("condition,treated,control", &sh.columns).unwrap()
}

fn matrix(sh: &SampleSheet) -> CountMatrix {
    CountMatrix {
        genes: vec!["g1".into(), "g2".into()],
        samples: sh.samples.iter().map(|s| s.name.clone()).collect(),
        counts: vec![vec![10, 12, 30, 28], vec![5, 6, 7, 8]],
    }
}

fn script(design: &str, shrink: Option<Shrink>) -> String {
    let sh = sheet();
    let d = Design::parse(design, &sh.columns).unwrap();
    let c = contrast(&sh);
    let cols: Vec<(String, ColType)> = sh
        .columns
        .iter()
        .filter(|(n, _)| d.terms().contains(n))
        .cloned()
        .collect();
    deseq2_script(&d, &c, &cols, shrink)
}

#[test]
fn script_default_has_no_shrink() {
    let s = script("~ condition", None);
    assert!(!s.contains("lfcShrink"), "{s}");
    assert!(s.contains("DESeq(dds"));
}

#[test]
fn script_with_shrink_applies_it_after_results_are_written() {
    let s = script("~ condition", Some(Shrink::Normal));
    let w = s.find("results.csv").unwrap();
    let l = s.find("lfcShrink").unwrap();
    assert!(w < l);
    assert!(s.contains("results.shrunk.csv"));
}

#[test]
#[allow(non_snake_case)]
fn script_uses_varianceStabilizingTransformation() {
    let s = script("~ condition", None);
    assert!(s.contains("varianceStabilizingTransformation("));
    assert!(!s.contains("vst("));
}

#[test]
fn script_sets_reference_level_explicitly() {
    let s = script("~ condition", None);
    assert!(s.contains("levels = c(\"control\""), "{s}");
    assert!(s.contains("contrast_spec <- c(\"condition\", \"treated\", \"control\")"), "{s}");
}

#[test]
fn numeric_column_not_converted_to_factor() {
    let s = script("~ age + condition", None);
    assert!(s.contains("\"age\" = \"numeric\""), "{s}");
    assert!(!s.contains("[[\"age\"]] <- factor"), "{s}");
    assert!(s.contains("[[\"condition\"]] <- factor"), "{s}");
}

#[test]
fn script_contains_no_user_text() {
    let sh = sheet();
    // Values Task 8 rejects, built by hand.
    let d = Design::parse("~ condition", &sh.columns).unwrap();
    let c = Contrast { factor: "condition".into(), level: "tr\"eat".into(), reference: "con\\trol\"x".into() };
    let cols = vec![("condition".to_string(), ColType::Factor)];
    let s = deseq2_script(&d, &c, &cols, None);
    assert!(s.contains("\"con\\\\trol\\\"x\""), "{s}");
    assert!(s.contains("\"tr\\\"eat\""), "{s}");
    assert!(!s.contains("con\\trol\"x"));
    // The script holds no sample names and no paths. Paths arrive through commandArgs.
    let s2 = script("~ condition", None);
    for bad in ["secret", ".fastq", "/data", "a1", "b2"] {
        assert!(!s2.contains(bad), "script contains {bad}: {s2}");
    }
    assert!(!s2.contains('/'));
    assert!(s2.contains("commandArgs(trailingOnly = TRUE)"));
}

#[test]
fn coldata_has_no_fastq_columns() {
    let sh = sheet();
    let d = Design::parse("~ batch + condition", &sh.columns).unwrap();
    let t = coldata_tsv(&sh, &d);
    assert_eq!(t.lines().next().unwrap(), "sample\tbatch\tcondition");
    assert!(!t.contains("fastq") && !t.contains("secret") && !t.contains("age"));
    assert!(t.contains("a1\tb1\tcontrol\n"));
    assert_eq!(t.lines().count(), 5);
}

#[test]
fn check_design_accepts_good_input() {
    let sh = sheet();
    check_design(&sh, &contrast(&sh), &matrix(&sh)).unwrap();
}

#[test]
fn one_replicate_group_refused() {
    let mut sh = sheet();
    sh.samples[1].covariates.insert("condition".into(), "other".into());
    let m = matrix(&sh);
    let e = check_design(&sh, &contrast(&sh), &m).unwrap_err();
    assert_eq!(e, DeError::TooFewReplicates { group: "control".into(), n: 1 });
}

#[test]
fn all_zero_sample_refused() {
    let sh = sheet();
    let mut m = matrix(&sh);
    for row in &mut m.counts {
        row[2] = 0;
    }
    assert_eq!(check_design(&sh, &contrast(&sh), &m).unwrap_err(), DeError::AllZeroSample("b1".into()));
}

#[test]
fn level_not_in_data_refused() {
    let sh = sheet();
    let c = Contrast { factor: "condition".into(), level: "ghost".into(), reference: "control".into() };
    assert_eq!(
        check_design(&sh, &c, &matrix(&sh)).unwrap_err(),
        DeError::LevelNotInData { factor: "condition".into(), level: "ghost".into() }
    );
    let c = Contrast { factor: "condition".into(), level: "treated".into(), reference: "ghost".into() };
    assert_eq!(
        check_design(&sh, &c, &matrix(&sh)).unwrap_err(),
        DeError::LevelNotInData { factor: "condition".into(), level: "ghost".into() }
    );
}

#[test]
fn sample_set_mismatch_refused() {
    let sh = sheet();
    let mut m = matrix(&sh);
    m.samples[3] = "zzz".into();
    assert!(matches!(check_design(&sh, &contrast(&sh), &m).unwrap_err(), DeError::SampleMismatch(_)));
    let mut m = matrix(&sh);
    m.samples.pop();
    for r in &mut m.counts {
        r.pop();
    }
    assert!(matches!(check_design(&sh, &contrast(&sh), &m).unwrap_err(), DeError::SampleMismatch(_)));
}

#[test]
fn reference_equals_level_refused() {
    let sh = sheet();
    let c = Contrast { factor: "condition".into(), level: "control".into(), reference: "control".into() };
    assert_eq!(check_design(&sh, &c, &matrix(&sh)).unwrap_err(), DeError::SameLevelAndReference("control".into()));
}

#[test]
fn numeric_contrast_factor_refused() {
    let sh = sheet();
    let c = Contrast { factor: "age".into(), level: "3".into(), reference: "4".into() };
    assert_eq!(check_design(&sh, &c, &matrix(&sh)).unwrap_err(), DeError::ContrastNotFactor("age".into()));
}

#[test]
fn apeglm_coefficient_is_matched_through_make_names() {
    let s = script("~ condition", Some(Shrink::Apeglm));
    assert!(s.contains("make.names(paste0(\"condition\", \"_\", \"treated\", \"_vs_\", \"control\"))"), "{s}");
    assert!(s.contains("resultsNames(dds)"));
    assert!(s.contains("paste(resultsNames(dds), collapse"), "{s}");
}

#[test]
fn script_records_vst_blind_false() {
    let s = script("~ condition", None);
    assert!(s.contains("vst: blind=FALSE"));
}

#[test]
fn count_above_i32_max_is_named_error() {
    // R29-9: R stores counts as 32-bit integers; a larger count is refused before R runs.
    let sh = sheet();
    let mut m = matrix(&sh);
    m.counts[1][2] = 2_147_483_647;
    check_design(&sh, &contrast(&sh), &m).unwrap();
    m.counts[1][2] = 2_147_483_648;
    let e = check_design(&sh, &contrast(&sh), &m).unwrap_err();
    assert_eq!(e, DeError::CountTooLarge { gene: "g2".into(), sample: "b1".into() });
    let t = e.to_string();
    assert!(t.contains("g2") && t.contains("b1") && t.contains("2147483647"), "{t}");
}
