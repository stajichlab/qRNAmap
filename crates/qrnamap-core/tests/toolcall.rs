#![allow(non_snake_case)] // test names F/R/FR/RF are pinned by the brief

use qrnamap_core::{
    argv_native, feature_counts, hisat2_align, hisat2_build, kallisto_index, kallisto_quant,
    ConfigError, Mode, Quantifier, RunConfig, Sample, Strandedness,
};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn s(v: &[&str]) -> Vec<String> {
    v.iter().map(|x| x.to_string()).collect()
}
fn p(v: &[&str]) -> Vec<PathBuf> {
    v.iter().map(PathBuf::from).collect()
}
fn sample(f1: &[&str], f2: &[&str]) -> Sample {
    Sample {
        name: "S1".into(),
        fastq1: s(f1),
        fastq2: s(f2),
        covariates: BTreeMap::new(),
    }
}
fn cfg(st: Strandedness) -> RunConfig {
    RunConfig { strandedness: st, ..RunConfig::default() }
}

#[test]
fn kallisto_index_argv() {
    let c = kallisto_index("tx.fa", "tx.idx");
    assert_eq!(argv_native(&c), s(&["kallisto", "index", "-i", "tx.idx", "tx.fa"]));
    assert_eq!(c.inputs, p(&["tx.fa"]));
    assert_eq!(c.outputs, p(&["tx.idx"]));
}

#[test]
fn kallisto_paired_unstranded() {
    let c = kallisto_quant(&sample(&["a_1.fq"], &["a_2.fq"]), "i.idx", "out", &cfg(Strandedness::Unstranded)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["kallisto", "quant", "-i", "i.idx", "-o", "out", "--plaintext", "-t", "1", "a_1.fq", "a_2.fq"])
    );
    assert_eq!(c.inputs, p(&["a_1.fq", "a_2.fq", "i.idx"]));
    assert_eq!(c.outputs, p(&["out/abundance.tsv", "out/run_info.json"]));
}

#[test]
fn kallisto_paired_reverse_adds_rf_stranded() {
    let c = kallisto_quant(&sample(&["a_1.fq"], &["a_2.fq"]), "i.idx", "out", &cfg(Strandedness::Reverse)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["kallisto", "quant", "-i", "i.idx", "-o", "out", "--plaintext", "-t", "1", "--rf-stranded", "a_1.fq", "a_2.fq"])
    );
}

#[test]
fn kallisto_forward_adds_fr_stranded() {
    let c = kallisto_quant(&sample(&["a_1.fq"], &["a_2.fq"]), "i.idx", "out", &cfg(Strandedness::Forward)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["kallisto", "quant", "-i", "i.idx", "-o", "out", "--plaintext", "-t", "1", "--fr-stranded", "a_1.fq", "a_2.fq"])
    );
}

#[test]
fn kallisto_single_end_needs_fragment_length() {
    let r = kallisto_quant(&sample(&["a.fq"], &[]), "i.idx", "out", &RunConfig::default());
    assert_eq!(r, Err(ConfigError::FragmentLengthRequired { sample: "S1".into() }));
}

#[test]
fn kallisto_single_end_uses_l_and_s() {
    let mut c0 = cfg(Strandedness::Unstranded);
    c0.fragment_len = Some((200.0, 20.0));
    let c = kallisto_quant(&sample(&["a.fq"], &[]), "i.idx", "out", &c0).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["kallisto", "quant", "-i", "i.idx", "-o", "out", "--plaintext", "-t", "1", "--single", "-l", "200", "-s", "20", "a.fq"])
    );
    c0.fragment_len = Some((187.5, 24.25));
    let c = kallisto_quant(&sample(&["a.fq"], &[]), "i.idx", "out", &c0).unwrap();
    let a = argv_native(&c);
    assert_eq!(&a[10..14], &s(&["-l", "187.5", "-s", "24.25"])[..]);
}

#[test]
fn hisat2_pair_reverse_is_RF() {
    let c = hisat2_align(&sample(&["a_1.fq"], &["a_2.fq"]), "ix", "ss.txt", "o.sam", &cfg(Strandedness::Reverse)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["hisat2", "-p", "1", "-x", "ix", "--known-splicesite-infile", "ss.txt", "--no-temp-splicesite",
            "--max-intronlen", "5000", "--rna-strandness", "RF", "-1", "a_1.fq", "-2", "a_2.fq", "-S", "o.sam"])
    );
    assert_eq!(c.inputs, p(&["a_1.fq", "a_2.fq", "ss.txt", "ix"]));
    assert_eq!(c.outputs, p(&["o.sam"]));
}

#[test]
fn hisat2_pair_forward_is_FR() {
    let c = hisat2_align(&sample(&["a_1.fq"], &["a_2.fq"]), "ix", "ss.txt", "o.sam", &cfg(Strandedness::Forward)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["hisat2", "-p", "1", "-x", "ix", "--known-splicesite-infile", "ss.txt", "--no-temp-splicesite",
            "--max-intronlen", "5000", "--rna-strandness", "FR", "-1", "a_1.fq", "-2", "a_2.fq", "-S", "o.sam"])
    );
}

#[test]
fn hisat2_single_forward_is_F() {
    let c = hisat2_align(&sample(&["a.fq"], &[]), "ix", "ss.txt", "o.sam", &cfg(Strandedness::Forward)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["hisat2", "-p", "1", "-x", "ix", "--known-splicesite-infile", "ss.txt", "--no-temp-splicesite",
            "--max-intronlen", "5000", "--rna-strandness", "F", "-U", "a.fq", "-S", "o.sam"])
    );
}

#[test]
fn hisat2_single_reverse_is_R() {
    let c = hisat2_align(&sample(&["a.fq"], &[]), "ix", "ss.txt", "o.sam", &cfg(Strandedness::Reverse)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["hisat2", "-p", "1", "-x", "ix", "--known-splicesite-infile", "ss.txt", "--no-temp-splicesite",
            "--max-intronlen", "5000", "--rna-strandness", "R", "-U", "a.fq", "-S", "o.sam"])
    );
}

#[test]
fn hisat2_has_no_temp_splicesite_and_known_sites() {
    let mut c0 = cfg(Strandedness::Unstranded);
    c0.threads = 4;
    c0.max_intron = 3000;
    let c = hisat2_align(&sample(&["a_1.fq"], &["a_2.fq"]), "ix", "ss.txt", "o.sam", &c0).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["hisat2", "-p", "4", "-x", "ix", "--known-splicesite-infile", "ss.txt", "--no-temp-splicesite",
            "--max-intronlen", "3000", "-1", "a_1.fq", "-2", "a_2.fq", "-S", "o.sam"])
    );
}

#[test]
fn hisat2_build_outputs_eight_files() {
    let c = hisat2_build("ref.fa", "idx/ref", &RunConfig::default()).unwrap();
    assert_eq!(argv_native(&c), s(&["hisat2-build", "-p", "1", "ref.fa", "idx/ref"]));
    assert_eq!(c.inputs, p(&["ref.fa"]));
    let want: Vec<PathBuf> = (1..=8).map(|i| PathBuf::from(format!("idx/ref.{i}.ht2"))).collect();
    assert_eq!(c.outputs, want);
    assert_eq!(c.outputs.len(), 8);
}

#[test]
fn featurecounts_paired_has_p_and_countreadpairs() {
    let c = feature_counts("a.sam", "g.saf", "counts.txt", true, &cfg(Strandedness::Unstranded)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["featureCounts", "-T", "1", "-F", "SAF", "-a", "g.saf", "-p", "--countReadPairs", "-s", "0", "-o", "counts.txt", "a.sam"])
    );
    assert_eq!(c.inputs, p(&["a.sam", "g.saf"]));
}

#[test]
fn featurecounts_single_end_has_no_p() {
    let c = feature_counts("a.sam", "g.saf", "counts.txt", false, &cfg(Strandedness::Forward)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["featureCounts", "-T", "1", "-F", "SAF", "-a", "g.saf", "-s", "1", "-o", "counts.txt", "a.sam"])
    );
}

#[test]
fn featurecounts_reverse_is_s2() {
    let c = feature_counts("a.sam", "g.saf", "counts.txt", true, &cfg(Strandedness::Reverse)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["featureCounts", "-T", "1", "-F", "SAF", "-a", "g.saf", "-p", "--countReadPairs", "-s", "2", "-o", "counts.txt", "a.sam"])
    );
}

#[test]
fn featurecounts_outputs_include_summary() {
    let c = feature_counts("a.sam", "g.saf", "counts.txt", true, &RunConfig::default()).unwrap();
    assert_eq!(c.outputs, p(&["counts.txt", "counts.txt.summary"]));
}

#[test]
fn lanes_kallisto_adjacent_pairs_hisat2_comma_lists() {
    let smp = sample(&["a1.fq", "b1.fq"], &["a2.fq", "b2.fq"]);
    let k = kallisto_quant(&smp, "i.idx", "out", &RunConfig::default()).unwrap();
    assert_eq!(
        argv_native(&k),
        s(&["kallisto", "quant", "-i", "i.idx", "-o", "out", "--plaintext", "-t", "1", "a1.fq", "a2.fq", "b1.fq", "b2.fq"])
    );
    let h = hisat2_align(&smp, "ix", "ss.txt", "o.sam", &RunConfig::default()).unwrap();
    assert_eq!(
        argv_native(&h),
        s(&["hisat2", "-p", "1", "-x", "ix", "--known-splicesite-infile", "ss.txt", "--no-temp-splicesite",
            "--max-intronlen", "5000", "-1", "a1.fq,b1.fq", "-2", "a2.fq,b2.fq", "-S", "o.sam"])
    );
}

#[test]
fn run_config_default_values() {
    let c = RunConfig::default();
    assert_eq!(c.mode, Mode::Quant);
    assert_eq!(c.quantifier, Quantifier::Kallisto);
    assert_eq!(c.strandedness, Strandedness::Unstranded);
    assert_eq!(c.reads_per_sample, 20_000_000);
    assert_eq!(c.fragment_len, None);
    assert_eq!(c.max_intron, 5000);
    assert_eq!(c.threads, 1);
}

// ---- R20: input validation ----

#[test]
fn fragment_length_nan_inf_zero_negative_rejected() {
    let smp = sample(&["a.fq"], &[]);
    let err = Err(ConfigError::InvalidFragmentLength { sample: "S1".into() });
    for bad in [f64::NAN, f64::INFINITY, 0.0, -5.0] {
        let mut c = cfg(Strandedness::Unstranded);
        c.fragment_len = Some((bad, 20.0));
        assert_eq!(kallisto_quant(&smp, "i.idx", "out", &c), err, "mean {bad}");
        c.fragment_len = Some((200.0, bad));
        assert_eq!(kallisto_quant(&smp, "i.idx", "out", &c), err, "sd {bad}");
    }
    // paired kallisto ignores fragment_len
    let mut c = cfg(Strandedness::Unstranded);
    c.fragment_len = Some((f64::NAN, -1.0));
    assert!(kallisto_quant(&sample(&["a_1.fq"], &["a_2.fq"]), "i.idx", "out", &c).is_ok());
}

#[test]
fn threads_zero_rejected_by_every_builder() {
    let mut c = cfg(Strandedness::Unstranded);
    c.threads = 0;
    c.fragment_len = Some((200.0, 20.0));
    let e = Err(ConfigError::InvalidThreads);
    assert_eq!(kallisto_quant(&sample(&["a_1.fq"], &["a_2.fq"]), "i", "o", &c), e);
    assert_eq!(kallisto_quant(&sample(&["a.fq"], &[]), "i", "o", &c), e);
    assert_eq!(hisat2_align(&sample(&["a_1.fq"], &["a_2.fq"]), "ix", "ss", "o.sam", &c), e);
    assert_eq!(hisat2_build("ref.fa", "idx/ref", &c), e);
    assert_eq!(feature_counts("a.sam", "g.saf", "c.txt", true, &c), e);
}

#[test]
fn kallisto_mate_count_mismatch_is_named_error() {
    let r = kallisto_quant(&sample(&["a1.fq", "b1.fq"], &["a2.fq"]), "i", "o", &RunConfig::default());
    assert_eq!(r, Err(ConfigError::MateCountMismatch { sample: "S1".into() }));
}

#[test]
fn hisat2_mate_count_mismatch_is_named_error() {
    let r = hisat2_align(&sample(&["a1.fq", "b1.fq"], &["a2.fq"]), "ix", "ss", "o.sam", &RunConfig::default());
    assert_eq!(r, Err(ConfigError::MateCountMismatch { sample: "S1".into() }));
}

#[test]
fn empty_fastq1_rejected_for_kallisto_and_hisat2() {
    let mut c = RunConfig::default();
    c.fragment_len = Some((200.0, 20.0));
    let e = Err(ConfigError::NoReads { sample: "S1".into() });
    assert_eq!(kallisto_quant(&sample(&[], &[]), "i", "o", &c), e);
    assert_eq!(hisat2_align(&sample(&[], &[]), "ix", "ss", "o.sam", &c), e);
    // empty fastq1 with non-empty fastq2 is also no usable R1 reads
    assert_eq!(kallisto_quant(&sample(&[], &["a2.fq"]), "i", "o", &c), e);
    assert_eq!(hisat2_align(&sample(&[], &["a2.fq"]), "ix", "ss", "o.sam", &c), e);
}

#[test]
fn hisat2_comma_in_path_rejected_kallisto_accepts_it() {
    let smp = sample(&["a,1.fq"], &["a2.fq"]);
    let r = hisat2_align(&smp, "ix", "ss", "o.sam", &RunConfig::default());
    assert_eq!(r, Err(ConfigError::CommaInPath { sample: "S1".into(), path: "a,1.fq".into() }));
    let smp2 = sample(&["a1.fq"], &["a,2.fq"]);
    let r = hisat2_align(&smp2, "ix", "ss", "o.sam", &RunConfig::default());
    assert_eq!(r, Err(ConfigError::CommaInPath { sample: "S1".into(), path: "a,2.fq".into() }));
    let k = kallisto_quant(&smp, "i.idx", "out", &RunConfig::default()).unwrap();
    assert_eq!(
        argv_native(&k),
        s(&["kallisto", "quant", "-i", "i.idx", "-o", "out", "--plaintext", "-t", "1", "a,1.fq", "a2.fq"])
    );
}

#[test]
fn config_error_display_names_sample() {
    let m = |e: ConfigError| e.to_string();
    assert!(m(ConfigError::FragmentLengthRequired { sample: "S1".into() }).contains("S1"));
    assert!(m(ConfigError::InvalidFragmentLength { sample: "S1".into() }).contains("S1"));
    assert!(m(ConfigError::MateCountMismatch { sample: "S1".into() }).contains("S1"));
    assert!(m(ConfigError::NoReads { sample: "S1".into() }).contains("S1"));
    let d = m(ConfigError::CommaInPath { sample: "S1".into(), path: "a,\"1\n.fq".into() });
    assert!(d.contains("S1"));
    assert!(d.contains("a,\\\"1\\n.fq"), "{d}");
    assert!(m(ConfigError::InvalidThreads).contains("threads"));
}
