#![allow(non_snake_case)] // test names F/R/FR/RF are pinned by the brief

use qrnamap_core::{
    argv_native, feature_counts, hisat2_align, hisat2_build, kallisto_index, kallisto_quant,
    salmon_index, salmon_quant, ConfigError, Mode, Quantifier, RunConfig, Sample, Strandedness,
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
    RunConfig {
        strandedness: st,
        ..RunConfig::default()
    }
}

#[test]
fn kallisto_index_argv() {
    let c = kallisto_index("tx.fa", "tx.idx");
    assert_eq!(
        argv_native(&c),
        s(&["kallisto", "index", "-i", "tx.idx", "tx.fa"])
    );
    assert_eq!(c.inputs, p(&["tx.fa"]));
    assert_eq!(c.outputs, p(&["tx.idx"]));
}

#[test]
fn kallisto_paired_unstranded() {
    let c = kallisto_quant(
        &sample(&["a_1.fq"], &["a_2.fq"]),
        "i.idx",
        "out",
        &cfg(Strandedness::Unstranded),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "kallisto",
            "quant",
            "-i",
            "i.idx",
            "-o",
            "out",
            "--plaintext",
            "-t",
            "1",
            "a_1.fq",
            "a_2.fq"
        ])
    );
    assert_eq!(c.inputs, p(&["a_1.fq", "a_2.fq", "i.idx"]));
    assert_eq!(c.outputs, p(&["out/abundance.tsv", "out/run_info.json"]));
}

#[test]
fn kallisto_paired_reverse_adds_rf_stranded() {
    let c = kallisto_quant(
        &sample(&["a_1.fq"], &["a_2.fq"]),
        "i.idx",
        "out",
        &cfg(Strandedness::Reverse),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "kallisto",
            "quant",
            "-i",
            "i.idx",
            "-o",
            "out",
            "--plaintext",
            "-t",
            "1",
            "--rf-stranded",
            "a_1.fq",
            "a_2.fq"
        ])
    );
}

#[test]
fn kallisto_forward_adds_fr_stranded() {
    let c = kallisto_quant(
        &sample(&["a_1.fq"], &["a_2.fq"]),
        "i.idx",
        "out",
        &cfg(Strandedness::Forward),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "kallisto",
            "quant",
            "-i",
            "i.idx",
            "-o",
            "out",
            "--plaintext",
            "-t",
            "1",
            "--fr-stranded",
            "a_1.fq",
            "a_2.fq"
        ])
    );
}

#[test]
fn kallisto_single_end_needs_fragment_length() {
    let r = kallisto_quant(
        &sample(&["a.fq"], &[]),
        "i.idx",
        "out",
        &RunConfig::default(),
    );
    assert_eq!(
        r,
        Err(ConfigError::FragmentLengthRequired {
            sample: "S1".into()
        })
    );
}

#[test]
fn kallisto_single_end_uses_l_and_s() {
    let mut c0 = cfg(Strandedness::Unstranded);
    c0.fragment_len = Some((200.0, 20.0));
    let c = kallisto_quant(&sample(&["a.fq"], &[]), "i.idx", "out", &c0).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "kallisto",
            "quant",
            "-i",
            "i.idx",
            "-o",
            "out",
            "--plaintext",
            "-t",
            "1",
            "--single",
            "-l",
            "200",
            "-s",
            "20",
            "a.fq"
        ])
    );
    c0.fragment_len = Some((187.5, 24.25));
    let c = kallisto_quant(&sample(&["a.fq"], &[]), "i.idx", "out", &c0).unwrap();
    let a = argv_native(&c);
    assert_eq!(&a[10..14], &s(&["-l", "187.5", "-s", "24.25"])[..]);
}

#[test]
fn kallisto_single_end_stranded_exact_argv() {
    // R29-10: single-end stranded kallisto, the exact full argv for each strandedness.
    for (st, flag) in [
        (Strandedness::Forward, "--fr-stranded"),
        (Strandedness::Reverse, "--rf-stranded"),
    ] {
        let c0 = RunConfig {
            strandedness: st,
            fragment_len: Some((180.0, 25.0)),
            ..RunConfig::default()
        };
        let c = kallisto_quant(&sample(&["a.fq", "b.fq"], &[]), "i.idx", "out", &c0).unwrap();
        assert_eq!(
            argv_native(&c),
            s(&[
                "kallisto",
                "quant",
                "-i",
                "i.idx",
                "-o",
                "out",
                "--plaintext",
                "-t",
                "1",
                flag,
                "--single",
                "-l",
                "180",
                "-s",
                "25",
                "a.fq",
                "b.fq"
            ])
        );
        assert_eq!(c.inputs, p(&["a.fq", "b.fq", "i.idx"]));
    }
}

#[test]
fn hisat2_pair_reverse_is_RF() {
    let c = hisat2_align(
        &sample(&["a_1.fq"], &["a_2.fq"]),
        "ix",
        "ss.txt",
        "o.sam",
        &cfg(Strandedness::Reverse),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "hisat2-align-s",
            "-p",
            "1",
            "-x",
            "ix",
            "--known-splicesite-infile",
            "ss.txt",
            "--no-temp-splicesite",
            "--max-intronlen",
            "5000",
            "--rna-strandness",
            "RF",
            "-1",
            "a_1.fq",
            "-2",
            "a_2.fq",
            "-S",
            "o.sam"
        ])
    );
    assert_eq!(c.inputs, p(&["a_1.fq", "a_2.fq", "ss.txt", "ix"]));
    assert_eq!(c.outputs, p(&["o.sam"]));
}

#[test]
fn hisat2_pair_forward_is_FR() {
    let c = hisat2_align(
        &sample(&["a_1.fq"], &["a_2.fq"]),
        "ix",
        "ss.txt",
        "o.sam",
        &cfg(Strandedness::Forward),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "hisat2-align-s",
            "-p",
            "1",
            "-x",
            "ix",
            "--known-splicesite-infile",
            "ss.txt",
            "--no-temp-splicesite",
            "--max-intronlen",
            "5000",
            "--rna-strandness",
            "FR",
            "-1",
            "a_1.fq",
            "-2",
            "a_2.fq",
            "-S",
            "o.sam"
        ])
    );
}

#[test]
fn hisat2_single_forward_is_F() {
    let c = hisat2_align(
        &sample(&["a.fq"], &[]),
        "ix",
        "ss.txt",
        "o.sam",
        &cfg(Strandedness::Forward),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "hisat2-align-s",
            "-p",
            "1",
            "-x",
            "ix",
            "--known-splicesite-infile",
            "ss.txt",
            "--no-temp-splicesite",
            "--max-intronlen",
            "5000",
            "--rna-strandness",
            "F",
            "-U",
            "a.fq",
            "-S",
            "o.sam"
        ])
    );
}

#[test]
fn hisat2_single_reverse_is_R() {
    let c = hisat2_align(
        &sample(&["a.fq"], &[]),
        "ix",
        "ss.txt",
        "o.sam",
        &cfg(Strandedness::Reverse),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "hisat2-align-s",
            "-p",
            "1",
            "-x",
            "ix",
            "--known-splicesite-infile",
            "ss.txt",
            "--no-temp-splicesite",
            "--max-intronlen",
            "5000",
            "--rna-strandness",
            "R",
            "-U",
            "a.fq",
            "-S",
            "o.sam"
        ])
    );
}

#[test]
fn hisat2_has_no_temp_splicesite_and_known_sites() {
    let mut c0 = cfg(Strandedness::Unstranded);
    c0.threads = 4;
    c0.max_intron = 3000;
    let c = hisat2_align(
        &sample(&["a_1.fq"], &["a_2.fq"]),
        "ix",
        "ss.txt",
        "o.sam",
        &c0,
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "hisat2-align-s",
            "-p",
            "4",
            "-x",
            "ix",
            "--known-splicesite-infile",
            "ss.txt",
            "--no-temp-splicesite",
            "--max-intronlen",
            "3000",
            "-1",
            "a_1.fq",
            "-2",
            "a_2.fq",
            "-S",
            "o.sam"
        ])
    );
}

#[test]
fn hisat2_build_outputs_eight_files() {
    let c = hisat2_build("ref.fa", "idx/ref", &RunConfig::default()).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&["hisat2-build", "-p", "1", "ref.fa", "idx/ref"])
    );
    assert_eq!(c.inputs, p(&["ref.fa"]));
    let want: Vec<PathBuf> = (1..=8)
        .map(|i| PathBuf::from(format!("idx/ref.{i}.ht2")))
        .collect();
    assert_eq!(c.outputs, want);
    assert_eq!(c.outputs.len(), 8);
}

#[test]
fn featurecounts_paired_has_p_and_countreadpairs() {
    let c = feature_counts(
        "a.sam",
        "g.saf",
        "counts.txt",
        true,
        &cfg(Strandedness::Unstranded),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "featureCounts",
            "-T",
            "1",
            "-F",
            "SAF",
            "-a",
            "g.saf",
            "-p",
            "--countReadPairs",
            "-s",
            "0",
            "-o",
            "counts.txt",
            "a.sam"
        ])
    );
    assert_eq!(c.inputs, p(&["a.sam", "g.saf"]));
}

#[test]
fn featurecounts_single_end_has_no_p() {
    let c = feature_counts(
        "a.sam",
        "g.saf",
        "counts.txt",
        false,
        &cfg(Strandedness::Forward),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "featureCounts",
            "-T",
            "1",
            "-F",
            "SAF",
            "-a",
            "g.saf",
            "-s",
            "1",
            "-o",
            "counts.txt",
            "a.sam"
        ])
    );
}

#[test]
fn featurecounts_reverse_is_s2() {
    let c = feature_counts(
        "a.sam",
        "g.saf",
        "counts.txt",
        true,
        &cfg(Strandedness::Reverse),
    )
    .unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "featureCounts",
            "-T",
            "1",
            "-F",
            "SAF",
            "-a",
            "g.saf",
            "-p",
            "--countReadPairs",
            "-s",
            "2",
            "-o",
            "counts.txt",
            "a.sam"
        ])
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
        s(&[
            "kallisto",
            "quant",
            "-i",
            "i.idx",
            "-o",
            "out",
            "--plaintext",
            "-t",
            "1",
            "a1.fq",
            "a2.fq",
            "b1.fq",
            "b2.fq"
        ])
    );
    let h = hisat2_align(&smp, "ix", "ss.txt", "o.sam", &RunConfig::default()).unwrap();
    assert_eq!(
        argv_native(&h),
        s(&[
            "hisat2-align-s",
            "-p",
            "1",
            "-x",
            "ix",
            "--known-splicesite-infile",
            "ss.txt",
            "--no-temp-splicesite",
            "--max-intronlen",
            "5000",
            "-1",
            "a1.fq,b1.fq",
            "-2",
            "a2.fq,b2.fq",
            "-S",
            "o.sam"
        ])
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
    let err = Err(ConfigError::InvalidFragmentLength {
        sample: "S1".into(),
    });
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
    assert_eq!(
        kallisto_quant(&sample(&["a_1.fq"], &["a_2.fq"]), "i", "o", &c),
        e
    );
    assert_eq!(kallisto_quant(&sample(&["a.fq"], &[]), "i", "o", &c), e);
    assert_eq!(
        hisat2_align(&sample(&["a_1.fq"], &["a_2.fq"]), "ix", "ss", "o.sam", &c),
        e
    );
    assert_eq!(hisat2_build("ref.fa", "idx/ref", &c), e);
    assert_eq!(feature_counts("a.sam", "g.saf", "c.txt", true, &c), e);
}

#[test]
fn kallisto_mate_count_mismatch_is_named_error() {
    let r = kallisto_quant(
        &sample(&["a1.fq", "b1.fq"], &["a2.fq"]),
        "i",
        "o",
        &RunConfig::default(),
    );
    assert_eq!(
        r,
        Err(ConfigError::MateCountMismatch {
            sample: "S1".into()
        })
    );
}

#[test]
fn hisat2_mate_count_mismatch_is_named_error() {
    let r = hisat2_align(
        &sample(&["a1.fq", "b1.fq"], &["a2.fq"]),
        "ix",
        "ss",
        "o.sam",
        &RunConfig::default(),
    );
    assert_eq!(
        r,
        Err(ConfigError::MateCountMismatch {
            sample: "S1".into()
        })
    );
}

#[test]
fn empty_fastq1_rejected_for_kallisto_and_hisat2() {
    let c = RunConfig {
        fragment_len: Some((200.0, 20.0)),
        ..RunConfig::default()
    };
    let e = Err(ConfigError::NoReads {
        sample: "S1".into(),
    });
    assert_eq!(kallisto_quant(&sample(&[], &[]), "i", "o", &c), e);
    assert_eq!(hisat2_align(&sample(&[], &[]), "ix", "ss", "o.sam", &c), e);
    // empty fastq1 with non-empty fastq2 is also no usable R1 reads
    assert_eq!(kallisto_quant(&sample(&[], &["a2.fq"]), "i", "o", &c), e);
    assert_eq!(
        hisat2_align(&sample(&[], &["a2.fq"]), "ix", "ss", "o.sam", &c),
        e
    );
}

#[test]
fn hisat2_comma_in_path_rejected_kallisto_accepts_it() {
    let smp = sample(&["a,1.fq"], &["a2.fq"]);
    let r = hisat2_align(&smp, "ix", "ss", "o.sam", &RunConfig::default());
    assert_eq!(
        r,
        Err(ConfigError::CommaInPath {
            sample: "S1".into(),
            path: "a,1.fq".into()
        })
    );
    let smp2 = sample(&["a1.fq"], &["a,2.fq"]);
    let r = hisat2_align(&smp2, "ix", "ss", "o.sam", &RunConfig::default());
    assert_eq!(
        r,
        Err(ConfigError::CommaInPath {
            sample: "S1".into(),
            path: "a,2.fq".into()
        })
    );
    let k = kallisto_quant(&smp, "i.idx", "out", &RunConfig::default()).unwrap();
    assert_eq!(
        argv_native(&k),
        s(&[
            "kallisto",
            "quant",
            "-i",
            "i.idx",
            "-o",
            "out",
            "--plaintext",
            "-t",
            "1",
            "a,1.fq",
            "a2.fq"
        ])
    );
}

#[test]
fn config_error_display_names_sample() {
    let m = |e: ConfigError| e.to_string();
    assert!(m(ConfigError::FragmentLengthRequired {
        sample: "S1".into()
    })
    .contains("S1"));
    assert!(m(ConfigError::InvalidFragmentLength {
        sample: "S1".into()
    })
    .contains("S1"));
    assert!(m(ConfigError::MateCountMismatch {
        sample: "S1".into()
    })
    .contains("S1"));
    assert!(m(ConfigError::NoReads {
        sample: "S1".into()
    })
    .contains("S1"));
    let d = m(ConfigError::CommaInPath {
        sample: "S1".into(),
        path: "a,\"1\n.fq".into(),
    });
    assert!(d.contains("S1"));
    assert!(d.contains("a,\\\"1\\n.fq"), "{d}");
    assert!(m(ConfigError::InvalidThreads).contains("threads"));
}

fn threads_cfg(st: Strandedness, threads: u32) -> RunConfig {
    RunConfig {
        strandedness: st,
        threads,
        ..RunConfig::default()
    }
}

#[test]
fn salmon_quant_paired_unstranded() {
    let c = salmon_quant(
        &sample(&["a"], &["b"]),
        "idx",
        "out",
        &cfg(Strandedness::Unstranded),
    )
    .unwrap();
    let a = argv_native(&c);
    assert_eq!(a[0], "salmon");
    let joined = a.join(" ");
    assert!(joined.contains("-l IU -1 a -2 b"), "{joined}");
    assert!(!joined.contains("--validateMappings"), "{joined}");
    assert!(joined.contains("-p 1"), "{joined}");
    assert!(joined.ends_with("-o out"), "{joined}");
    assert_eq!(&a[1], "quant");
    assert_eq!(
        c.outputs,
        p(&["out/quant.sf", "out/aux_info/meta_info.json"])
    );
}

#[test]
fn salmon_quant_two_lanes_paired() {
    let c = salmon_quant(
        &sample(&["l1_1.fq", "l2_1.fq"], &["l1_2.fq", "l2_2.fq"]),
        "idx",
        "out",
        &cfg(Strandedness::Unstranded),
    )
    .unwrap();
    let joined = argv_native(&c).join(" ");
    assert!(
        joined.contains("-1 l1_1.fq l2_1.fq -2 l1_2.fq l2_2.fq"),
        "{joined}"
    );
}

#[test]
fn salmon_quant_strandedness_paired() {
    let f = salmon_quant(
        &sample(&["a"], &["b"]),
        "idx",
        "out",
        &cfg(Strandedness::Forward),
    )
    .unwrap();
    assert!(argv_native(&f).join(" ").contains("-l ISF "));
    let r = salmon_quant(
        &sample(&["a"], &["b"]),
        "idx",
        "out",
        &cfg(Strandedness::Reverse),
    )
    .unwrap();
    assert!(argv_native(&r).join(" ").contains("-l ISR "));
}

#[test]
fn salmon_quant_single_end() {
    let mut c0 = cfg(Strandedness::Unstranded);
    c0.fragment_len = Some((200.0, 20.0));
    let c = salmon_quant(&sample(&["a.fq"], &[]), "idx", "out", &c0).unwrap();
    let joined = argv_native(&c).join(" ");
    assert!(joined.contains("-l U -r a.fq"), "{joined}");
    assert!(joined.contains("--fldMean 200 --fldSD 20"), "{joined}");
    let mut c1 = cfg(Strandedness::Forward);
    c1.fragment_len = Some((200.0, 20.0));
    let f = salmon_quant(&sample(&["a.fq"], &[]), "idx", "out", &c1).unwrap();
    assert!(argv_native(&f).join(" ").contains("-l SF "));
    c1.strandedness = Strandedness::Reverse;
    let r = salmon_quant(&sample(&["a.fq"], &[]), "idx", "out", &c1).unwrap();
    assert!(argv_native(&r).join(" ").contains("-l SR "));
}

#[test]
fn salmon_quant_single_end_needs_fragment_length() {
    let e = salmon_quant(
        &sample(&["a.fq"], &[]),
        "idx",
        "out",
        &cfg(Strandedness::Unstranded),
    )
    .unwrap_err();
    assert_eq!(
        e,
        ConfigError::FragmentLengthRequired {
            sample: "S1".into()
        }
    );
    for bad in [(f64::NAN, 20.0), (200.0, 0.0), (0.0, 20.0)] {
        let mut c0 = cfg(Strandedness::Unstranded);
        c0.fragment_len = Some(bad);
        let e = salmon_quant(&sample(&["a.fq"], &[]), "idx", "out", &c0).unwrap_err();
        assert_eq!(
            e,
            ConfigError::InvalidFragmentLength {
                sample: "S1".into()
            }
        );
    }
}

#[test]
fn salmon_quant_rejects_bad_input() {
    let ok = sample(&["a"], &["b"]);
    assert_eq!(
        salmon_quant(&ok, "idx", "out", &threads_cfg(Strandedness::Unstranded, 0)).unwrap_err(),
        ConfigError::InvalidThreads
    );
    assert_eq!(
        salmon_quant(
            &sample(&[], &[]),
            "idx",
            "out",
            &cfg(Strandedness::Unstranded)
        )
        .unwrap_err(),
        ConfigError::NoReads {
            sample: "S1".into()
        }
    );
    assert_eq!(
        salmon_quant(
            &sample(&["a", "c"], &["b"]),
            "idx",
            "out",
            &cfg(Strandedness::Unstranded)
        )
        .unwrap_err(),
        ConfigError::MateCountMismatch {
            sample: "S1".into()
        }
    );
}

#[test]
fn salmon_index_args() {
    let c = salmon_index("t.fa", "idx", &threads_cfg(Strandedness::Unstranded, 4)).unwrap();
    assert_eq!(
        argv_native(&c),
        s(&[
            "salmon",
            "index",
            "-t",
            "t.fa",
            "-i",
            "idx",
            "-p",
            "4",
            "--keepDuplicates"
        ])
    );
    assert_eq!(c.inputs, p(&["t.fa"]));
    assert_eq!(c.outputs, p(&["idx"]));
    assert_eq!(
        salmon_index("t.fa", "idx", &threads_cfg(Strandedness::Unstranded, 0)).unwrap_err(),
        ConfigError::InvalidThreads
    );
}

#[test]
fn quantifier_name() {
    assert_eq!(Quantifier::Kallisto.name(), "kallisto");
    assert_eq!(Quantifier::Salmon.name(), "salmon");
}

fn salmon_opt_cfg(softclip: bool, frac: Option<f64>) -> RunConfig {
    RunConfig {
        salmon_softclip_overhangs: softclip,
        salmon_min_score_fraction: frac,
        ..RunConfig::default()
    }
}

fn salmon_argv(c: &RunConfig) -> Result<Vec<String>, ConfigError> {
    salmon_quant(&sample(&["a"], &["b"]), "idx", "out", c).map(|t| argv_native(&t))
}

#[test]
fn salmon_quant_default_argv_unchanged() {
    // Same argv as main (fe1ce5e) for a paired unstranded sample.
    assert_eq!(
        salmon_argv(&cfg(Strandedness::Unstranded)).unwrap(),
        s(&[
            "salmon", "quant", "-i", "idx", "-l", "IU", "-1", "a", "-2", "b", "-p", "1", "-o",
            "out"
        ])
    );
}

#[test]
fn salmon_quant_softclip_overhangs_flag() {
    assert_eq!(
        salmon_argv(&salmon_opt_cfg(true, None)).unwrap(),
        s(&[
            "salmon",
            "quant",
            "-i",
            "idx",
            "-l",
            "IU",
            "-1",
            "a",
            "-2",
            "b",
            "--softclipOverhangs",
            "-p",
            "1",
            "-o",
            "out"
        ])
    );
}

#[test]
fn salmon_quant_min_score_fraction_flag() {
    let a = salmon_argv(&salmon_opt_cfg(false, Some(0.5))).unwrap();
    let i = a.iter().position(|x| x == "--minScoreFraction").unwrap();
    assert_eq!(a[i + 1], "0.5");
    assert!(!a.contains(&"--softclipOverhangs".to_string()));
    let a = salmon_argv(&salmon_opt_cfg(false, Some(1.0))).unwrap();
    let i = a.iter().position(|x| x == "--minScoreFraction").unwrap();
    assert_eq!(a[i + 1], "1");
}

#[test]
fn salmon_quant_both_flags_order() {
    let a = salmon_argv(&salmon_opt_cfg(true, Some(0.5))).unwrap();
    let pos = |k: &str| a.iter().position(|x| x == k).unwrap();
    assert!(pos("--softclipOverhangs") < pos("--minScoreFraction"));
    assert!(pos("--minScoreFraction") < pos("-p"));
}

#[test]
fn salmon_quant_min_score_fraction_rejects_bad_values() {
    for v in [0.0, -0.1, 1.5, f64::NAN, f64::INFINITY] {
        let e = salmon_argv(&salmon_opt_cfg(false, Some(v))).unwrap_err();
        assert!(
            matches!(e, ConfigError::InvalidMinScoreFraction { .. }),
            "{v}: {e:?}"
        );
        assert!(e.to_string().contains("--salmon-min-score-fraction"));
    }
}

#[test]
fn kallisto_quant_ignores_salmon_options() {
    let smp = sample(&["a"], &["b"]);
    let base = kallisto_quant(&smp, "i.idx", "out", &RunConfig::default()).unwrap();
    let set = kallisto_quant(&smp, "i.idx", "out", &salmon_opt_cfg(true, Some(0.5))).unwrap();
    assert_eq!(argv_native(&base), argv_native(&set));
}

#[test]
fn salmon_quant_single_end_with_overhang_flags() {
    let mut c = salmon_opt_cfg(true, Some(0.5));
    c.fragment_len = Some((200.0, 20.0));
    let t = salmon_quant(&sample(&["a.fq"], &[]), "idx", "out", &c).unwrap();
    assert_eq!(
        argv_native(&t),
        s(&[
            "salmon",
            "quant",
            "-i",
            "idx",
            "-l",
            "U",
            "-r",
            "a.fq",
            "--fldMean",
            "200",
            "--fldSD",
            "20",
            "--softclipOverhangs",
            "--minScoreFraction",
            "0.5",
            "-p",
            "1",
            "-o",
            "out"
        ])
    );
}
