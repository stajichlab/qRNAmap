//! QC numbers on a count matrix: log-CPM, sample correlation, PCA.

use qrnamap_core::qc::{log_cpm, pca, sample_correlation, QcError};
use qrnamap_core::CountMatrix;

fn matrix(samples: &[&str], rows: &[(&str, &[u64])]) -> CountMatrix {
    CountMatrix {
        genes: rows.iter().map(|(g, _)| g.to_string()).collect(),
        samples: samples.iter().map(|s| s.to_string()).collect(),
        counts: rows.iter().map(|(_, c)| c.to_vec()).collect(),
    }
}

/// Library sizes are 1e6 for every sample, so CPM equals the count for the small genes.
fn small() -> CountMatrix {
    matrix(
        &["a", "b", "c"],
        &[
            ("g1", &[1, 7, 0]),
            ("g2", &[3, 0, 0]),
            ("g3", &[999_996, 999_993, 1_000_000]),
        ],
    )
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn log_cpm_values() {
    let l = log_cpm(&small()).unwrap();
    // cpm 1 -> log2(2) = 1; cpm 3 -> log2(4) = 2; cpm 7 -> log2(8) = 3; cpm 0 -> 0.
    assert!(close(l[0][0], 1.0));
    assert!(close(l[0][1], 3.0));
    assert!(close(l[0][2], 0.0));
    assert!(close(l[1][0], 2.0));
    assert!(close(l[1][1], 0.0));
    assert!(close(l[2][0], 999_997f64.log2()));
    assert!(close(l[2][1], 999_994f64.log2()));
    assert!(close(l[2][2], 1_000_001f64.log2()));
    assert_eq!(l.len(), 3);
    assert_eq!(l[0].len(), 3);
}

#[test]
fn log_cpm_zero_library_names_the_sample() {
    let m = matrix(&["a", "empty"], &[("g1", &[5, 0]), ("g2", &[5, 0])]);
    let e = log_cpm(&m).unwrap_err();
    assert_eq!(e, QcError::ZeroLibrary("empty".into()));
    assert!(e.to_string().contains("empty"));
}

#[test]
fn correlation_diagonal_one_and_symmetric() {
    let c = sample_correlation(&small()).unwrap();
    assert_eq!(c.len(), 3);
    for i in 0..3 {
        assert!(close(c[i][i], 1.0));
        for j in 0..3 {
            assert!(close(c[i][j], c[j][i]));
            assert!(c[i][j] <= 1.0 + 1e-12 && c[i][j] >= -1.0 - 1e-12);
        }
    }
}

#[test]
fn correlation_of_identical_samples_is_one() {
    let m = matrix(
        &["a", "b", "c"],
        &[
            ("g1", &[10, 10, 40]),
            ("g2", &[200, 200, 10]),
            ("g3", &[5, 5, 90]),
            ("g4", &[0, 0, 0]),
        ],
    );
    let c = sample_correlation(&m).unwrap();
    assert!(close(c[0][1], 1.0));
    assert!(c[0][2] < 1.0);
}

#[test]
fn correlation_ignores_genes_with_zero_total() {
    let with_zero = matrix(
        &["a", "b"],
        &[("g1", &[10, 30]), ("g2", &[50, 20]), ("g3", &[0, 0])],
    );
    let without = matrix(&["a", "b"], &[("g1", &[10, 30]), ("g2", &[50, 20])]);
    let x = sample_correlation(&with_zero).unwrap();
    let y = sample_correlation(&without).unwrap();
    assert!(close(x[0][1], y[0][1]));
}

#[test]
fn correlation_needs_two_samples() {
    let m = matrix(&["a"], &[("g1", &[5]), ("g2", &[6])]);
    assert_eq!(
        sample_correlation(&m).unwrap_err(),
        QcError::TooFewSamples { need: 2, have: 1 }
    );
}

fn two_groups() -> CountMatrix {
    // Samples a1, a2 are high in g1 and g2; b1, b2 are high in g3 and g4.
    matrix(
        &["a1", "a2", "b1", "b2"],
        &[
            ("g1", &[900, 880, 100, 120]),
            ("g2", &[500, 520, 50, 40]),
            ("g3", &[60, 70, 800, 790]),
            ("g4", &[30, 35, 600, 610]),
            ("g5", &[300, 300, 300, 300]),
        ],
    )
}

#[test]
fn pca_separates_two_groups_on_pc1() {
    let p = pca(&two_groups(), 500).unwrap();
    assert_eq!(p.scores.len(), 4);
    let pc1: Vec<f64> = p.scores.iter().map(|s| s[0]).collect();
    assert!(pc1[0].signum() == pc1[1].signum());
    assert!(pc1[2].signum() == pc1[3].signum());
    assert!(pc1[0].signum() != pc1[2].signum());
    assert!(p.variance_explained[0] > 0.8, "{:?}", p.variance_explained);
    assert!(p.variance_explained[0] >= p.variance_explained[1]);
    assert!(p.variance_explained[0] + p.variance_explained[1] <= 1.0 + 1e-9);
}

#[test]
fn pca_scores_are_centred() {
    let p = pca(&two_groups(), 500).unwrap();
    for k in 0..2 {
        let sum: f64 = p.scores.iter().map(|s| s[k]).sum();
        assert!(sum.abs() < 1e-9, "PC{} sum {sum}", k + 1);
    }
}

#[test]
fn pca_top_genes_limit_is_used() {
    let all = pca(&two_groups(), 500).unwrap();
    let two = pca(&two_groups(), 2).unwrap();
    assert_ne!(all.scores, two.scores);
    assert_eq!(two.scores.len(), 4);
}

#[test]
fn pca_needs_three_samples() {
    let m = matrix(&["a", "b"], &[("g1", &[5, 9]), ("g2", &[6, 1])]);
    assert_eq!(
        pca(&m, 500).unwrap_err(),
        QcError::TooFewSamples { need: 3, have: 2 }
    );
}

#[test]
fn pca_zero_library_names_the_sample() {
    let m = matrix(
        &["a", "b", "dead"],
        &[("g1", &[5, 9, 0]), ("g2", &[6, 1, 0])],
    );
    assert_eq!(
        pca(&m, 500).unwrap_err(),
        QcError::ZeroLibrary("dead".into())
    );
}
