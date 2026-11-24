//! The shared gene-by-sample matrix step used by the CLI and the browser page.

use std::collections::BTreeMap;

use qrnamap_core::matrix::{count_matrix, genes_from_tx2gene, MatrixError};

fn m(pairs: &[(&str, u64)]) -> BTreeMap<String, u64> {
    pairs.iter().map(|(g, c)| (g.to_string(), *c)).collect()
}

#[test]
fn genes_from_tx2gene_keeps_first_appearance_order() {
    let t = "t1\tgB\nt2\tgA\nt3\tgB\n\nt4\tgC\r\n";
    assert_eq!(genes_from_tx2gene(t).unwrap(), ["gB", "gA", "gC"]);
}

#[test]
fn genes_from_tx2gene_reports_the_bad_line() {
    let e = genes_from_tx2gene("t1\tg1\nonly_one_column\n").unwrap_err();
    assert_eq!(e.line, 2);
    let e = genes_from_tx2gene("t1\t\n").unwrap_err();
    assert_eq!(e.line, 1);
}

#[test]
fn matrix_with_gene_list_fixes_rows_and_warns() {
    let per = vec![
        ("S1".to_string(), m(&[("g1", 3), ("extra", 2)])),
        ("S2".to_string(), m(&[])),
    ];
    let genes = vec!["g2".to_string(), "g1".to_string()];
    let (mat, warnings) = count_matrix(&per, Some(&genes)).unwrap();
    assert_eq!(mat.to_tsv(), "gene\tS1\tS2\ng2\t0\t0\ng1\t3\t0\n");
    assert!(warnings.iter().any(|w| w.contains("extra")), "{warnings:?}");
    assert!(
        warnings
            .iter()
            .any(|w| w.contains("S2") && w.contains("all-zero")),
        "{warnings:?}"
    );
}

#[test]
fn matrix_without_gene_list_uses_sorted_union() {
    let per = vec![
        ("S1".to_string(), m(&[("b", 1)])),
        ("S2".to_string(), m(&[("a", 2)])),
    ];
    let (mat, _) = count_matrix(&per, None).unwrap();
    assert_eq!(mat.genes, ["a", "b"]);
}

#[test]
fn matrix_refuses_labels_with_tab_newline_cr() {
    for bad in ["x\ty", "x\ny", "x\ry"] {
        let per = vec![("S1".to_string(), m(&[(bad, 1)]))];
        assert_eq!(
            count_matrix(&per, None).unwrap_err(),
            MatrixError::BadGeneLabel(bad.to_string())
        );
        let genes = vec![bad.to_string()];
        assert!(count_matrix(&[], Some(&genes)).is_err());
    }
}
