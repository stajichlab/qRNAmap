//! Matrix assembly checks and the `qrnamap` binary's warnings. No external tool is needed.

use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

use qrnamap_cli::pipeline::{build_matrix, write_gene_counts, PipelineError};

fn m(pairs: &[(&str, u64)]) -> BTreeMap<String, u64> {
    pairs.iter().map(|(g, c)| (g.to_string(), *c)).collect()
}

#[test]
fn matrix_rejects_gene_label_with_tab() {
    for bad in ["gene\tX", "gene\nX", "gene\rX"] {
        let per = vec![("S1".to_string(), m(&[("geneA", 3), (bad, 1)]))];
        let err = build_matrix(&per, None).unwrap_err();
        assert!(matches!(err, PipelineError::BadGeneLabel(ref l) if l == bad), "{err:?}");
        // Also when the label comes only from the gene list.
        let per = vec![("S1".to_string(), m(&[("geneA", 3)]))];
        let genes = vec!["geneA".to_string(), bad.to_string()];
        assert!(matches!(build_matrix(&per, Some(&genes)), Err(PipelineError::BadGeneLabel(_))));
    }
    // The per-sample writer refuses it too.
    let dir = tempfile::tempdir().unwrap();
    let err = write_gene_counts(&dir.path().join("x.tsv"), "S1", &m(&[("a\tb", 1)])).unwrap_err();
    assert!(matches!(err, PipelineError::BadGeneLabel(_)));
}

#[test]
fn build_matrix_warns_unlisted_and_all_zero() {
    let per = vec![
        ("S1".to_string(), m(&[("geneA", 3), ("geneZ", 9)])),
        ("S2".to_string(), m(&[("geneA", 0)])),
        ("S3".to_string(), m(&[])),
    ];
    let genes = vec!["geneA".to_string(), "geneB".to_string()];
    let (mat, warnings) = build_matrix(&per, Some(&genes)).unwrap();
    assert_eq!(mat.genes, genes);
    assert_eq!(mat.counts, vec![vec![3, 0, 0], vec![0, 0, 0]]);
    let all = warnings.join("\n");
    assert!(all.contains("geneZ"), "{all}");
    assert!(all.contains("'S2'") && all.contains("all-zero sample"), "{all}");
    assert!(all.contains("'S3'"), "{all}");
    assert!(!all.contains("'S1': all-zero"), "{all}");
}

/// Write per-sample gene count files and a sheet; return the sheet path.
fn setup(dir: &Path, samples: &[(&str, BTreeMap<String, u64>)]) -> std::path::PathBuf {
    let mut sheet = String::from("sample\tfastq1\n");
    for (name, counts) in samples {
        write_gene_counts(&dir.join(format!("{name}.genecounts.tsv.zst")), name, counts).unwrap();
        sheet.push_str(&format!("{name}\t/x/{name}.fq\n"));
    }
    let p = dir.join("sheet.tsv");
    std::fs::write(&p, sheet).unwrap();
    p
}

fn qrnamap(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qrnamap")).args(args).output().unwrap()
}

#[test]
fn unlisted_genes_warning_is_emitted() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sheet = setup(d, &[("S1", m(&[("geneA", 5), ("geneQ", 2)])), ("S2", m(&[("geneA", 4)]))]);
    let t2g = d.join("tx2gene.tsv");
    std::fs::write(&t2g, "tA1\tgeneA\ntA2\tgeneA\ntB1\tgeneB\n").unwrap();
    let out = qrnamap(&["matrix", "--out", d.to_str().unwrap(), "--samples", sheet.to_str().unwrap(), "--genes-from", "tx2gene", t2g.to_str().unwrap()]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "{stderr}");
    assert!(stderr.contains("warning") && stderr.contains("geneQ"), "{stderr}");
    let text = qrnamap_cli::pipeline::read_table(&d.join("counts.tsv.zst")).unwrap();
    assert_eq!(text, "gene\tS1\tS2\ngeneA\t5\t4\ngeneB\t0\t0\n");
}

#[test]
fn all_zero_sample_warns() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sheet = setup(d, &[("S1", m(&[("geneA", 5)])), ("Zero", m(&[("geneA", 0)]))]);
    let out = qrnamap(&["matrix", "--out", d.to_str().unwrap(), "--samples", sheet.to_str().unwrap()]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "{stderr}");
    assert!(stderr.contains("all-zero sample") && stderr.contains("Zero"), "{stderr}");
}

#[test]
fn matrix_missing_sample_file_is_named_error() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sheet = setup(d, &[("S1", m(&[("geneA", 5)]))]);
    std::fs::write(&sheet, "sample\tfastq1\nS1\t/x/a.fq\nS2\t/x/b.fq\n").unwrap();
    let out = qrnamap(&["matrix", "--out", d.to_str().unwrap(), "--samples", sheet.to_str().unwrap()]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(stderr.contains("S2.genecounts.tsv.zst"), "{stderr}");
}

#[test]
fn transcripts_subcommand_writes_compressed_tables() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let fx = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/mini");
    let (fa, t2g, saf, ss) = (d.join("tx.fa.zst"), d.join("t2g.tsv.gz"), d.join("g.saf"), d.join("ss.txt.zst"));
    let out = qrnamap(&[
        "transcripts", "--reference", fx.join("genome.fa").to_str().unwrap(), "--gff", fx.join("genes.gff3").to_str().unwrap(),
        "--out", fa.to_str().unwrap(), "--tx2gene", t2g.to_str().unwrap(), "--saf", saf.to_str().unwrap(), "--splice-sites", ss.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    // Magic bytes: zstd, gzip, plain.
    assert_eq!(&std::fs::read(&fa).unwrap()[..4], &[0x28, 0xb5, 0x2f, 0xfd]);
    assert_eq!(&std::fs::read(&t2g).unwrap()[..2], &[0x1f, 0x8b]);
    assert!(std::fs::read_to_string(&saf).unwrap().starts_with("GeneID\t"));
    let t = qrnamap_cli::pipeline::read_table(&t2g).unwrap();
    assert!(t.starts_with("tA1\tgeneA\n"), "{t}");
    let f = qrnamap_cli::pipeline::read_table(&fa).unwrap();
    assert_eq!(f.matches('>').count(), 6);
    let s = qrnamap_cli::pipeline::read_table(&ss).unwrap();
    assert_eq!(s.lines().count(), 5);
}
