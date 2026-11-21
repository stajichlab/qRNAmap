//! Count parser and matrix tests. Expected values are literals written by hand.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use qrnamap_core::{
    assemble, build_transcripts, gene_counts_from_abundance, gene_counts_from_salmon, parse_featurecounts, saf, tx2gene_tsv,
    unlisted_genes, CountError, CountMatrix,
};
use qseq_io::{parse_gff3, read_fasta};

const GENOME_FA: &str = "../../tests/fixtures/mini/genome.fa";
const GENES_GFF3: &str = "../../tests/fixtures/mini/genes.gff3";

fn t2g(pairs: &[(&str, &str)]) -> HashMap<String, String> {
    pairs
        .iter()
        .map(|(t, g)| (t.to_string(), g.to_string()))
        .collect()
}

fn abundance(rows: &[(&str, &str)]) -> String {
    let mut s = String::from("target_id\tlength\teff_length\test_counts\ttpm\n");
    for (id, est) in rows {
        s.push_str(&format!("{id}\t100\t50.0\t{est}\t1.0\n"));
    }
    s
}

fn map(pairs: &[(&str, u64)]) -> BTreeMap<String, u64> {
    pairs.iter().map(|(g, c)| (g.to_string(), *c)).collect()
}

fn names(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn sum_est_counts_per_gene_rounds() {
    let m = t2g(&[("t1", "g1"), ("t2", "g1"), ("t3", "g2"), ("t4", "g3")]);
    let tsv = abundance(&[("t1", "1.4"), ("t2", "1.2"), ("t3", "2.5"), ("t4", "2.4999")]);
    let got = gene_counts_from_abundance(&tsv, &m).unwrap();
    // g1: 1.4 + 1.2 = 2.6 -> 3 (rounding each transcript first would give 1 + 1 = 2)
    assert_eq!(got, map(&[("g1", 3), ("g2", 3), ("g3", 2)]));
}

#[test]
fn half_rounds_away_from_zero() {
    let m = t2g(&[("a", "ga"), ("b", "gb"), ("c", "gc")]);
    let tsv = abundance(&[("a", "0.5"), ("b", "1.5"), ("c", "2.5")]);
    let got = gene_counts_from_abundance(&tsv, &m).unwrap();
    assert_eq!(got, map(&[("ga", 1), ("gb", 2), ("gc", 3)]));
}

#[test]
fn unknown_transcript_is_error() {
    let m = t2g(&[("t1", "g1")]);
    let tsv = abundance(&[("t1", "1"), ("tX", "2")]);
    assert_eq!(
        gene_counts_from_abundance(&tsv, &m),
        Err(CountError::UnknownTranscript("tX".into()))
    );
}

#[test]
fn zero_count_gene_present_in_output() {
    let m = t2g(&[("t1", "g1"), ("t2", "g2"), ("t3", "g3")]);
    let tsv = abundance(&[("t1", "0"), ("t2", "4")]);
    let got = gene_counts_from_abundance(&tsv, &m).unwrap();
    // g1 is listed with 0 and appears; g3 is absent from the table and gets no entry.
    assert_eq!(got, map(&[("g1", 0), ("g2", 4)]));
}

#[test]
fn abundance_missing_column_is_named_error() {
    let m = t2g(&[("t1", "g1")]);
    let tsv = "target_id\tlength\teff_length\ttpm\nt1\t100\t50\t1\n";
    assert_eq!(
        gene_counts_from_abundance(tsv, &m),
        Err(CountError::MissingColumn("est_counts".into()))
    );
}

#[test]
fn abundance_columns_found_by_name() {
    let m = t2g(&[("t1", "g1")]);
    let tsv = "tpm\test_counts\ttarget_id\n9\t7\tt1\n";
    assert_eq!(gene_counts_from_abundance(tsv, &m).unwrap(), map(&[("g1", 7)]));
}

#[test]
fn abundance_negative_est_counts_is_bad_row_with_line() {
    let m = t2g(&[("t1", "g1"), ("t2", "g2")]);
    let tsv = abundance(&[("t1", "1"), ("t2", "-0.5")]);
    assert_eq!(
        gene_counts_from_abundance(&tsv, &m),
        Err(CountError::BadRow { line: 3 })
    );
}

#[test]
fn abundance_non_finite_and_short_rows_are_bad_rows() {
    let m = t2g(&[("t1", "g1")]);
    assert_eq!(
        gene_counts_from_abundance(&abundance(&[("t1", "NaN")]), &m),
        Err(CountError::BadRow { line: 2 })
    );
    assert_eq!(
        gene_counts_from_abundance(&abundance(&[("t1", "inf")]), &m),
        Err(CountError::BadRow { line: 2 })
    );
    let short = "target_id\tlength\teff_length\test_counts\ttpm\nt1\t100\n";
    assert_eq!(
        gene_counts_from_abundance(short, &m),
        Err(CountError::BadRow { line: 2 })
    );
}

#[test]
fn abundance_crlf_and_blank_trailing_line() {
    let m = t2g(&[("t1", "g1")]);
    let tsv = "target_id\tlength\teff_length\test_counts\ttpm\r\nt1\t100\t50\t3.2\t1\r\n\r\n";
    assert_eq!(gene_counts_from_abundance(tsv, &m).unwrap(), map(&[("g1", 3)]));
}

const FC: &str = "# Program:featureCounts v2.1.1; Command:\"featureCounts\" \"-F\" \"SAF\"\n\
Geneid\tChr\tStart\tEnd\tStrand\tLength\t/path/s1.sam\n\
geneA\tchrA;chrA\t101;301\t160;350\t+;+\t110\t12\n\
geneB\tchrA\t601\t870\t-\t270\t0\n";

#[test]
fn featurecounts_skips_comment_and_header() {
    let got = parse_featurecounts(FC).unwrap();
    assert_eq!(got.len(), 2);
    assert!(!got.contains_key("Geneid"));
}

#[test]
fn featurecounts_reads_last_column() {
    let got = parse_featurecounts(FC).unwrap();
    assert_eq!(got, map(&[("geneA", 12), ("geneB", 0)]));
}

#[test]
fn featurecounts_crlf_and_blank_trailing_line() {
    let tsv = "# Program:featureCounts\r\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ts.bam\r\ng1\tc\t1\t2\t+\t2\t5\r\n\r\n";
    assert_eq!(parse_featurecounts(tsv).unwrap(), map(&[("g1", 5)]));
}

#[test]
fn featurecounts_non_integer_count_is_bad_row_with_line() {
    let tsv = "# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ts.bam\ng1\tc\t1\t2\t+\t2\t5\ng2\tc\t1\t2\t+\t2\t1.5\n";
    assert_eq!(parse_featurecounts(tsv), Err(CountError::BadRow { line: 4 }));
    let neg = "# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ts.bam\ng1\tc\t1\t2\t+\t2\t-3\n";
    assert_eq!(parse_featurecounts(neg), Err(CountError::BadRow { line: 3 }));
}

const FC_HDR: &str = "Geneid\tChr\tStart\tEnd\tStrand\tLength\ts.bam\n";

#[test]
fn featurecounts_without_comment_line_keeps_first_gene() {
    let tsv = format!("{FC_HDR}g1\tc\t1\t2\t+\t2\t5\ng2\tc\t1\t2\t+\t2\t0\ng3\tc\t1\t2\t+\t2\t9\n");
    assert_eq!(
        parse_featurecounts(&tsv).unwrap(),
        map(&[("g1", 5), ("g2", 0), ("g3", 9)])
    );
}

#[test]
fn featurecounts_first_cell_not_geneid_is_missing_header() {
    let tsv = "g1\tc\t1\t2\t+\t2\t5\ng2\tc\t1\t2\t+\t2\t7\n";
    assert_eq!(parse_featurecounts(tsv), Err(CountError::MissingHeader));
}

#[test]
fn featurecounts_empty_and_comment_only_are_missing_header() {
    assert_eq!(parse_featurecounts(""), Err(CountError::MissingHeader));
    assert_eq!(
        parse_featurecounts("# Program:featureCounts\n"),
        Err(CountError::MissingHeader)
    );
}

#[test]
fn featurecounts_short_row_is_bad_row_with_line() {
    let tsv = format!("# Program:featureCounts\n{FC_HDR}g1\tc\t1\t2\t+\t2\t5\ng2\tc\t1\t2\t+\t2\n");
    assert_eq!(parse_featurecounts(&tsv), Err(CountError::BadRow { line: 4 }));
}

#[test]
fn featurecounts_long_row_is_bad_row() {
    let tsv = format!("# Program:featureCounts\n{FC_HDR}g1\tc\t1\t2\t+\t2\t5\t6\n");
    assert_eq!(parse_featurecounts(&tsv), Err(CountError::BadRow { line: 3 }));
}

#[test]
fn featurecounts_blank_trailing_and_crlf_still_ok() {
    let tsv = "# Program:featureCounts\r\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ts.bam\r\ng1\tc\t1\t2\t+\t2\t5\r\n\r\ng2\tc\t1\t2\t+\t2\t6\r\n\r\n";
    assert_eq!(parse_featurecounts(tsv).unwrap(), map(&[("g1", 5), ("g2", 6)]));
}

#[test]
fn featurecounts_duplicate_gene_is_error() {
    let tsv = format!("{FC_HDR}g1\tc\t1\t2\t+\t2\t5\ng1\tc\t1\t2\t+\t2\t6\n");
    assert_eq!(
        parse_featurecounts(&tsv),
        Err(CountError::DuplicateGene("g1".into()))
    );
}

#[test]
fn assemble_fills_zero() {
    let s1 = ("s1".to_string(), map(&[("g1", 5), ("g2", 7)]));
    let s2 = ("s2".to_string(), map(&[("g2", 9)]));
    let m = assemble(&[s1, s2], &names(&["g2", "g1", "g3"]));
    assert_eq!(m.genes, names(&["g2", "g1", "g3"]));
    assert_eq!(m.samples, names(&["s1", "s2"]));
    assert_eq!(m.counts, vec![vec![7, 9], vec![5, 0], vec![0, 0]]);
}

#[test]
fn unlisted_genes_reported() {
    let s1 = ("s1".to_string(), map(&[("g1", 5), ("gX", 1)]));
    let s2 = ("s2".to_string(), map(&[("gY", 2), ("gX", 3)]));
    let genes = names(&["g1"]);
    assert_eq!(
        unlisted_genes(&[s1.clone(), s2.clone()], &genes),
        names(&["gX", "gY"])
    );
    let m = assemble(&[s1, s2], &genes);
    assert_eq!(m.counts, vec![vec![5, 0]]); // unlisted genes are dropped
}

#[test]
fn matrix_tsv_header_is_gene_then_samples() {
    let m = CountMatrix {
        genes: names(&["g1"]),
        samples: names(&["a", "b"]),
        counts: vec![vec![1, 2]],
    };
    assert_eq!(m.to_tsv().lines().next().unwrap(), "gene\ta\tb");
}

#[test]
fn matrix_tsv_exact_text() {
    let m = CountMatrix {
        genes: names(&["g1", "g2"]),
        samples: names(&["a", "b"]),
        counts: vec![vec![1, 2], vec![30, 0]],
    };
    assert_eq!(m.to_tsv(), "gene\ta\tb\ng1\t1\t2\ng2\t30\t0\n");
}

#[test]
fn gene_rows_identical_in_both_modes() {
    let genome = read_fasta(&std::fs::read(GENOME_FA).unwrap()[..]).unwrap();
    let model = parse_gff3(&std::fs::read(GENES_GFF3).unwrap()[..]).unwrap();
    let records = build_transcripts(&genome, &model).unwrap();

    let mut tx2gene: HashMap<String, String> = HashMap::new();
    for line in tx2gene_tsv(&records).lines() {
        let (t, g) = line.split_once('\t').unwrap();
        tx2gene.insert(t.to_string(), g.to_string());
    }
    let genes_tx: BTreeSet<String> = tx2gene.values().cloned().collect();
    let genes_saf: BTreeSet<String> = saf(&model)
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(genes_tx, genes_saf);
    assert!(genes_tx.contains("geneC"), "CDS-only gene");
    assert!(genes_tx.contains("geneR"), "rRNA gene");
    assert_eq!(genes_tx.len(), 5); // geneA geneB geneC geneR geneU in the fixture

    // Simulated tables: kallisto rows for every transcript, featureCounts rows for every SAF gene.
    let mut tx_ids: Vec<&String> = tx2gene.keys().collect();
    tx_ids.sort();
    let rows: Vec<(&str, &str)> = tx_ids.iter().map(|t| (t.as_str(), "1")).collect();
    let kal = gene_counts_from_abundance(&abundance(&rows), &tx2gene).unwrap();
    let mut fc =
        String::from("# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ts.bam\n");
    for g in &genes_saf {
        fc.push_str(&format!("{g}\tc\t1\t2\t+\t2\t4\n"));
    }
    let fcm = parse_featurecounts(&fc).unwrap();

    let genes: Vec<String> = genes_tx.iter().cloned().collect();
    let a = assemble(&[("s".to_string(), kal)], &genes);
    let b = assemble(&[("s".to_string(), fcm)], &genes);
    let mut ga = a.genes.clone();
    let mut gb = b.genes.clone();
    ga.sort();
    gb.sort();
    assert_eq!(ga, gb);
    // No gene is lost: every row is non-zero in both simulated modes.
    assert!(a.counts.iter().all(|r| r[0] >= 1));
    assert!(b.counts.iter().all(|r| r[0] == 4));
}

#[test]
fn featurecounts_gene_id_starting_with_hash_after_header_is_a_gene() {
    // R29-7: '#' lines are skipped only before the header.
    let tsv = format!("# Program:featureCounts\n{FC_HDR}#g1\tc\t1\t2\t+\t2\t5\ng2\tc\t1\t2\t+\t2\t3\n");
    assert_eq!(parse_featurecounts(&tsv).unwrap(), map(&[("#g1", 5), ("g2", 3)]));
}

const SF_HEAD: &str = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";

#[test]
fn salmon_gene_sum_rounds_once() {
    let m = t2g(&[("t1", "g1"), ("t2", "g1")]);
    let a = format!("{SF_HEAD}t1\t100\t90\t1\t0.4\nt2\t100\t90\t1\t0.4\n");
    assert_eq!(gene_counts_from_salmon(&a, &m).unwrap(), map(&[("g1", 1)]));
    let b = format!("{SF_HEAD}t1\t100\t90\t1\t0.25\nt2\t100\t90\t1\t0.25\n");
    assert_eq!(gene_counts_from_salmon(&b, &m).unwrap(), map(&[("g1", 1)]));
}

#[test]
fn salmon_unknown_transcript() {
    let m = t2g(&[("t1", "g1")]);
    let a = format!("{SF_HEAD}tX\t100\t90\t1\t3\n");
    assert_eq!(
        gene_counts_from_salmon(&a, &m),
        Err(CountError::UnknownTranscript("tX".into()))
    );
}

#[test]
fn salmon_missing_column() {
    let m = t2g(&[("t1", "g1")]);
    let a = "Name\tLength\tTPM\nt1\t100\t1\n";
    assert_eq!(
        gene_counts_from_salmon(a, &m),
        Err(CountError::MissingColumn("NumReads".into()))
    );
}

#[test]
fn salmon_bad_numreads() {
    let m = t2g(&[("t1", "g1")]);
    for bad in ["-1", "nan", "inf", "abc"] {
        let a = format!("{SF_HEAD}t1\t100\t90\t1\t{bad}\n");
        assert_eq!(
            gene_counts_from_salmon(&a, &m),
            Err(CountError::BadRow { line: 2 }),
            "{bad}"
        );
    }
}

#[test]
fn salmon_zero_gene_is_present_with_zero() {
    let m = t2g(&[("t1", "g1")]);
    let a = format!("{SF_HEAD}t1\t100\t90\t0\t0\n");
    assert_eq!(gene_counts_from_salmon(&a, &m).unwrap(), map(&[("g1", 0)]));
}

#[test]
fn salmon_and_kallisto_agree_on_same_numbers() {
    let m = t2g(&[("t1", "g1"), ("t2", "g1"), ("t3", "g2")]);
    let k = "target_id\tlength\teff_length\test_counts\ttpm\n\
             t1\t100\t90\t2.5\t1\nt2\t100\t90\t3.2\t1\nt3\t100\t90\t7\t1\n";
    let s = format!(
        "{SF_HEAD}t1\t100\t90\t1\t2.5\nt2\t100\t90\t1\t3.2\nt3\t100\t90\t1\t7\n"
    );
    assert_eq!(
        gene_counts_from_salmon(&s, &m).unwrap(),
        gene_counts_from_abundance(k, &m).unwrap()
    );
}
