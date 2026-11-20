//! Transcriptome builder, SAF and splice-site tests.
//! Expected sequences are literals from `tests/fixtures/mini/make_fixture.py` (an independent
//! Python implementation) or are derived by hand in the comments; none come from the code under test.

use qrnamap_core::{
    build_transcripts, saf, splice_sites, tx2gene_tsv, write_fasta, BuildError, TxRecord,
};
use qseq_io::{parse_gff3, read_fasta, Gene, GeneModel, Strand, Transcript};

const GENOME_FA: &[u8] = include_bytes!("../../../tests/fixtures/mini/genome.fa");
const GENES_GFF3: &[u8] = include_bytes!("../../../tests/fixtures/mini/genes.gff3");

const TA1: &str = "AAGTGCCTAGCATTACGGTCAAATCTGGTTATACTCGTTGCGTGTACTTCACGGATTATGCATGTGTCGCTTACGATCGATCGTTCATCGTGTGCTCGTATGCCAACGAT";
const TB1: &str = "TCCTCTACATGGCAATTGCCGTTTCAGCTTTCGCTGATCTTGCCGCATCTTCTTGAGCAATATAGCAGTGAAAAAACCGATCCCCAAAGAGTACAAGTGAACCATTATCAGAGAGATGCT";
const TB2: &str = "GAGACTGGAGCGTAGCCCACAGCGCAAGTATCCTCTACATGGCAATTGCCTGCCGCATCTTCTTGAGCAATATAGCAGTGAAAAAACCGATCCCCAAAGAGTACAAGTGA";
const TR1: &str =
    "TCGAAGCGCGTTTTAAATTCTGTTGGCGCGCCGACACGGTACCGAAACGTCCAAATTGAGCAATAATCCTCTGATCCAAG";
const TC1: &str = "ATAAAGTAACCACGGAGTAACCGTTCCGCCCCCGCGCGTGGTTCAACTCCGATGCGAAGGTTGGCGGCGT";
const TU1: &str =
    "GTAATAAACGCTTCAGCACTGTTGACAACCACCTCGATTGGCACATGGCGCAAATTAGGCGCGATGTCCCCCGGTGGGGC";

fn genome() -> Vec<(String, Vec<u8>)> {
    read_fasta(GENOME_FA).unwrap()
}

fn model() -> GeneModel {
    parse_gff3(GENES_GFF3).unwrap()
}

fn fixture_records() -> Vec<TxRecord> {
    build_transcripts(&genome(), &model()).unwrap()
}

fn rec<'a>(recs: &'a [TxRecord], id: &str) -> &'a TxRecord {
    recs.iter()
        .find(|r| r.id == id)
        .unwrap_or_else(|| panic!("no record {id}"))
}

fn seq(r: &TxRecord) -> &str {
    std::str::from_utf8(&r.seq).unwrap()
}

// Small in-memory genome. Positions: 1A 2C 3G 4T 5T 6G 7C 8A 9A 10G 11G 12C 13C 14T 15T 16A 17A 18C 19G 20T
const CTG: &str = "ACGTTGCAAGGCCTTAACGT";

fn small_genome() -> Vec<(String, Vec<u8>)> {
    vec![("ctg".to_string(), CTG.as_bytes().to_vec())]
}

fn tx(id: &str, exons: &[(u64, u64)]) -> Transcript {
    Transcript {
        id: id.to_string(),
        kind: "mRNA".to_string(),
        exons: exons.to_vec(),
        cds: vec![],
        utr5: vec![],
        utr3: vec![],
    }
}

fn gene(label: &str, contig: &str, strand: Strand, transcripts: Vec<Transcript>) -> Gene {
    Gene {
        label: label.to_string(),
        contig: contig.to_string(),
        strand,
        transcripts,
    }
}

// ---- brief tests on the fixture ----

#[test]
fn plus_strand_two_exons_concatenate() {
    let recs = fixture_records();
    let r = rec(&recs, "tA1");
    assert_eq!(r.gene, "geneA");
    // exons 101-160 + 301-350 = 110 bases; exon 1 is soft-masked at 141-160 in genome.fa
    assert_eq!(r.seq.len(), 110);
    assert_eq!(seq(r), TA1);
}

#[test]
fn minus_strand_reverse_complemented() {
    let recs = fixture_records();
    let r = rec(&recs, "tB1");
    assert_eq!(r.gene, "geneB");
    assert_eq!(seq(r), TB1);
}

#[test]
fn both_isoforms_emitted() {
    let recs = fixture_records();
    let b: Vec<&TxRecord> = recs.iter().filter(|r| r.gene == "geneB").collect();
    assert_eq!(b.len(), 2);
    assert_eq!(seq(rec(&recs, "tB1")), TB1);
    assert_eq!(seq(rec(&recs, "tB2")), TB2);
}

#[test]
fn cds_only_gene_uses_cds() {
    let recs = fixture_records();
    let r = rec(&recs, "tC1");
    assert_eq!(r.gene, "geneC");
    assert_eq!(seq(r), TC1);
}

#[test]
fn rrna_in_index_with_biotype() {
    let recs = fixture_records();
    let r = rec(&recs, "tR1");
    assert_eq!(r.biotype, "rRNA");
    assert_eq!(r.gene, "geneR");
    assert_eq!(seq(r), TR1);
    assert_eq!(rec(&recs, "tA1").biotype, "mRNA");
}

#[test]
fn missing_contig_is_named_error() {
    let g: Vec<(String, Vec<u8>)> = genome().into_iter().filter(|(n, _)| n != "chrB").collect();
    let err = build_transcripts(&g, &model()).unwrap_err();
    assert_eq!(
        err,
        BuildError::MissingContig {
            gene: "geneR".to_string(),
            contig: "chrB".to_string()
        }
    );
    assert!(err.to_string().contains("chrB"));
    assert!(err.to_string().contains("geneR"));
}

#[test]
fn saf_is_exon_union_per_gene() {
    // geneB: tB1 801-840 and tB2 821-870 overlap and merge into 801-870; 621-640 lies inside 601-640.
    // geneC is CDS-only; geneU is UTR+CDS (adjacent, merged); geneR is rRNA and is kept.
    let expected = "GeneID\tChr\tStart\tEnd\tStrand\n\
geneA\tchrA\t101\t160\t+\n\
geneA\tchrA\t301\t350\t+\n\
geneB\tchrA\t601\t640\t-\n\
geneB\tchrA\t701\t740\t-\n\
geneB\tchrA\t801\t870\t-\n\
geneR\tchrB\t101\t180\t+\n\
geneC\tchrB\t401\t430\t-\n\
geneC\tchrB\t501\t540\t-\n\
geneU\tchrB\t701\t780\t+\n";
    assert_eq!(saf(&model()), expected);
}

#[test]
fn tx2gene_has_one_row_per_transcript() {
    let expected = "tA1\tgeneA\ntB1\tgeneB\ntB2\tgeneB\ntR1\tgeneR\ntC1\tgeneC\ntU1\tgeneU\n";
    assert_eq!(tx2gene_tsv(&fixture_records()), expected);
}

#[test]
fn saf_and_tx2gene_name_the_same_genes() {
    let recs = fixture_records();
    let mut from_tx: Vec<String> = recs.iter().map(|r| r.gene.clone()).collect();
    from_tx.sort();
    from_tx.dedup();
    let s = saf(&model());
    let mut from_saf: Vec<String> = s
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    from_saf.sort();
    from_saf.dedup();
    assert_eq!(from_tx, from_saf);
}

#[test]
fn utr_merged_with_cds_when_no_exons_fixture() {
    let recs = fixture_records();
    assert_eq!(seq(rec(&recs, "tU1")), TU1);
}

#[test]
fn fixture_splice_sites() {
    // Hand-derived: left = last exon base - 1, right = first exon base - 1 (0-based).
    // tA1 160|301; tB1 640|701, 740|801; tB2 640|701 (shared, once), 740|821; tC1 (CDS) 430|501.
    let expected = "chrA\t159\t300\t+\n\
chrA\t639\t700\t-\n\
chrA\t739\t800\t-\n\
chrA\t739\t820\t-\n\
chrB\t429\t500\t-\n";
    assert_eq!(splice_sites(&model()), expected);
}

// ---- extra rules on small in-memory models ----

#[test]
fn utr_merged_with_cds_when_no_exons() {
    let mut t = tx("t1", &[]);
    t.utr5 = vec![(1, 4)];
    t.cds = vec![(8, 12)];
    t.utr3 = vec![(13, 15)];
    let m = GeneModel {
        genes: vec![gene("g1", "ctg", Strand::Plus, vec![t])],
    };
    let recs = build_transcripts(&small_genome(), &m).unwrap();
    // 1-4 "ACGT" + 8-15 "AAGGCCTT"
    assert_eq!(seq(&recs[0]), "ACGTAAGGCCTT");
    assert_eq!(
        saf(&m),
        "GeneID\tChr\tStart\tEnd\tStrand\ng1\tctg\t1\t4\t+\ng1\tctg\t8\t15\t+\n"
    );
}

#[test]
fn exons_take_precedence_over_cds() {
    let mut t = tx("t1", &[(1, 4)]);
    t.cds = vec![(2, 3), (10, 12)];
    let m = GeneModel {
        genes: vec![gene("g1", "ctg", Strand::Plus, vec![t])],
    };
    let recs = build_transcripts(&small_genome(), &m).unwrap();
    assert_eq!(seq(&recs[0]), "ACGT");
}

#[test]
fn exon_set_merges_overlapping_and_adjacent_intervals() {
    // (1,5)+(3,8) overlap, (9,12) is adjacent: one interval 1-12; then 16-18.
    let t = tx("t1", &[(9, 12), (1, 5), (3, 8), (16, 18)]);
    let m = GeneModel {
        genes: vec![gene("g1", "ctg", Strand::Plus, vec![t])],
    };
    let recs = build_transcripts(&small_genome(), &m).unwrap();
    // 1-12 "ACGTTGCAAGGC" + 16-18 "AAC"
    assert_eq!(seq(&recs[0]), "ACGTTGCAAGGCAAC");
    assert_eq!(
        saf(&m),
        "GeneID\tChr\tStart\tEnd\tStrand\ng1\tctg\t1\t12\t+\ng1\tctg\t16\t18\t+\n"
    );
}

#[test]
fn duplicate_transcript_is_named_error() {
    let m = GeneModel {
        genes: vec![
            gene("g1", "ctg", Strand::Plus, vec![tx("t1", &[(1, 4)])]),
            gene("g2", "ctg", Strand::Minus, vec![tx("t1", &[(5, 8)])]),
        ],
    };
    let err = build_transcripts(&small_genome(), &m).unwrap_err();
    assert_eq!(err, BuildError::DuplicateTranscript("t1".to_string()));
    assert!(err.to_string().contains("t1"));
}

#[test]
fn out_of_range_is_named_error() {
    for iv in [(15, 21), (0, 3), (5, 3)] {
        let m = GeneModel {
            genes: vec![gene(
                "g1",
                "ctg",
                Strand::Plus,
                vec![tx("tbad", &[(1, 4), iv])],
            )],
        };
        let err = build_transcripts(&small_genome(), &m).unwrap_err();
        assert_eq!(
            err,
            BuildError::OutOfRange {
                transcript: "tbad".to_string()
            },
            "interval {iv:?}"
        );
        assert!(err.to_string().contains("tbad"));
    }
    // The last base of the contig is in range.
    let m = GeneModel {
        genes: vec![gene("g1", "ctg", Strand::Plus, vec![tx("t1", &[(18, 20)])])],
    };
    assert_eq!(
        seq(&build_transcripts(&small_genome(), &m).unwrap()[0]),
        "CGT"
    );
}

#[test]
fn empty_interval_transcript_skipped() {
    let m = GeneModel {
        genes: vec![
            gene("g0", "ctg", Strand::Plus, vec![tx("t0", &[])]),
            gene(
                "g1",
                "ctg",
                Strand::Plus,
                vec![tx("t1", &[]), tx("t2", &[(1, 4)])],
            ),
            gene("g2", "ctg", Strand::Plus, vec![]),
            // a gene with no intervals on a contig that is not in the genome needs no sequence
            gene("g3", "absent", Strand::Plus, vec![tx("t3", &[])]),
        ],
    };
    let recs = build_transcripts(&small_genome(), &m).unwrap();
    assert_eq!(recs.len(), 1);
    assert_eq!(recs[0].id, "t2");
    assert_eq!(tx2gene_tsv(&recs), "t2\tg1\n");
    assert_eq!(
        saf(&m),
        "GeneID\tChr\tStart\tEnd\tStrand\ng1\tctg\t1\t4\t+\n"
    );
    assert_eq!(splice_sites(&m), "");
}

#[test]
fn lowercase_genome_gives_uppercase_transcript() {
    let g = vec![("ctg".to_string(), b"acgtn".to_vec())];
    let m = GeneModel {
        genes: vec![
            gene("gp", "ctg", Strand::Plus, vec![tx("tp", &[(1, 5)])]),
            gene("gm", "ctg", Strand::Minus, vec![tx("tm", &[(1, 5)])]),
        ],
    };
    let recs = build_transcripts(&g, &m).unwrap();
    assert_eq!(seq(&recs[0]), "ACGTN");
    // reverse complement of ACGTN
    assert_eq!(seq(&recs[1]), "NACGT");
}

#[test]
fn iupac_codes_complemented_and_unknown_becomes_n() {
    // complements: A>T C>G R>Y Y>R K>M M>K S>S W>W B>V D>H H>D V>B N>N X>N, then reversed
    let g = vec![("ctg".to_string(), b"ACRYKMSWBDHVNx".to_vec())];
    let m = GeneModel {
        genes: vec![gene("gm", "ctg", Strand::Minus, vec![tx("tm", &[(1, 14)])])],
    };
    let recs = build_transcripts(&g, &m).unwrap();
    assert_eq!(seq(&recs[0]), "NNBDHVWSKMRYGT");
}

#[test]
fn splice_sites_deduplicated_across_isoforms() {
    let g = vec![("ctg".to_string(), [CTG, CTG].concat().into_bytes())];
    let m = GeneModel {
        genes: vec![gene(
            "g1",
            "ctg",
            Strand::Plus,
            vec![
                tx("t1", &[(1, 4), (11, 14), (21, 24)]),
                tx("t2", &[(2, 4), (11, 14)]),
            ],
        )],
    };
    assert_eq!(splice_sites(&m), "ctg\t3\t10\t+\nctg\t13\t20\t+\n");
    assert_eq!(build_transcripts(&g, &m).unwrap().len(), 2);
}

#[test]
fn single_exon_transcript_gives_no_splice_site() {
    let m = GeneModel {
        genes: vec![gene("g1", "ctg", Strand::Minus, vec![tx("t1", &[(1, 20)])])],
    };
    assert_eq!(splice_sites(&m), "");
}

#[test]
fn intron_of_five_bases_or_less_gives_no_splice_site() {
    // hisat2_extract_splice_sites.py joins exons whose start - previous end <= 5.
    let m = GeneModel {
        genes: vec![gene(
            "g1",
            "ctg",
            Strand::Plus,
            vec![tx("t1", &[(1, 4), (9, 12), (20, 22)])],
        )],
    };
    // 9 - 4 = 5: joined; 20 - 12 = 8: junction left 11, right 19
    assert_eq!(splice_sites(&m), "ctg\t11\t19\t+\n");
}

#[test]
fn splice_sites_sorted_by_chrom_then_position() {
    let m = GeneModel {
        genes: vec![
            gene(
                "g1",
                "chrZ",
                Strand::Plus,
                vec![tx("t1", &[(1, 4), (21, 24)])],
            ),
            gene(
                "g2",
                "chrA",
                Strand::Minus,
                vec![tx("t2", &[(31, 34), (51, 54)])],
            ),
            gene(
                "g3",
                "chrA",
                Strand::Plus,
                vec![tx("t3", &[(1, 4), (21, 24)])],
            ),
        ],
    };
    assert_eq!(
        splice_sites(&m),
        "chrA\t3\t20\t+\nchrA\t33\t50\t-\nchrZ\t3\t20\t+\n"
    );
}

#[test]
fn write_fasta_wraps_at_60() {
    let long: Vec<u8> = vec![b'A'; 61];
    let recs = vec![
        TxRecord {
            id: "t1".into(),
            gene: "g1".into(),
            biotype: "mRNA".into(),
            seq: b"ACGT".to_vec(),
        },
        TxRecord {
            id: "t2".into(),
            gene: "g1".into(),
            biotype: "mRNA".into(),
            seq: long,
        },
    ];
    let mut out = Vec::new();
    write_fasta(&recs, &mut out).unwrap();
    let expected = format!(">t1\nACGT\n>t2\n{}\nA\n", "A".repeat(60));
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn duplicate_contig_name_is_named_error() {
    // R29-5: a genome FASTA with two records of one name is an error, not "first wins".
    let mut genome = small_genome();
    genome.push(("ctg".to_string(), b"TTTTTTTT".to_vec()));
    let m = GeneModel { genes: vec![gene("g1", "ctg", Strand::Plus, vec![tx("t1", &[(1, 4)])])] };
    let err = build_transcripts(&genome, &m).unwrap_err();
    assert_eq!(err, BuildError::DuplicateContig("ctg".to_string()));
    assert!(err.to_string().contains("'ctg'"), "{err}");
}
