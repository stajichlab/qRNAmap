use qrnamap_core::genes::{gene_table, gene_table_tsv, GeneClass};
use qseq_io::{Gene, GeneModel, Strand, Transcript};
use std::collections::BTreeSet;

fn tx(kind: &str) -> Transcript {
    Transcript {
        id: format!("t_{kind}"),
        kind: kind.to_string(),
        exons: vec![(1, 10)],
        cds: vec![],
        utr5: vec![],
        utr3: vec![],
        proteins: vec![],
    }
}

fn gene(label: &str, contig: &str, strand: Strand, biotype: Option<&str>, kinds: &[&str]) -> Gene {
    Gene {
        label: label.to_string(),
        contig: contig.to_string(),
        strand,
        biotype: biotype.map(|s| s.to_string()),
        symbol: None,
        gene_id: None,
        dbxref: vec![],
        old_locus_tag: None,
        description: None,
        transcripts: kinds.iter().map(|k| tx(k)).collect(),
    }
}

fn model(genes: Vec<Gene>, mito: &[&str]) -> GeneModel {
    GeneModel {
        genes,
        mito_contigs: mito.iter().map(|s| s.to_string()).collect::<BTreeSet<_>>(),
    }
}

#[test]
fn classes_for_the_yeast_biotypes() {
    let table = [
        ("protein_coding", GeneClass::ProteinCoding),
        ("rRNA", GeneClass::Rrna),
        ("tRNA", GeneClass::Trna),
        ("snoRNA", GeneClass::OtherNcrna),
        ("snRNA", GeneClass::OtherNcrna),
        ("ncRNA", GeneClass::OtherNcrna),
        ("antisense_RNA", GeneClass::OtherNcrna),
        ("misc_RNA", GeneClass::OtherNcrna),
        ("RNase_P_RNA", GeneClass::OtherNcrna),
        ("SRP_RNA", GeneClass::OtherNcrna),
        ("telomerase_RNA", GeneClass::OtherNcrna),
        ("pseudogene", GeneClass::Pseudogene),
        ("lnc_RNA", GeneClass::OtherNcrna),
        ("V_gene_segment", GeneClass::Other),
    ];
    for (bt, want) in table {
        let m = model(vec![gene("g", "c", Strand::Plus, Some(bt), &["mRNA"])], &[]);
        let rows = gene_table(&m, &[]);
        assert_eq!(rows[0].class, want, "biotype {bt}");
        assert_eq!(rows[0].biotype, bt);
    }
    assert_eq!(GeneClass::ProteinCoding.as_str(), "protein_coding");
    assert_eq!(GeneClass::Rrna.as_str(), "rRNA");
    assert_eq!(GeneClass::Trna.as_str(), "tRNA");
    assert_eq!(GeneClass::OtherNcrna.as_str(), "other_ncRNA");
    assert_eq!(GeneClass::Pseudogene.as_str(), "pseudogene");
    assert_eq!(GeneClass::Other.as_str(), "other");
}

#[test]
fn biotype_falls_back_to_transcript_kind() {
    let m = model(
        vec![
            gene("a", "c", Strand::Plus, None, &["lnc_RNA"]),
            gene("b", "c", Strand::Plus, None, &["transcript"]),
            gene("c", "c", Strand::Plus, None, &["primary_transcript"]),
            gene("d", "c", Strand::Plus, None, &["V_gene_segment"]),
            gene("e", "c", Strand::Plus, None, &["tRNA"]),
        ],
        &[],
    );
    let rows = gene_table(&m, &[]);
    let got: Vec<(&str, GeneClass)> = rows.iter().map(|r| (r.biotype.as_str(), r.class)).collect();
    assert_eq!(
        got,
        vec![
            ("lnc_RNA", GeneClass::OtherNcrna),
            ("transcript", GeneClass::Other),
            ("primary_transcript", GeneClass::Other),
            ("V_gene_segment", GeneClass::Other),
            ("tRNA", GeneClass::Trna),
        ]
    );
}

#[test]
fn mrna_wins_over_other_transcript_kinds() {
    let m = model(
        vec![gene("a", "c", Strand::Plus, None, &["lnc_RNA", "mRNA"])],
        &[],
    );
    let rows = gene_table(&m, &[]);
    assert_eq!(rows[0].biotype, "mRNA");
    assert_eq!(rows[0].class, GeneClass::ProteinCoding);
    assert_eq!(rows[0].n_transcripts, 2);
}

#[test]
fn unknown_when_no_biotype_and_no_transcripts() {
    let m = model(vec![gene("a", "c", Strand::Minus, None, &[])], &[]);
    let rows = gene_table(&m, &[]);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].biotype, "unknown");
    assert_eq!(rows[0].class, GeneClass::Other);
    assert_eq!(rows[0].n_transcripts, 0);
    assert_eq!(rows[0].strand, '-');
}

#[test]
fn mito_from_region_and_from_extra_contig() {
    let m = model(
        vec![
            gene("n", "chrI", Strand::Plus, Some("protein_coding"), &["mRNA"]),
            gene("m", "chrM", Strand::Plus, Some("protein_coding"), &["mRNA"]),
            gene(
                "x",
                "scaf9",
                Strand::Plus,
                Some("protein_coding"),
                &["mRNA"],
            ),
        ],
        &["chrM"],
    );
    let rows = gene_table(&m, &["scaf9".to_string()]);
    let flags: Vec<bool> = rows.iter().map(|r| r.mito).collect();
    assert_eq!(flags, vec![false, true, true]);
    let rows = gene_table(&m, &[]);
    let flags: Vec<bool> = rows.iter().map(|r| r.mito).collect();
    assert_eq!(flags, vec![false, true, false]);
}

#[test]
fn extra_mito_contig_without_genes_is_fine() {
    let m = model(
        vec![gene(
            "n",
            "chrI",
            Strand::Plus,
            Some("protein_coding"),
            &["mRNA"],
        )],
        &["regionOnly"],
    );
    let rows = gene_table(&m, &["nothing_here".to_string()]);
    assert_eq!(rows.len(), 1);
    assert!(!rows[0].mito);
}

#[test]
fn tsv_header_and_rows_exact() {
    let m = model(
        vec![
            gene(
                "YAL001C",
                "chrI",
                Strand::Minus,
                Some("protein_coding"),
                &["mRNA", "mRNA"],
            ),
            gene("Q0020", "chrM", Strand::Plus, Some("rRNA"), &["rRNA"]),
            gene("odd", "chrII", Strand::Plus, None, &[]),
        ],
        &["chrM"],
    );
    let tsv = gene_table_tsv(&gene_table(&m, &[]));
    let want = "gene\tcontig\tstrand\tbiotype\tclass\tmito\tn_transcripts\n\
                YAL001C\tchrI\t-\tprotein_coding\tprotein_coding\tno\t2\n\
                Q0020\tchrM\t+\trRNA\trRNA\tyes\t1\n\
                odd\tchrII\t+\tunknown\tother\tno\t0\n";
    assert_eq!(tsv, want);
}

#[test]
fn tsv_does_not_panic_on_odd_labels() {
    let m = model(vec![gene("a\tb\nc\r", "c", Strand::Plus, None, &[])], &[]);
    let _ = gene_table_tsv(&gene_table(&m, &[]));
}

#[test]
fn rows_follow_model_order() {
    let m = model(
        vec![
            gene("z", "c", Strand::Plus, None, &[]),
            gene("a", "c", Strand::Plus, None, &[]),
            gene("m", "c", Strand::Plus, None, &[]),
        ],
        &[],
    );
    let names: Vec<String> = gene_table(&m, &[]).into_iter().map(|r| r.gene).collect();
    assert_eq!(names, vec!["z", "a", "m"]);
}

#[test]
fn gene_info_tsv_exact_bytes_for_snippet() {
    let s = "##gff-version 3\n\
chr1\tsrc\tgene\t100\t1000\t.\t+\t.\tID=g1;locus_tag=LOC1;Name=SYM%3B1;Dbxref=GeneID:1001,SGD:S001;old_locus_tag=OLD%2C1\n\
chr1\tsrc\tmRNA\t100\t1000\t.\t+\t.\tID=t1;Parent=g1;product=heat%09shock protein\n\
chr1\tsrc\texon\t100\t500\t.\t+\t.\tID=e1;Parent=t1\n\
chr1\tsrc\tCDS\t200\t500\t.\t+\t0\tID=c1;Parent=t1;protein_id=PROT1\n\
chr1\tsrc\tmRNA\t100\t800\t.\t+\t.\tID=t2;Parent=g1\n\
chr1\tsrc\texon\t100\t800\t.\t+\t.\tID=e2;Parent=t2\n\
chr1\tsrc\tCDS\t200\t800\t.\t+\t0\tID=c2;Parent=t2;protein_id=PROT2\n\
chr1\tsrc\tgene\t2000\t2500\t.\t+\t.\tID=g2;Name=NAME2;Dbxref=GeneID:1002\n\
chr1\tsrc\ttRNA\t2000\t2500\t.\t+\t.\tID=t3;Parent=g2\n\
chr1\tsrc\texon\t2000\t2500\t.\t+\t.\tID=e3;Parent=t3\n\
chr1\tsrc\tgene\t3000\t3500\t.\t-\t.\tID=g3;gene_biotype=pseudogene\n\
chr1\tsrc\tpseudogenic_transcript\t3000\t3500\t.\t-\t.\tID=t4;Parent=g3\n\
chr1\tsrc\texon\t3000\t3500\t.\t-\t.\tID=e4;Parent=t4\n";
    let m = qseq_io::parse_gff3(std::io::Cursor::new(s.as_bytes())).unwrap();
    let tsv = qrnamap_core::gene_info_tsv(&m, &[]);
    let want = "gene\tsymbol\tgene_id\told_locus_tag\ttranscripts\tproteins\tbiotype\tclass\tcontig\tstrand\tdbxref\tdescription\n\
LOC1\tSYM;1\t1001\tOLD,1\tt1,t2\tPROT1,PROT2\tmRNA\tprotein_coding\tchr1\t+\tGeneID:1001,SGD:S001\theat shock protein\n\
NAME2\tNAME2\t1002\t\tt3\t\ttRNA\ttRNA\tchr1\t+\tGeneID:1002\t\n\
g3\t\t\t\tt4\t\tpseudogene\tpseudogene\tchr1\t-\t\t\n";
    assert_eq!(tsv, want);
}

#[test]
fn gene_info_tsv_keeps_sgd_only_dbxref() {
    // SGD-native GFFs carry `SGD:S000...` xrefs and no NCBI GeneID; the link-out needs the value.
    let s = "##gff-version 3\n\
chrI\tSGD\tgene\t335\t649\t.\t+\t.\tID=S000006318;Name=YAL069W;Dbxref=SGD:S000006318\n\
chrI\tSGD\tmRNA\t335\t649\t.\t+\t.\tID=t1;Parent=S000006318\n";
    let m = qseq_io::parse_gff3(std::io::Cursor::new(s.as_bytes())).unwrap();
    let tsv = qrnamap_core::gene_info_tsv(&m, &[]);
    let mut lines = tsv.lines();
    let header = lines.next().unwrap();
    let col = header.split('\t').position(|h| h == "dbxref").unwrap();
    let row = lines.next().unwrap();
    assert_eq!(row.split('\t').nth(col).unwrap(), "SGD:S000006318");
}
