use qseq_io::{parse_gff3, GffError, Strand};
use std::io::Cursor;

fn parse(s: &str) -> Result<qseq_io::GeneModel, GffError> {
    parse_gff3(Cursor::new(s.as_bytes().to_vec()))
}

const TWO_ISOFORMS: &str = "##gff-version 3\n\
chr1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1;Name=G1\n\
chr1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\tsrc\texon\t100\t300\t.\t+\t.\tID=e1;Parent=t1\n\
chr1\tsrc\texon\t500\t900\t.\t+\t.\tID=e2;Parent=t1\n\
chr1\tsrc\tCDS\t150\t300\t.\t+\t0\tID=c1;Parent=t1\n\
chr1\tsrc\tCDS\t500\t800\t.\t+\t0\tID=c1;Parent=t1\n\
chr1\tsrc\tmRNA\t100\t700\t.\t+\t.\tID=t2;Parent=g1\n\
chr1\tsrc\texon\t500\t700\t.\t+\t.\tID=e3;Parent=t2\n\
chr1\tsrc\texon\t100\t300\t.\t+\t.\tID=e4;Parent=t2\n";

#[test]
fn keeps_all_isoforms() {
    let m = parse(TWO_ISOFORMS).unwrap();
    assert_eq!(m.genes.len(), 1);
    let g = &m.genes[0];
    assert_eq!(g.label, "G1");
    assert_eq!(g.contig, "chr1");
    assert_eq!(g.strand, Strand::Plus);
    assert_eq!(g.transcripts.len(), 2);
    assert_eq!(g.transcripts[0].id, "t1");
    assert_eq!(g.transcripts[0].kind, "mRNA");
    assert_eq!(g.transcripts[0].exons, vec![(100, 300), (500, 900)]);
    assert_eq!(g.transcripts[0].cds, vec![(150, 300), (500, 800)]);
    assert_eq!(g.transcripts[1].id, "t2");
    // exons given out of order are sorted ascending
    assert_eq!(g.transcripts[1].exons, vec![(100, 300), (500, 700)]);
}

#[test]
fn stops_at_fasta_section() {
    let s = format!("{TWO_ISOFORMS}##FASTA\n>chr1\nACGT\nchr1\tx\tgene\t1\t2\t.\t?\t.\tID=bad\n");
    let m = parse(&s).unwrap();
    assert_eq!(m.genes.len(), 1);
}

#[test]
fn crlf_lines_parse() {
    let s = TWO_ISOFORMS.replace('\n', "\r\n");
    let m = parse(&s).unwrap();
    assert_eq!(m.genes.len(), 1);
    assert_eq!(m.genes[0].label, "G1");
    assert_eq!(m.genes[0].transcripts.len(), 2);
    assert_eq!(m.genes[0].strand, Strand::Plus);
}

#[test]
fn strand_dot_is_error() {
    let s = "##gff-version 3\n\
chr1\ts\tgene\t1\t100\t.\t.\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t.\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t.\t.\tParent=t1\n";
    match parse(s) {
        Err(GffError::BadStrand { line, value }) => {
            assert_eq!(value, ".");
            assert_eq!(line, 2); // features are checked in file order; the gene line comes first
        }
        other => panic!("expected BadStrand, got {other:?}"),
    }
}

#[test]
fn strand_on_gene_line_reported_with_gene_line() {
    let s = "chr1\ts\tgene\t1\t100\t.\t?\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n";
    match parse(s) {
        Err(GffError::BadStrand { line, value }) => {
            assert_eq!((line, value.as_str()), (1, "?"));
        }
        other => panic!("expected BadStrand, got {other:?}"),
    }
}

#[test]
fn short_line_is_error() {
    let s = "##gff-version 3\nchr1\ts\tgene\t1\t100\t.\t+\n";
    match parse(s) {
        Err(GffError::BadLine { line }) => assert_eq!(line, 2),
        other => panic!("expected BadLine, got {other:?}"),
    }
}

#[test]
fn rrna_trna_kept_with_kind() {
    let s = "chr1\ts\tgene\t10\t200\t.\t-\t.\tID=r1g;locus_tag=RRN1\n\
chr1\ts\trRNA\t10\t200\t.\t-\t.\tID=r1;Parent=r1g\n\
chr1\ts\texon\t10\t200\t.\t-\t.\tParent=r1\n\
chr1\ts\tgene\t300\t380\t.\t+\t.\tID=tg;locus_tag=TRN1\n\
chr1\ts\ttRNA\t300\t380\t.\t+\t.\tID=tr;Parent=tg\n\
chr1\ts\texon\t300\t380\t.\t+\t.\tParent=tr\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 2);
    assert_eq!(m.genes[0].transcripts[0].kind, "rRNA");
    assert_eq!(m.genes[0].strand, Strand::Minus);
    assert_eq!(m.genes[1].transcripts[0].kind, "tRNA");
    assert_eq!(m.genes[1].transcripts[0].exons, vec![(300, 380)]);
}

#[test]
fn utr_features_collected() {
    let s = "chr1\ts\tgene\t1\t1000\t.\t+\t.\tID=g\n\
chr1\ts\tmRNA\t1\t1000\t.\t+\t.\tID=t;Parent=g\n\
chr1\ts\tfive_prime_UTR\t1\t99\t.\t+\t.\tParent=t\n\
chr1\ts\texon\t1\t1000\t.\t+\t.\tParent=t\n\
chr1\ts\tCDS\t100\t900\t.\t+\t0\tParent=t\n\
chr1\ts\tthree_prime_UTR\t901\t1000\t.\t+\t.\tParent=t\n";
    let t = &parse(s).unwrap().genes[0].transcripts[0];
    assert_eq!(t.utr5, vec![(1, 99)]);
    assert_eq!(t.utr3, vec![(901, 1000)]);
    assert_eq!(t.cds, vec![(100, 900)]);
    assert_eq!(t.exons, vec![(1, 1000)]);
}

#[test]
fn duplicate_gene_label_is_error() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1;Name=SAME\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n\
chr2\ts\tgene\t1\t100\t.\t+\t.\tID=g2;Name=SAME\n\
chr2\ts\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g2\n\
chr2\ts\texon\t1\t100\t.\t+\t.\tParent=t2\n";
    match parse(s) {
        Err(GffError::DuplicateGeneLabel(l)) => assert_eq!(l, "SAME"),
        other => panic!("expected DuplicateGeneLabel, got {other:?}"),
    }
}

#[test]
fn label_prefers_locus_tag_over_name_over_id() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1;Name=n1;locus_tag=LT1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n\
chr1\ts\tgene\t200\t300\t.\t+\t.\tID=g2;Name=n2\n\
chr1\ts\tmRNA\t200\t300\t.\t+\t.\tID=t2;Parent=g2\n\
chr1\ts\texon\t200\t300\t.\t+\t.\tParent=t2\n\
chr1\ts\tgene\t400\t500\t.\t+\t.\tID=g3\n\
chr1\ts\tmRNA\t400\t500\t.\t+\t.\tID=t3;Parent=g3\n\
chr1\ts\texon\t400\t500\t.\t+\t.\tParent=t3\n";
    let m = parse(s).unwrap();
    let labels: Vec<&str> = m.genes.iter().map(|g| g.label.as_str()).collect();
    assert_eq!(labels, vec!["LT1", "n2", "g3"]);
}

#[test]
fn multi_parent_cds_and_shared_exon() {
    let s = "chr1\ts\tgene\t1\t300\t.\t+\t.\tID=g\n\
chr1\ts\tmRNA\t1\t300\t.\t+\t.\tID=a;Parent=g\n\
chr1\ts\tmRNA\t1\t300\t.\t+\t.\tID=b;Parent=g\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=a,b\n\
chr1\ts\texon\t200\t300\t.\t+\t.\tParent=a\n\
chr1\ts\tCDS\t10\t90\t.\t+\t0\tParent=a,b\n";
    let g = &parse(s).unwrap().genes[0];
    assert_eq!(g.transcripts.len(), 2);
    assert_eq!(g.transcripts[0].exons, vec![(1, 100), (200, 300)]);
    assert_eq!(g.transcripts[1].exons, vec![(1, 100)]);
    assert_eq!(g.transcripts[0].cds, vec![(10, 90)]);
    assert_eq!(g.transcripts[1].cds, vec![(10, 90)]);
}

#[test]
fn noncoding_exon_and_pseudogenic_exon_are_exons() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr1\ts\ttRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\tnoncoding_exon\t1\t40\t.\t+\t.\tParent=t1\n\
chr1\ts\tnoncoding_exon\t60\t100\t.\t+\t.\tParent=t1\n\
chr1\ts\tpseudogene\t200\t300\t.\t-\t.\tID=p1\n\
chr1\ts\tpseudogenic_transcript\t200\t300\t.\t-\t.\tID=pt1;Parent=p1\n\
chr1\ts\tpseudogenic_exon\t200\t300\t.\t-\t.\tParent=pt1\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 2);
    assert_eq!(m.genes[0].transcripts[0].exons, vec![(1, 40), (60, 100)]);
    assert_eq!(m.genes[1].label, "p1");
    assert_eq!(m.genes[1].transcripts[0].kind, "pseudogenic_transcript");
    assert_eq!(m.genes[1].transcripts[0].exons, vec![(200, 300)]);
}

#[test]
fn childless_gene_is_kept_with_no_transcripts() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1;Name=lonely\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 1);
    assert_eq!(m.genes[0].label, "lonely");
    assert_eq!(m.genes[0].strand, Strand::Plus);
    assert!(m.genes[0].transcripts.is_empty());
}

#[test]
fn childless_transcript_is_kept_with_empty_intervals() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr1\ts\tncRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 1);
    let t = &m.genes[0].transcripts;
    assert_eq!(t.len(), 1);
    assert_eq!(t[0].id, "t1");
    assert_eq!(t[0].kind, "ncRNA");
    assert!(t[0].exons.is_empty() && t[0].cds.is_empty() && t[0].utr5.is_empty() && t[0].utr3.is_empty());
}

#[test]
fn exon_with_unknown_parent_is_error() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=nope\n";
    match parse(s) {
        Err(GffError::UnresolvedParent { line, parent }) => {
            assert_eq!((line, parent.as_str()), (3, "nope"));
        }
        other => panic!("expected UnresolvedParent, got {other:?}"),
    }
}

#[test]
fn exon_without_parent_is_error() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tID=e1\n";
    match parse(s) {
        Err(GffError::UnresolvedParent { line, parent }) => {
            assert_eq!((line, parent.as_str()), (2, ""));
        }
        other => panic!("expected UnresolvedParent, got {other:?}"),
    }
}

#[test]
fn child_line_before_parent_line_parses() {
    let s = "chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 1);
    assert_eq!(m.genes[0].label, "g1");
    assert_eq!(m.genes[0].transcripts[0].exons, vec![(1, 100)]);
}

#[test]
fn transcript_with_region_parent_is_own_gene() {
    let s = "chr1\ts\tregion\t1\t5000\t.\t.\t.\tID=chr1reg\n\
chr1\ts\tlnc_RNA\t10\t100\t.\t+\t.\tID=r1;Parent=chr1reg\n\
chr1\ts\texon\t10\t100\t.\t+\t.\tParent=r1\n\
chr1\ts\tlnc_RNA\t200\t300\t.\t-\t.\tID=r2;Parent=chr1reg\n\
chr1\ts\texon\t200\t300\t.\t-\t.\tParent=r2\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 2);
    assert_eq!(m.genes[0].label, "r1");
    assert_eq!(m.genes[0].strand, Strand::Plus);
    assert_eq!(m.genes[1].label, "r2");
    assert_eq!(m.genes[1].strand, Strand::Minus);
    assert_eq!(m.genes[1].transcripts.len(), 1);
}

#[test]
fn two_transcripts_one_gene_no_duplicate_label_error() {
    let m = parse(TWO_ISOFORMS).unwrap();
    assert_eq!(m.genes.len(), 1);
    assert_eq!(m.genes[0].transcripts.len(), 2);
}

#[test]
fn percent_decoding_and_orphan_transcript_becomes_gene() {
    let s = "chr1\ts\tlnc_RNA\t1\t100\t.\t-\t.\tID=rna1;Name=abc%2Cdef%20x\n\
chr1\ts\texon\t1\t100\t.\t-\t.\tParent=rna1\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 1);
    assert_eq!(m.genes[0].label, "abc,def x");
    assert_eq!(m.genes[0].strand, Strand::Minus);
    assert_eq!(m.genes[0].transcripts[0].kind, "lnc_RNA");
    assert_eq!(m.genes[0].transcripts[0].id, "rna1");
}

// R29-6: strand, contig and repeated-ID checks.

#[test]
fn transcript_strand_differs_from_gene_is_error() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t-\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t-\t.\tParent=t1\n";
    match parse(s) {
        Err(GffError::StrandMismatch { line, .. }) => assert_eq!(line, 2),
        other => panic!("expected StrandMismatch, got {other:?}"),
    }
    let e = parse(s).unwrap_err().to_string();
    assert!(e.contains("line 2") && e.contains("strand"), "{e}");
}

#[test]
fn transcript_contig_differs_from_gene_is_error() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr2\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr2\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n";
    match parse(s) {
        Err(GffError::ContigMismatch { line, .. }) => assert_eq!(line, 2),
        other => panic!("expected ContigMismatch, got {other:?}"),
    }
}

#[test]
fn child_contig_differs_from_its_transcript_is_error() {
    let s = "chr1\ts\tgene\t1\t100\t.\t+\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t50\t.\t+\t.\tParent=t1\n\
chr2\ts\texon\t60\t100\t.\t+\t.\tParent=t1\n";
    match parse(s) {
        Err(GffError::ContigMismatch { line, .. }) => assert_eq!(line, 4),
        other => panic!("expected ContigMismatch, got {other:?}"),
    }
}

#[test]
fn repeated_id_with_different_type_or_place_is_error() {
    for (second, what) in [
        ("chr1\ts\tncRNA\t200\t300\t.\t+\t.\tID=t1;Parent=g1\n", "type"),
        ("chr1\ts\tmRNA\t200\t300\t.\t-\t.\tID=t1;Parent=g1\n", "strand"),
        ("chr2\ts\tmRNA\t200\t300\t.\t+\t.\tID=t1;Parent=g1\n", "contig"),
    ] {
        let s = format!(
            "chr1\ts\tgene\t1\t300\t.\t+\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n{second}\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n"
        );
        match parse(&s) {
            Err(GffError::DuplicateId { id, line }) => assert_eq!((id.as_str(), line), ("t1", 3), "{what}"),
            other => panic!("{what}: expected DuplicateId, got {other:?}"),
        }
    }
}

#[test]
fn multi_line_feature_with_one_id_still_parses() {
    // GFF3 allows one discontinuous feature on several lines with the same ID.
    let s = "chr1\ts\tgene\t1\t300\t.\t+\t.\tID=g1\n\
chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\tmRNA\t200\t300\t.\t+\t.\tID=t1;Parent=g1\n\
chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1\n\
chr1\ts\texon\t200\t300\t.\t+\t.\tParent=t1\n";
    let m = parse(s).unwrap();
    assert_eq!(m.genes.len(), 1);
    assert_eq!(m.genes[0].transcripts.len(), 1);
    assert_eq!(m.genes[0].transcripts[0].exons, vec![(1, 100), (200, 300)]);
}

#[test]
fn gene_biotype_attribute_is_kept() {
    let m = parse("chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;gene_biotype=rRNA\n").unwrap();
    assert_eq!(m.genes[0].biotype, Some("rRNA".to_string()));
}

#[test]
fn missing_gene_biotype_is_none() {
    let m = parse("chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n").unwrap();
    assert_eq!(m.genes[0].biotype, None);
    let m = parse("chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;gene_biotype=\n").unwrap();
    assert_eq!(m.genes[0].biotype, None);
}

#[test]
fn percent_encoded_biotype_is_decoded() {
    let m = parse("chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;gene_biotype=lnc%5FRNA%3Bx\n").unwrap();
    assert_eq!(m.genes[0].biotype, Some("lnc_RNA;x".to_string()));
}

#[test]
fn mito_region_contig_is_collected() {
    let m = parse(
        "NC_001224.1\tRefSeq\tregion\t1\t85779\t.\t+\t.\tID=NC_001224.1:1..85779;Dbxref=taxon:559292;Is_circular=true;Name=MT;gbkey=Src;genome=mitochondrion;mol_type=genomic DNA;strain=S288C;substrain=FY1679\n\
         chrM2\t.\tregion\t1\t10\t.\t+\t.\tID=r2;genome=mitochondrial\n",
    )
    .unwrap();
    assert!(m.mito_contigs.contains("NC_001224.1"));
    assert!(m.mito_contigs.contains("chrM2"));
    assert_eq!(m.mito_contigs.len(), 2);
}

#[test]
fn chromosome_region_is_not_mito() {
    let m = parse("NC_001133.9\tRefSeq\tregion\t1\t230218\t.\t+\t.\tID=NC_001133.9:1..230218;Name=I;gbkey=Src;genome=chromosome;mol_type=genomic DNA\n").unwrap();
    assert!(m.mito_contigs.is_empty());
}

#[test]
fn mito_region_without_genes_is_not_an_error() {
    let m = parse("chrM\t.\tregion\t1\t10\t.\t+\t.\tID=r;genome=mitochondrion\n").unwrap();
    assert!(m.genes.is_empty());
    assert_eq!(m.mito_contigs.len(), 1);
}
