//! Gene table: class, biotype and mitochondrial flag per gene. Pure; no I/O.

use qseq_io::{Gene, GeneModel, Strand};
use std::fmt::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeneClass {
    ProteinCoding,
    Rrna,
    Trna,
    OtherNcrna,
    Pseudogene,
    Other,
}

impl GeneClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            GeneClass::ProteinCoding => "protein_coding",
            GeneClass::Rrna => "rRNA",
            GeneClass::Trna => "tRNA",
            GeneClass::OtherNcrna => "other_ncRNA",
            GeneClass::Pseudogene => "pseudogene",
            GeneClass::Other => "other",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneRow {
    pub gene: String,
    pub contig: String,
    pub strand: char,
    pub biotype: String,
    pub class: GeneClass,
    pub mito: bool,
    pub n_transcripts: usize,
}

fn biotype_of(gene: &Gene) -> String {
    if let Some(b) = &gene.biotype {
        return b.clone();
    }
    if gene.transcripts.iter().any(|t| t.kind == "mRNA") {
        return "mRNA".to_string();
    }
    match gene.transcripts.first() {
        Some(t) => t.kind.clone(),
        None => "unknown".to_string(),
    }
}

fn class_of(biotype: &str) -> GeneClass {
    let b = biotype.to_lowercase();
    if b == "protein_coding" || b == "mrna" {
        GeneClass::ProteinCoding
    } else if b == "rrna" {
        GeneClass::Rrna
    } else if b == "trna" {
        GeneClass::Trna
    } else if b.contains("pseudogene") || b.contains("pseudogenic") {
        GeneClass::Pseudogene
    } else if b.ends_with("rna") {
        GeneClass::OtherNcrna
    } else {
        GeneClass::Other
    }
}

/// One row per gene, in model order.
pub fn gene_table(model: &GeneModel, extra_mito_contigs: &[String]) -> Vec<GeneRow> {
    model
        .genes
        .iter()
        .map(|g| {
            let biotype = biotype_of(g);
            let class = class_of(&biotype);
            let mito =
                model.mito_contigs.contains(&g.contig) || extra_mito_contigs.contains(&g.contig);
            GeneRow {
                gene: g.label.clone(),
                contig: g.contig.clone(),
                strand: match g.strand {
                    Strand::Plus => '+',
                    Strand::Minus => '-',
                },
                biotype,
                class,
                mito,
                n_transcripts: g.transcripts.len(),
            }
        })
        .collect()
}

/// Tab-separated table with a header line; every line ends with a newline.
///
/// Labels and contig names are written as they are: this function does not escape or check
/// them. A gene label with a tab, newline or carriage return gives a broken table. Callers must
/// refuse such labels first (the CLI does this in `write_genes_table`, `BadGeneLabel`); a caller
/// without that check (for example a browser build) must add it.
pub fn gene_table_tsv(rows: &[GeneRow]) -> String {
    let mut s = String::from("gene\tcontig\tstrand\tbiotype\tclass\tmito\tn_transcripts\n");
    for r in rows {
        let _ = writeln!(
            s,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.gene,
            r.contig,
            r.strand,
            r.biotype,
            r.class.as_str(),
            if r.mito { "yes" } else { "no" },
            r.n_transcripts
        );
    }
    s
}

/// Header and rows of the gene alias table `gene_info.tsv`.
///
/// Header: `gene\tsymbol\tgene_id\told_locus_tag\ttranscripts\tproteins\tbiotype\tclass\tcontig\tstrand\n`.
/// Unknown attributes are ignored; a parse error in one alias never fails the transcriptome build
/// (the alias cell is empty).
pub fn gene_info_tsv(model: &GeneModel, _extra_mito_contigs: &[String]) -> String {
    let mut s = String::from(
        "gene\tsymbol\tgene_id\told_locus_tag\ttranscripts\tproteins\tbiotype\tclass\tcontig\tstrand\n",
    );
    for g in &model.genes {
        let biotype = biotype_of(g);
        let class = class_of(&biotype);
        let strand = match g.strand {
            Strand::Plus => '+',
            Strand::Minus => '-',
        };
        let symbol = g.symbol.as_deref().unwrap_or("");
        let gene_id = g.gene_id.as_deref().unwrap_or("");
        let old_locus_tag = g.old_locus_tag.as_deref().unwrap_or("");

        let mut tx_ids: Vec<&str> = Vec::new();
        let mut proteins: Vec<&str> = Vec::new();
        for t in &g.transcripts {
            if !tx_ids.contains(&t.id.as_str()) {
                tx_ids.push(&t.id);
            }
            for p in &t.proteins {
                if !proteins.contains(&p.as_str()) {
                    proteins.push(p);
                }
            }
        }
        let transcripts_str = tx_ids.join(",");
        let proteins_str = proteins.join(",");

        let _ = writeln!(
            s,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.label,
            symbol,
            gene_id,
            old_locus_tag,
            transcripts_str,
            proteins_str,
            biotype,
            class.as_str(),
            g.contig,
            strand
        );
    }
    s
}
