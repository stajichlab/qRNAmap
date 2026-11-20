//! featureCounts SAF annotation from the gene model.

use qseq_io::GeneModel;

use crate::transcripts::{exon_set, merge, strand_char};

/// SAF text: header `GeneID\tChr\tStart\tEnd\tStrand`, then one row per merged interval of the union of
/// the exon sets (ruling R4) of all transcripts of a gene. 1-based inclusive. `GeneID` is the gene label,
/// the same name `tx2gene` uses. A gene with no interval in any transcript has no row (ruling R15).
pub fn saf(model: &GeneModel) -> String {
    let mut s = String::from("GeneID\tChr\tStart\tEnd\tStrand\n");
    for gene in &model.genes {
        let all: Vec<(u64, u64)> = gene.transcripts.iter().flat_map(exon_set).collect();
        let strand = strand_char(gene.strand);
        for (start, end) in merge(all) {
            s.push_str(&format!(
                "{}\t{}\t{start}\t{end}\t{strand}\n",
                gene.label, gene.contig
            ));
        }
    }
    s
}
