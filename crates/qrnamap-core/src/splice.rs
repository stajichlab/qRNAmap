//! Known splice sites for HISAT2 (`--known-splicesite-infile`), the Rust equivalent of
//! `hisat2_extract_splice_sites.py`.

use std::collections::BTreeSet;

use qseq_io::GeneModel;

use crate::transcripts::{exon_set, strand_char};

/// One line per distinct intron: `chrom\tleft\tright\tstrand\n`. `left` is the 0-based position of the
/// last base of the upstream exon, `right` the 0-based position of the first base of the downstream exon.
/// Introns come from the exon set of each transcript (ruling R4). As in `hisat2_extract_splice_sites.py`,
/// two exons whose start minus previous end is 5 or less are treated as one exon (no junction).
/// Sorted by chrom (byte order), left, right, strand.
pub fn splice_sites(model: &GeneModel) -> String {
    let mut sites: BTreeSet<(&str, u64, u64, char)> = BTreeSet::new();
    for gene in &model.genes {
        let strand = strand_char(gene.strand);
        for t in &gene.transcripts {
            let ex = exon_set(t);
            for w in ex.windows(2) {
                let (prev_end, next_start) = (w[0].1, w[1].0);
                if next_start - prev_end <= 5 {
                    continue;
                }
                // next_start > prev_end + 5 here, so next_start - 1 cannot underflow; prev_end is 0 only
                // for a start-0 interval, which the GFF3 parser rejects.
                sites.insert((
                    gene.contig.as_str(),
                    prev_end.saturating_sub(1),
                    next_start - 1,
                    strand,
                ));
            }
        }
    }
    let mut s = String::new();
    for (chrom, left, right, strand) in sites {
        s.push_str(&format!("{chrom}\t{left}\t{right}\t{strand}\n"));
    }
    s
}
