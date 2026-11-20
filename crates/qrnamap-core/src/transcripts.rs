//! Transcript FASTA and `tx2gene` from a genome and a GFF3 gene model.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::io::{self, Write};

use qseq_io::{GeneModel, Strand, Transcript};

/// One transcript sequence, in transcript orientation, uppercase.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TxRecord {
    pub id: String,
    pub gene: String,
    /// GFF type of the transcript feature (mRNA, rRNA, tRNA, ...).
    pub biotype: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    /// The gene's contig is not in the genome FASTA.
    MissingContig { gene: String, contig: String },
    /// An interval has start 0, start > end, or end beyond the contig length.
    OutOfRange { transcript: String },
    /// Two transcripts have the same id.
    DuplicateTranscript(String),
    /// Two genome FASTA records have the same name.
    DuplicateContig(String),
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BuildError::MissingContig { gene, contig } => {
                write!(f, "gene '{gene}' is on contig '{contig}', which is not in the genome FASTA")
            }
            BuildError::OutOfRange { transcript } => write!(
                f,
                "transcript '{transcript}' has an interval with start 0, start > end, or end beyond the contig"
            ),
            BuildError::DuplicateTranscript(id) => write!(f, "two transcripts have the id '{id}'"),
            BuildError::DuplicateContig(name) => {
                write!(f, "the genome FASTA has two records named '{name}'")
            }
        }
    }
}

impl std::error::Error for BuildError {}

/// The intervals a transcript is made of, before merging (ruling R4):
/// exons if present; else CDS together with the explicit UTR features; else CDS only.
fn source_intervals(t: &Transcript) -> Vec<(u64, u64)> {
    if !t.exons.is_empty() {
        t.exons.clone()
    } else {
        t.cds
            .iter()
            .chain(&t.utr5)
            .chain(&t.utr3)
            .copied()
            .collect()
    }
}

/// The exon set of a transcript (ruling R4), shared by the transcript builder, SAF and splice sites.
/// Overlapping or adjacent intervals are merged; the result is sorted ascending. Empty when the
/// transcript has no exon, CDS or UTR interval (ruling R15).
pub(crate) fn exon_set(t: &Transcript) -> Vec<(u64, u64)> {
    merge(source_intervals(t))
}

/// Sort and merge overlapping or adjacent 1-based inclusive intervals.
pub(crate) fn merge(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    v.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(v.len());
    for (s, e) in v {
        match out.last_mut() {
            Some(last) if s <= last.1.saturating_add(1) => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

pub(crate) fn strand_char(s: Strand) -> char {
    match s {
        Strand::Plus => '+',
        Strand::Minus => '-',
    }
}

/// Uppercase base on the plus strand. IUPAC letters are kept; U becomes T; anything else becomes N.
fn plus_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        c @ (b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y' | b'K' | b'M' | b'S' | b'W' | b'B'
        | b'D' | b'H' | b'V') => c,
        b'U' => b'T',
        _ => b'N',
    }
}

/// Uppercase complement, including the IUPAC ambiguity letters; anything else becomes N.
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'S' => b'S',
        b'W' => b'W',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

/// Build one record per transcript with a non-empty exon set, in gene-model order. All isoforms are
/// kept; rRNA and other non-coding transcripts are kept, with `biotype` set to the transcript kind.
pub fn build_transcripts(
    genome: &[(String, Vec<u8>)],
    model: &GeneModel,
) -> Result<Vec<TxRecord>, BuildError> {
    let mut contigs: HashMap<&str, &[u8]> = HashMap::with_capacity(genome.len());
    for (name, seq) in genome {
        if contigs.insert(name.as_str(), seq.as_slice()).is_some() {
            return Err(BuildError::DuplicateContig(name.clone()));
        }
    }
    let mut seen: HashSet<&str> = HashSet::new();
    let mut out = Vec::new();
    for gene in &model.genes {
        for t in &gene.transcripts {
            let raw = source_intervals(t);
            if raw.is_empty() {
                continue; // R15: no sequence, not an error
            }
            let contig =
                *contigs
                    .get(gene.contig.as_str())
                    .ok_or_else(|| BuildError::MissingContig {
                        gene: gene.label.clone(),
                        contig: gene.contig.clone(),
                    })?;
            let len = contig.len() as u64;
            if raw.iter().any(|&(s, e)| s == 0 || s > e || e > len) {
                return Err(BuildError::OutOfRange {
                    transcript: t.id.clone(),
                });
            }
            if !seen.insert(t.id.as_str()) {
                return Err(BuildError::DuplicateTranscript(t.id.clone()));
            }
            let ivs = merge(raw);
            let total: u64 = ivs.iter().map(|&(s, e)| e - s + 1).sum();
            let mut seq = Vec::with_capacity(total as usize);
            for &(s, e) in &ivs {
                seq.extend_from_slice(&contig[(s - 1) as usize..e as usize]);
            }
            match gene.strand {
                Strand::Plus => seq.iter_mut().for_each(|b| *b = plus_base(*b)),
                Strand::Minus => {
                    seq.reverse();
                    seq.iter_mut().for_each(|b| *b = complement(*b));
                }
            }
            out.push(TxRecord {
                id: t.id.clone(),
                gene: gene.label.clone(),
                biotype: t.kind.clone(),
                seq,
            });
        }
    }
    Ok(out)
}

/// Write records as FASTA: header `>id`, sequence lines of 60 bases.
pub fn write_fasta<W: Write>(records: &[TxRecord], mut w: W) -> io::Result<()> {
    for r in records {
        w.write_all(b">")?;
        w.write_all(r.id.as_bytes())?;
        w.write_all(b"\n")?;
        for line in r.seq.chunks(60) {
            w.write_all(line)?;
            w.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// `transcript\tgene\n` rows, one per record, no header.
pub fn tx2gene_tsv(records: &[TxRecord]) -> String {
    let mut s = String::new();
    for r in records {
        s.push_str(&r.id);
        s.push('\t');
        s.push_str(&r.gene);
        s.push('\n');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(
        exons: &[(u64, u64)],
        cds: &[(u64, u64)],
        utr5: &[(u64, u64)],
        utr3: &[(u64, u64)],
    ) -> Transcript {
        Transcript {
            id: "t".into(),
            kind: "mRNA".into(),
            exons: exons.to_vec(),
            cds: cds.to_vec(),
            utr5: utr5.to_vec(),
            utr3: utr3.to_vec(),
        }
    }

    #[test]
    fn exon_set_merges_overlapping_and_adjacent() {
        assert_eq!(
            exon_set(&t(&[(9, 12), (1, 5), (3, 8), (16, 18)], &[], &[], &[])),
            vec![(1, 12), (16, 18)]
        );
        assert_eq!(
            exon_set(&t(&[(1, 10), (2, 3)], &[], &[], &[])),
            vec![(1, 10)]
        );
        assert_eq!(
            exon_set(&t(&[(1, 4), (6, 8)], &[], &[], &[])),
            vec![(1, 4), (6, 8)]
        );
    }

    #[test]
    fn exon_set_rule_order() {
        // exons win over CDS and UTRs
        assert_eq!(
            exon_set(&t(&[(1, 4)], &[(10, 20)], &[(30, 31)], &[])),
            vec![(1, 4)]
        );
        // CDS merged with UTRs
        assert_eq!(
            exon_set(&t(&[], &[(10, 20)], &[(5, 9)], &[(21, 25), (40, 45)])),
            vec![(5, 25), (40, 45)]
        );
        // CDS only
        assert_eq!(
            exon_set(&t(&[], &[(10, 20), (30, 40)], &[], &[])),
            vec![(10, 20), (30, 40)]
        );
        // nothing
        assert!(exon_set(&t(&[], &[], &[], &[])).is_empty());
    }
}
