//! Larger synthetic fixture, generated deterministically (ruling R24 (3)).
//!
//! The generator owns the gene structure. It writes the genome FASTA and the GFF3 from that
//! structure, and it builds the TRUE transcript sequences itself (exon slices of its own
//! genome, joined in genomic order, reverse-complemented for minus-strand genes). It does not
//! call `build_transcripts`, so a bug in the builder shows up as a difference from the truth.
//!
//! Layout: 2 contigs of 20,000 bases. 12 genes, 6 per contig, alternating strands. Exons are
//! 150 to 450 bases, introns 80 to 300 bases, with canonical motifs written into the introns
//! (GT..AG on plus, CT..AC in plus-strand view on minus). Genes:
//! * g02 (minus): 2 isoforms, the second skips exon 2 (shared exons 1, 3, 4).
//! * g07 (plus): 3 isoforms: all 6 exons; skip exon 3; exons 1 to 4 only.
//! * g05 (plus): one rRNA, a single 1,500-base exon.
//! * g08 (minus): CDS-only mRNA (CDS lines, no exon lines).
//! * the rest: one mRNA with 2 to 6 exons.
//!
//! True pair counts per isoform are in `COUNTS_LARGE` (sum 40,000). g11 has 150 pairs, below
//! the 200-pair cut used by the tolerance checks.

use std::fmt::Write as _;

use super::sim::{revcomp, Lcg};
use qrnamap_core::TxRecord;

pub const CONTIG_LEN: usize = 20_000;

#[derive(Debug, Clone)]
pub struct Iso {
    pub id: String,
    pub kind: &'static str,
    /// 1-based inclusive genomic intervals, ascending.
    pub exons: Vec<(usize, usize)>,
    pub cds_only: bool,
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub id: String,
    pub contig: String,
    pub plus: bool,
    pub isos: Vec<Iso>,
}

pub struct Large {
    pub genome: Vec<(String, Vec<u8>)>,
    pub genes: Vec<Gene>,
    /// True transcripts in gene order, isoforms in order. Built by the generator.
    pub truth: Vec<TxRecord>,
    /// True pairs per entry of `truth`.
    pub counts: Vec<u64>,
}

/// (gene id, contig index, plus strand, exon count, kind). Kinds: "mrna", "skip2", "g07",
/// "rrna", "cds".
const PLAN: [(&str, usize, bool, usize, &str); 12] = [
    ("g01", 0, true, 3, "mrna"),
    ("g02", 0, false, 4, "skip2"),
    ("g03", 0, true, 2, "mrna"),
    ("g04", 0, false, 5, "mrna"),
    ("g05", 0, true, 1, "rrna"),
    ("g06", 0, false, 3, "mrna"),
    ("g07", 1, true, 6, "g07"),
    ("g08", 1, false, 3, "cds"),
    ("g09", 1, true, 4, "mrna"),
    ("g10", 1, false, 2, "mrna"),
    ("g11", 1, true, 3, "mrna"),
    ("g12", 1, false, 5, "mrna"),
];

/// True pairs per isoform, in the order of `Large::truth`. Sum 40,000.
pub const COUNTS_LARGE: [u64; 15] = [
    4000, // g01
    2500, 1500, // g02 a, b
    3000, // g03
    2000, // g04
    5000, // g05 rRNA
    1000, // g06
    3000, 2000, 1000, // g07 a, b, c
    2500, // g08 CDS-only
    6000, // g09
    3000, // g10
    150,  // g11
    3350, // g12
];

fn range(rng: &mut Lcg, lo: usize, hi: usize) -> usize {
    lo + rng.below(hi - lo + 1)
}

pub fn generate() -> Large {
    let mut rng = Lcg::new(0x5eed_1a7e);
    let mut genome: Vec<(String, Vec<u8>)> = ["chr1", "chr2"]
        .iter()
        .map(|n| {
            let seq = (0..CONTIG_LEN).map(|_| b"ACGT"[rng.below(4)]).collect();
            (n.to_string(), seq)
        })
        .collect();
    let mut cursor = [500usize, 500usize];
    let mut genes = Vec::new();
    for (id, c, plus, n_exons, kind) in PLAN {
        let mut exons = Vec::new();
        let mut pos = cursor[c];
        for i in 0..n_exons {
            let len = if kind == "rrna" {
                1500
            } else {
                range(&mut rng, 150, 450)
            };
            exons.push((pos, pos + len - 1));
            pos += len;
            if i + 1 < n_exons {
                pos += range(&mut rng, 80, 300);
            }
        }
        assert!(pos + 500 < CONTIG_LEN, "gene {id} does not fit");
        cursor[c] = pos + range(&mut rng, 600, 1200);
        // Canonical intron motifs between consecutive exons.
        let seq = &mut genome[c].1;
        for w in exons.windows(2) {
            let (s, e) = (w[0].1 + 1, w[1].0 - 1); // intron, 1-based inclusive
            let (donor, acceptor) = if plus { (b"GT", b"AG") } else { (b"CT", b"AC") };
            seq[s - 1..s + 1].copy_from_slice(donor);
            seq[e - 2..e].copy_from_slice(acceptor);
        }
        let iso = |suffix: &str, ex: Vec<(usize, usize)>, kind: &'static str, cds_only: bool| Iso {
            id: format!("{id}{suffix}"),
            kind,
            exons: ex,
            cds_only,
        };
        let isos = match kind {
            "skip2" => vec![
                iso(".a", exons.clone(), "mRNA", false),
                iso(
                    ".b",
                    exons
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != 1)
                        .map(|(_, e)| *e)
                        .collect(),
                    "mRNA",
                    false,
                ),
            ],
            "g07" => vec![
                iso(".a", exons.clone(), "mRNA", false),
                iso(
                    ".b",
                    exons
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != 2)
                        .map(|(_, e)| *e)
                        .collect(),
                    "mRNA",
                    false,
                ),
                iso(".c", exons[..4].to_vec(), "mRNA", false),
            ],
            "rrna" => vec![iso(".r", exons.clone(), "rRNA", false)],
            "cds" => vec![iso(".a", exons.clone(), "mRNA", true)],
            _ => vec![iso(".a", exons.clone(), "mRNA", false)],
        };
        genes.push(Gene {
            id: id.to_string(),
            contig: genome[c].0.clone(),
            plus,
            isos,
        });
    }

    // True transcripts from the generator's own structure.
    let mut truth = Vec::new();
    for g in &genes {
        let contig = &genome.iter().find(|(n, _)| *n == g.contig).unwrap().1;
        for iso in &g.isos {
            let mut s: Vec<u8> = Vec::new();
            for &(a, b) in &iso.exons {
                s.extend_from_slice(&contig[a - 1..b]);
            }
            if !g.plus {
                s = revcomp(&s);
            }
            truth.push(TxRecord {
                id: iso.id.clone(),
                gene: g.id.clone(),
                biotype: iso.kind.into(),
                seq: s,
            });
        }
    }
    assert_eq!(truth.len(), COUNTS_LARGE.len());
    Large {
        genome,
        genes,
        truth,
        counts: COUNTS_LARGE.to_vec(),
    }
}

impl Large {
    pub fn fasta(&self) -> String {
        let mut s = String::new();
        for (n, seq) in &self.genome {
            writeln!(s, ">{n}").unwrap();
            for line in seq.chunks(60) {
                s.push_str(std::str::from_utf8(line).unwrap());
                s.push('\n');
            }
        }
        s
    }

    /// GFF3 written from the generator's gene structure.
    pub fn gff3(&self) -> String {
        let mut s = String::from("##gff-version 3\n");
        for (n, seq) in &self.genome {
            writeln!(s, "##sequence-region {n} 1 {}", seq.len()).unwrap();
        }
        for g in &self.genes {
            let st = if g.plus { '+' } else { '-' };
            let lo = g.isos.iter().map(|i| i.exons[0].0).min().unwrap();
            let hi = g
                .isos
                .iter()
                .map(|i| i.exons.last().unwrap().1)
                .max()
                .unwrap();
            writeln!(
                s,
                "{}\tgen\tgene\t{lo}\t{hi}\t.\t{st}\t.\tID={};Name={}",
                g.contig, g.id, g.id
            )
            .unwrap();
            for iso in &g.isos {
                let (a, b) = (iso.exons[0].0, iso.exons.last().unwrap().1);
                writeln!(
                    s,
                    "{}\tgen\t{}\t{a}\t{b}\t.\t{st}\t.\tID={};Parent={}",
                    g.contig, iso.kind, iso.id, g.id
                )
                .unwrap();
                for &(x, y) in &iso.exons {
                    let (ty, phase) = if iso.cds_only {
                        ("CDS", "0")
                    } else {
                        ("exon", ".")
                    };
                    writeln!(
                        s,
                        "{}\tgen\t{ty}\t{x}\t{y}\t.\t{st}\t{phase}\tParent={}",
                        g.contig, iso.id
                    )
                    .unwrap();
                }
            }
        }
        s
    }

    /// True pairs per gene.
    pub fn truth_per_gene(&self) -> std::collections::BTreeMap<String, u64> {
        super::sim::truth_per_gene(&self.truth, &self.counts)
    }
}
