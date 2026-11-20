//! Error-free paired-end read simulator with known truth.
//!
//! Each pair comes from one transcript. Fragment length is `min(FRAGMENT_LEN, transcript length)`:
//! when a transcript is shorter than 200 bases (all mini fixture transcripts are, 70 to 120
//! bases) the fragment is the whole transcript. Read length is `min(READ_LEN, fragment length)`.
//! Read 1 is the first bases of the fragment in transcript orientation, read 2 is the reverse
//! complement of its last bases (a "forward" library: read 1 on the transcript strand).
//! A reverse (dUTP-like) library swaps the two reads. An unstranded library picks one of the two
//! orientations per pair at random. The order of pairs is shuffled. The random numbers come from
//! a 64-bit LCG with a fixed seed, so the output is deterministic.

use std::io::Write;
use std::path::Path;

use qrnamap_core::{Strandedness, TxRecord};

pub const READ_LEN: usize = 75;
pub const FRAGMENT_LEN: usize = 200;

pub struct Lcg(u64);

impl Lcg {
    pub fn new(seed: u64) -> Self {
        Lcg(seed)
    }
    /// Knuth MMIX constants. The high 32 bits are returned.
    pub fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }
    /// Uniform in `0..n` (n > 0). The modulo bias is below 1e-5 for the small n used here.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u32() as u64 % n as u64) as usize
    }
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// One simulated pair: name, read 1, read 2, and the source transcript index.
pub struct Pair {
    pub name: String,
    pub r1: Vec<u8>,
    pub r2: Vec<u8>,
    pub tx: usize,
}

/// Read and fragment lengths for `simulate_with`.
#[derive(Debug, Clone, Copy)]
pub struct SimParams {
    pub read_len: usize,
    /// Fragment length is uniform in `frag_mean - frag_spread ..= frag_mean + frag_spread`,
    /// then capped at the transcript length.
    pub frag_mean: usize,
    pub frag_spread: usize,
}

/// Simulate `counts[i]` pairs from `records[i]` with read length 75 and fragment length 200
/// (the mini fixture settings).
pub fn simulate(records: &[TxRecord], counts: &[u64], lib: Strandedness, seed: u64) -> Vec<Pair> {
    simulate_with(
        records,
        counts,
        lib,
        seed,
        SimParams { read_len: READ_LEN, frag_mean: FRAGMENT_LEN, frag_spread: 0 },
    )
}

/// Simulate `counts[i]` pairs from `records[i]`. The fragment start is uniform over the
/// positions where the fragment fits in the transcript.
pub fn simulate_with(
    records: &[TxRecord],
    counts: &[u64],
    lib: Strandedness,
    seed: u64,
    p: SimParams,
) -> Vec<Pair> {
    assert_eq!(records.len(), counts.len());
    let mut rng = Lcg::new(seed);
    let mut order: Vec<usize> = Vec::new();
    for (i, &n) in counts.iter().enumerate() {
        order.extend(std::iter::repeat_n(i, n as usize));
    }
    // Fisher-Yates shuffle.
    for i in (1..order.len()).rev() {
        let j = rng.below(i + 1);
        order.swap(i, j);
    }
    let mut out = Vec::with_capacity(order.len());
    for (k, &t) in order.iter().enumerate() {
        let seq = &records[t].seq;
        // No random draw when the spread is 0, so the mini fixture output is unchanged.
        let want = if p.frag_spread == 0 {
            p.frag_mean
        } else {
            p.frag_mean - p.frag_spread + rng.below(2 * p.frag_spread + 1)
        };
        let flen = want.min(seq.len());
        let start = rng.below(seq.len() - flen + 1);
        let frag = &seq[start..start + flen];
        let rlen = p.read_len.min(flen);
        let fwd1 = frag[..rlen].to_vec();
        let fwd2 = revcomp(&frag[flen - rlen..]);
        let reverse = match lib {
            Strandedness::Forward => false,
            Strandedness::Reverse => true,
            Strandedness::Unstranded => rng.next_u32() & 1 == 1,
        };
        let (r1, r2) = if reverse { (fwd2, fwd1) } else { (fwd1, fwd2) };
        out.push(Pair { name: format!("sim{k}_{}", records[t].id), r1, r2, tx: t });
    }
    out
}

/// Write pairs `range` as two plain FASTQ files.
pub fn write_fastq(pairs: &[Pair], r1: &Path, r2: &Path) {
    let mut a = std::io::BufWriter::new(std::fs::File::create(r1).unwrap());
    let mut b = std::io::BufWriter::new(std::fs::File::create(r2).unwrap());
    for p in pairs {
        writeln!(a, "@{}/1\n{}\n+\n{}", p.name, std::str::from_utf8(&p.r1).unwrap(), "I".repeat(p.r1.len())).unwrap();
        writeln!(b, "@{}/2\n{}\n+\n{}", p.name, std::str::from_utf8(&p.r2).unwrap(), "I".repeat(p.r2.len())).unwrap();
    }
    a.flush().unwrap();
    b.flush().unwrap();
}

/// True pair count per gene.
pub fn truth_per_gene(records: &[TxRecord], counts: &[u64]) -> std::collections::BTreeMap<String, u64> {
    let mut m = std::collections::BTreeMap::new();
    for (r, &n) in records.iter().zip(counts) {
        *m.entry(r.gene.clone()).or_insert(0) += n;
    }
    m
}
