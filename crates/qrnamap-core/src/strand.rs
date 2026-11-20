//! Strandedness estimate from featureCounts runs on a read subsample. Pure arithmetic, no I/O.
//!
//! How the CLI (`qrnamap strand-check`) gets the three counts:
//! 1. Take the first N read pairs of the sample (default 100,000) and map them with HISAT2,
//!    without `--rna-strandness`.
//! 2. Run featureCounts three times on that SAM with the SAF of the gene model:
//!    `-s 1`, `-s 2` and `-s 0` (paired samples also get `-p --countReadPairs`, so the unit
//!    is a fragment).
//! 3. `sense` is the `Assigned` count of the `-s 1` run: fragments whose first read lies on
//!    the strand of a gene.
//!    `antisense` is the `Assigned` count of the `-s 2` run: fragments whose first read lies
//!    on the strand opposite to a gene.
//!    `ambiguous` is the `Unassigned_Ambiguity` count of the `-s 0` run: fragments that
//!    overlap two or more genes when strand is ignored. These fragments give no strand
//!    information. A fragment that overlaps a plus-strand gene and a minus-strand gene can be
//!    in `sense` and in `antisense` too, because each stranded run sees only one of the genes.
//!
//! The suggestion uses only `sense` and `antisense`. `ambiguous` is reported.

use crate::config::Strandedness;

/// A strand wins when its fraction of the stranded-assigned fragments (`sense + antisense`)
/// is greater than this value. 0.8 is a chosen value. UNMEASURED: not calibrated on real
/// libraries.
pub const STRAND_THRESHOLD: f64 = 0.8;

/// Result of [`assigned_fraction`].
#[derive(Debug, Clone, PartialEq)]
pub struct StrandReport {
    pub sense: u64,
    pub antisense: u64,
    pub ambiguous: u64,
    /// `sense / (sense + antisense)`; 0 when both are 0.
    pub sense_fraction: f64,
    /// `antisense / (sense + antisense)`; 0 when both are 0.
    pub antisense_fraction: f64,
    /// `ambiguous / (sense + antisense + ambiguous)`; 0 when all are 0.
    pub ambiguous_fraction: f64,
    /// Forward when `sense_fraction > STRAND_THRESHOLD`, Reverse when
    /// `antisense_fraction > STRAND_THRESHOLD`, otherwise Unstranded.
    pub suggested: Strandedness,
    /// True when the larger of `sense_fraction` and `antisense_fraction` is below
    /// `STRAND_THRESHOLD`, or no fragment was assigned. For a truly unstranded library this
    /// is expected; for a library said to be stranded it means the setting is wrong.
    pub low_fraction_warning: bool,
}

fn frac(n: u64, d: u64) -> f64 {
    if d == 0 {
        0.0
    } else {
        n as f64 / d as f64
    }
}

/// Fractions and a suggested strandedness from `(sense, antisense, ambiguous)` fragment
/// counts. See the module documentation for how the counts are made.
pub fn assigned_fraction(strand_counts: (u64, u64, u64)) -> StrandReport {
    let (sense, antisense, ambiguous) = strand_counts;
    let stranded = sense + antisense;
    let sense_fraction = frac(sense, stranded);
    let antisense_fraction = frac(antisense, stranded);
    let ambiguous_fraction = frac(ambiguous, stranded + ambiguous);
    let suggested = if sense_fraction > STRAND_THRESHOLD {
        Strandedness::Forward
    } else if antisense_fraction > STRAND_THRESHOLD {
        Strandedness::Reverse
    } else {
        Strandedness::Unstranded
    };
    let winning = sense_fraction.max(antisense_fraction);
    StrandReport {
        sense,
        antisense,
        ambiguous,
        sense_fraction,
        antisense_fraction,
        ambiguous_fraction,
        suggested,
        low_fraction_warning: stranded == 0 || winning < STRAND_THRESHOLD,
    }
}
