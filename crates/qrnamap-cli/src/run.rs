//! `qrnamap run` support: the reads-per-sample check.

use std::fmt;
use std::io::{self, BufRead};
use std::path::Path;

use qrnamap_core::SampleSheet;
use qseq_io::FastqReader;

use crate::io::{open_input, InputError};

/// Turn a read error into an `InputError`, keeping a wrapped `Truncated` or `Corrupt` error.
fn input_error(path: &Path, e: io::Error) -> InputError {
    if e.get_ref().is_some_and(|inner| inner.is::<InputError>()) {
        return *e.into_inner().and_then(|b| b.downcast::<InputError>().ok()).expect("checked with is::<InputError>()");
    }
    InputError::Io { path: path.to_path_buf(), source: e }
}

/// Count the FASTQ records in `reader`. With `limit`, stop as soon as the count is above the
/// limit and return `limit + 1`, so a large file is not read to its end. A truncated or bad
/// record is an error. `path` is used in error messages.
pub fn count_records(reader: Box<dyn BufRead>, path: &Path, limit: Option<u64>) -> Result<u64, InputError> {
    let mut n = 0u64;
    for rec in FastqReader::new(reader) {
        rec.map_err(|e| input_error(path, e))?;
        n += 1;
        if limit.is_some_and(|l| n > l) {
            break;
        }
    }
    Ok(n)
}

/// Count the records (read pairs, when `path` is a mate 1 file) of a plain, gzip or zstd FASTQ.
/// See `count_records` for `limit`.
pub fn count_pairs(path: &Path, limit: Option<u64>) -> Result<u64, InputError> {
    count_records(open_input(path)?, path, limit)
}

/// A sample with more read pairs than `reads_per_sample`. In the CLI this is a warning, never
/// an error (spec 11a).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Warning {
    pub sample: String,
    pub limit: u64,
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "sample '{}' has more than {} read pairs (reads_per_sample = {}); the run continues. Raise --reads-per-sample to silence this, or check the SLURM resources",
            self.sample, self.limit, self.limit
        )
    }
}

/// A warning for each sample whose `fastq1` files (all lanes, summed) hold more than `limit`
/// records. Reading stops as soon as a sample is over the limit.
pub fn check_reads_per_sample(sheet: &SampleSheet, limit: u64) -> Result<Vec<Warning>, InputError> {
    let mut out = Vec::new();
    for s in &sheet.samples {
        let mut total = 0u64;
        for lane in &s.fastq1 {
            total += count_pairs(Path::new(lane), Some(limit - total))?;
            if total > limit {
                out.push(Warning { sample: s.name.clone(), limit });
                break;
            }
        }
    }
    Ok(out)
}
