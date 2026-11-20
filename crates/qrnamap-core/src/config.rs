//! Run configuration shared by the CLI and the web front end. Plain data, no I/O.

/// Library strandedness. Each tool maps it to its own flags in `toolcall`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strandedness {
    Unstranded,
    Forward,
    Reverse,
}

/// Quantification route.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Quant,
    Align,
}

/// Pseudo-alignment quantifier. Salmon is added in a later plan.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantifier {
    Kallisto,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RunConfig {
    pub mode: Mode,
    pub quantifier: Quantifier,
    pub strandedness: Strandedness,
    /// Read pairs per sample (default 20,000,000).
    pub reads_per_sample: u64,
    /// Fragment length mean and SD. Required for single-end kallisto.
    pub fragment_len: Option<(f64, f64)>,
    /// HISAT2 `--max-intronlen`. The default of 5000 is a chosen value.
    /// UNVERIFIED for fungal genomes. Revisit with data.
    pub max_intron: u32,
    pub threads: u32,
}

impl Default for RunConfig {
    fn default() -> Self {
        RunConfig {
            mode: Mode::Quant,
            quantifier: Quantifier::Kallisto,
            strandedness: Strandedness::Unstranded,
            reads_per_sample: 20_000_000,
            fragment_len: None,
            max_intron: 5000,
            threads: 1,
        }
    }
}
