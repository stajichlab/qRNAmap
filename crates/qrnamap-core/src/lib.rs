//! qRNAmap core. No file or network I/O; builds for wasm32.

pub mod config;
pub mod design;
pub mod saf;
pub mod sheet;
pub mod splice;
pub mod toolcall;
pub mod transcripts;

pub use config::{Mode, Quantifier, RunConfig, Strandedness};
pub use toolcall::{
    argv_native, feature_counts, hisat2_align, hisat2_build, kallisto_index, kallisto_quant,
    ConfigError, Tool, ToolCall,
};

pub use design::{Contrast, Design, DesignError};
pub use saf::saf;
pub use sheet::{ColType, Sample, SampleSheet, SheetError};
pub use splice::splice_sites;
pub use transcripts::{build_transcripts, tx2gene_tsv, write_fasta, BuildError, TxRecord};
