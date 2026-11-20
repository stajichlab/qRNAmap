//! qRNAmap core. No file or network I/O; builds for wasm32.

pub mod design;
pub mod saf;
pub mod sheet;
pub mod splice;
pub mod transcripts;

pub use design::{Contrast, Design, DesignError};
pub use saf::saf;
pub use sheet::{ColType, Sample, SampleSheet, SheetError};
pub use splice::splice_sites;
pub use transcripts::{build_transcripts, tx2gene_tsv, write_fasta, BuildError, TxRecord};
