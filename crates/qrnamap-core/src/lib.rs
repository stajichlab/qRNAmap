//! qRNAmap core. No file or network I/O; builds for wasm32.

pub mod saf;
pub mod splice;
pub mod transcripts;

pub use saf::saf;
pub use splice::splice_sites;
pub use transcripts::{build_transcripts, tx2gene_tsv, write_fasta, BuildError, TxRecord};
