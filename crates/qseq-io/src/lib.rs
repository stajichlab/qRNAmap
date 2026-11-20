//! Sequence readers (FASTA, FASTQ) over in-memory bytes. No file or network I/O; builds for wasm32.

pub mod fasta;
pub mod fastq;

pub use fasta::read_fasta;
pub use fastq::{FastqChunker, FastqReader, DEFAULT_MAX_RECORD_BYTES};
