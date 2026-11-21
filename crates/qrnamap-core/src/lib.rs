//! qRNAmap core. No file or network I/O; builds for wasm32.

pub mod config;
pub mod counts;
pub mod de;
pub mod design;
pub mod genes;
pub mod saf;
pub mod sheet;
pub mod splice;
pub mod strand;
pub mod toolcall;
pub mod transcripts;

pub use config::{Mode, Quantifier, RunConfig, Strandedness};
pub use toolcall::{
    argv_native, feature_counts, hisat2_align, hisat2_build, kallisto_index, kallisto_quant,
    salmon_index, salmon_quant, ConfigError, Tool, ToolCall,
};

pub use counts::{
    assemble, gene_counts_from_abundance, gene_counts_from_salmon, parse_featurecounts,
    parse_featurecounts_summary, unlisted_genes, CountError, CountMatrix,
};
pub use design::{Contrast, Design, DesignError};
pub use genes::{gene_table, gene_table_tsv, GeneClass, GeneRow};
pub use saf::saf;
pub use sheet::{ColType, Sample, SampleSheet, SheetError};
pub use splice::splice_sites;
pub use strand::{assigned_fraction, StrandReport, STRAND_THRESHOLD};
pub use transcripts::{build_transcripts, tx2gene_tsv, write_fasta, BuildError, TxRecord};
