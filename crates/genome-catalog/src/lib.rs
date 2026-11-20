//! Catalog of reference genomes (FASTA and GFF3 URLs, sizes, annotation facts).
//!
//! This crate does no file or network I/O. The snapshot in `data/genomes.json` is
//! embedded at build time. `scripts/make_genome_list.py` regenerates it.

use serde::Deserialize;
use std::fmt;

/// The only schema version this crate reads.
pub const SCHEMA_VERSION: u32 = 1;

/// Genomes larger than this (in bases) get a size warning.
pub const SIZE_WARNING_BASES: u64 = 50_000_000;

const BUILTIN_JSON: &str = include_str!("../data/genomes.json");

/// A downloadable file.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct FileRef {
    pub name: String,
    pub url: String,
    pub bytes: u64,
}

/// One genome. Unknown JSON fields are ignored.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct Genome {
    pub accession: String,
    pub source: String,
    pub label: String,
    pub organism: String,
    #[serde(default)]
    pub strain: Option<String>,
    pub level: String,
    pub assembly: String,
    pub genome_size: u64,
    #[serde(default)]
    pub genes: Option<u32>,
    pub fasta: FileRef,
    /// `None` when NCBI has no GFF3 for the assembly.
    #[serde(default)]
    pub gff: Option<FileRef>,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
pub struct Catalog {
    pub schema_version: u32,
    pub generated: String,
    pub genomes: Vec<Genome>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CatalogError {
    /// Malformed JSON, or JSON that does not match the schema (including a missing `schema_version`).
    Json(String),
    /// `schema_version` is present but is not 1.
    UnknownSchemaVersion(u32),
}

impl fmt::Display for CatalogError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CatalogError::Json(msg) => write!(f, "invalid genome catalog JSON: {msg}"),
            CatalogError::UnknownSchemaVersion(v) => write!(
                f,
                "unknown genome catalog schema_version {v} (this build reads version {SCHEMA_VERSION})"
            ),
        }
    }
}

impl std::error::Error for CatalogError {}

/// Parse a catalog from JSON text.
pub fn parse(json: &str) -> Result<Catalog, CatalogError> {
    let catalog: Catalog =
        serde_json::from_str(json).map_err(|e| CatalogError::Json(e.to_string()))?;
    if catalog.schema_version != SCHEMA_VERSION {
        return Err(CatalogError::UnknownSchemaVersion(catalog.schema_version));
    }
    Ok(catalog)
}

/// The catalog embedded in this build.
///
/// Panics only if the embedded snapshot is invalid. The test `builtin_parses`
/// guards this, so a bad snapshot fails in CI and cannot reach a release.
pub fn builtin() -> Catalog {
    parse(BUILTIN_JSON).expect("embedded genomes.json must parse (guarded by test builtin_parses)")
}

impl Catalog {
    pub fn find(&self, accession: &str) -> Option<&Genome> {
        self.genomes.iter().find(|g| g.accession == accession)
    }

    /// Genomes that have a GFF3, which RNA-seq quantification needs.
    pub fn usable_for_rna(&self) -> Vec<&Genome> {
        self.genomes.iter().filter(|g| g.gff.is_some()).collect()
    }
}

/// A warning text when the genome is strictly larger than 50,000,000 bases.
pub fn size_warning(g: &Genome) -> Option<String> {
    if g.genome_size > SIZE_WARNING_BASES {
        Some(format!(
            "{} has a genome of {:.1} Mb, which is above the {} Mb size limit for a quick run.",
            g.organism,
            g.genome_size as f64 / 1_000_000.0,
            SIZE_WARNING_BASES / 1_000_000
        ))
    } else {
        None
    }
}
