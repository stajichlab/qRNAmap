//! `qrnamap fetch`: download a genome FASTA and its GFF3 annotation.
//!
//! * `ncbi`: `--id` is an accession from the genome catalog (`genome_catalog::builtin()`). The
//!   files keep the catalog names and are checked against the catalog byte lengths.
//! * `fungidb`: `--id` is an organism name. The API key comes from `--eupathdb-key-file` or from
//!   `$EUPATHDB_KEY`; the key file wins when both are given (spec 5). An empty key is a named
//!   error before any request. The key is never printed.

use std::fmt;
use std::path::{Path, PathBuf};

use genome_catalog::{size_warning, Catalog, CatalogError};
use genome_fetch::{fetch_fungidb, fetch_ncbi, FetchError, Fetched, KeySource, KEY_VARIABLE};

/// The FungiDB site. `--base-url` replaces it (used by the tests with a local server).
pub const FUNGIDB_BASE_URL: &str = "https://fungidb.org";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Ncbi,
    Fungidb,
}

#[derive(Debug, Clone)]
pub struct FetchArgs {
    pub source: Source,
    /// Catalog accession (ncbi) or organism name (fungidb).
    pub id: String,
    pub out: PathBuf,
    pub key_file: Option<PathBuf>,
    /// FungiDB only. `None` means `FUNGIDB_BASE_URL`.
    pub base_url: Option<String>,
    /// A catalog file instead of the built-in one (test-only option `--catalog`).
    pub catalog: Option<PathBuf>,
}

/// What was fetched, and the warnings for the user.
#[derive(Debug, Clone)]
pub struct FetchReport {
    pub fetched: Fetched,
    pub warnings: Vec<String>,
}

#[derive(Debug)]
pub enum FetchCmdError {
    /// The accession is not in the genome catalog.
    UnknownAccession(String),
    /// The API key from this source is empty (`what` names the source, never the key).
    EmptyKey {
        what: String,
    },
    /// `--base-url` was given for a source that does not use it.
    BaseUrlNotUsed,
    /// The `--catalog` file cannot be read or parsed.
    Catalog {
        path: PathBuf,
        detail: String,
    },
    /// The key file cannot be read.
    KeyFile {
        path: PathBuf,
        source: std::io::Error,
    },
    Fetch(FetchError),
}

impl fmt::Display for FetchCmdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FetchCmdError::UnknownAccession(a) => write!(
                f,
                "accession {a:?} is not in the genome catalog built into qrnamap; only catalog genomes can be fetched from NCBI"
            ),
            FetchCmdError::EmptyKey { what } => write!(f, "the FungiDB API key from {what} is empty"),
            FetchCmdError::BaseUrlNotUsed => write!(f, "--base-url applies to --source fungidb only"),
            FetchCmdError::Catalog { path, detail } => write!(f, "catalog {}: {detail}", path.display()),
            FetchCmdError::KeyFile { path, source } => write!(f, "key file {}: {source}", path.display()),
            FetchCmdError::Fetch(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for FetchCmdError {}

/// True when the key file text holds no key: no line other than blanks and `#` comments, or an
/// `EUPATHDB_KEY=` assignment with an empty value. The genome-fetch parser would otherwise read
/// the assignment line itself as the key.
fn key_file_is_empty(text: &str) -> bool {
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    if lines.is_empty() {
        return true;
    }
    lines.iter().any(|l| {
        let rest = l.strip_prefix("export ").map(str::trim_start).unwrap_or(l);
        match rest
            .strip_prefix(KEY_VARIABLE)
            .and_then(|r| r.strip_prefix('='))
        {
            Some(v) => {
                let v = v.trim();
                v.is_empty() || v == "\"\"" || v == "''"
            }
            None => false,
        }
    })
}

/// Pick exactly one key source. The key file wins over the environment variable. An empty key
/// (empty key file, empty assignment in it, or an empty or blank `$EUPATHDB_KEY`) is an error.
/// `env_value` is the value of `$EUPATHDB_KEY`, `None` when it is not set.
pub fn choose_key(
    key_file: Option<&Path>,
    env_value: Option<&str>,
) -> Result<KeySource, FetchCmdError> {
    if let Some(p) = key_file {
        let text = std::fs::read_to_string(p).map_err(|source| FetchCmdError::KeyFile {
            path: p.to_path_buf(),
            source,
        })?;
        if key_file_is_empty(&text) {
            return Err(FetchCmdError::EmptyKey {
                what: format!("key file {}", p.display()),
            });
        }
        return Ok(KeySource::File(p.to_path_buf()));
    }
    match env_value {
        Some(v) if v.trim().is_empty() => Err(FetchCmdError::EmptyKey {
            what: format!("the environment variable {KEY_VARIABLE}"),
        }),
        Some(_) => Ok(KeySource::Env),
        None => Ok(KeySource::None),
    }
}

fn load_catalog(path: Option<&Path>) -> Result<Catalog, FetchCmdError> {
    match path {
        None => Ok(genome_catalog::builtin()),
        Some(p) => {
            let text = std::fs::read_to_string(p).map_err(|e| FetchCmdError::Catalog {
                path: p.to_path_buf(),
                detail: e.to_string(),
            })?;
            genome_catalog::parse(&text).map_err(|e: CatalogError| FetchCmdError::Catalog {
                path: p.to_path_buf(),
                detail: e.to_string(),
            })
        }
    }
}

/// Run one fetch. All checks (catalog lookup, key choice) happen before any request.
pub fn run_fetch(args: &FetchArgs) -> Result<FetchReport, FetchCmdError> {
    match args.source {
        Source::Ncbi => {
            if args.base_url.is_some() {
                return Err(FetchCmdError::BaseUrlNotUsed);
            }
            let catalog = load_catalog(args.catalog.as_deref())?;
            let g = catalog
                .find(&args.id)
                .ok_or_else(|| FetchCmdError::UnknownAccession(args.id.clone()))?;
            let mut warnings = Vec::new();
            if g.gff.is_none() {
                warnings.push(format!(
                    "{} ({}) has no GFF3 annotation in the catalog: only the FASTA is fetched, and it cannot be used for RNA-seq counting",
                    g.accession, g.organism
                ));
            }
            if let Some(w) = size_warning(g) {
                warnings.push(w);
            }
            let fetched = fetch_ncbi(g, &args.out).map_err(FetchCmdError::Fetch)?;
            Ok(FetchReport { fetched, warnings })
        }
        Source::Fungidb => {
            let env = std::env::var(KEY_VARIABLE).ok();
            let key = choose_key(args.key_file.as_deref(), env.as_deref())?;
            let base = args.base_url.as_deref().unwrap_or(FUNGIDB_BASE_URL);
            let fetched =
                fetch_fungidb(&args.id, &key, base, &args.out).map_err(FetchCmdError::Fetch)?;
            Ok(FetchReport {
                fetched,
                warnings: Vec::new(),
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn key_file_wins_and_env_is_fallback() {
        let d = tempfile::tempdir().unwrap();
        let kf = d.path().join("k");
        std::fs::write(&kf, "KEY\n").unwrap();
        assert_eq!(
            choose_key(Some(&kf), Some("ENV")).unwrap(),
            KeySource::File(kf.clone())
        );
        assert_eq!(choose_key(None, Some("ENV")).unwrap(), KeySource::Env);
        assert_eq!(choose_key(None, None).unwrap(), KeySource::None);
    }

    #[test]
    fn empty_keys_are_errors() {
        assert!(matches!(
            choose_key(None, Some("")),
            Err(FetchCmdError::EmptyKey { .. })
        ));
        assert!(matches!(
            choose_key(None, Some(" \t")),
            Err(FetchCmdError::EmptyKey { .. })
        ));
        for body in [
            "",
            "\n# c\n",
            "EUPATHDB_KEY=\n",
            "export EUPATHDB_KEY=''\n",
            "EUPATHDB_KEY=\"\"\n",
        ] {
            assert!(key_file_is_empty(body), "{body:?}");
        }
        for body in ["K\n", "fungidb.org K\n", "EUPATHDB_KEY=abc\n"] {
            assert!(!key_file_is_empty(body), "{body:?}");
        }
    }
}
