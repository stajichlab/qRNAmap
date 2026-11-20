//! Native genome and annotation downloads (NCBI and FungiDB). Ensembl is not built.
//!
//! This crate does network and file I/O, so it is native only. Every download goes to a `.part` file in the output
//! directory. The file gets its final name only after the download is complete and (for NCBI) the size is checked.

mod fungidb;
mod ncbi;

pub use fungidb::fetch_fungidb;
pub use ncbi::fetch_ncbi;

use std::fmt;
use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

/// Name of the environment variable that holds the FungiDB API key.
pub const KEY_VARIABLE: &str = "EUPATHDB_KEY";

/// Paths of the downloaded files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fetched {
    pub fasta: PathBuf,
    /// `None` when the source has no annotation file for this genome.
    pub gff: Option<PathBuf>,
}

/// Where the FungiDB API key comes from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KeySource {
    /// The environment variable `EUPATHDB_KEY`.
    Env,
    /// A key file. Lines are `KEY` or `HOST KEY`; `#` starts a comment; `export EUPATHDB_KEY="..."` is accepted.
    /// The key file wins over the environment variable (the caller chooses which variant to pass).
    File(PathBuf),
    /// No key. Public files still download; a protected file gives `FetchError::NeedsKey`.
    None,
}

/// Download errors. No variant holds an API key.
#[derive(Debug)]
pub enum FetchError {
    /// The downloaded file has a different byte length than the catalog says.
    SizeMismatch {
        /// File name from the catalog.
        file: String,
        /// Length in the catalog.
        expected: u64,
        /// Length received.
        actual: u64,
    },
    /// The server answered 401 and no key was given.
    NeedsKey {
        /// The URL that needs the key.
        url: String,
    },
    /// An HTTP or connection failure, a rejected key, or an unusable server answer.
    Http {
        url: String,
        /// HTTP status when the server answered, else `None`.
        status: Option<u16>,
        message: String,
    },
    /// A local file problem (create, write, rename, key file).
    Io(io::Error),
    /// The API key is empty. `source` names where it came from (never the key).
    EmptyKey { source: String },
}

impl fmt::Display for FetchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FetchError::SizeMismatch { file, expected, actual } => {
                write!(f, "{file}: downloaded {actual} bytes, the catalog says {expected} bytes. The file was not kept")
            }
            FetchError::NeedsKey { url } => write!(
                f,
                "{url} needs an API key. Set the environment variable {KEY_VARIABLE} (get the key from your FungiDB account page, Web Services Access)"
            ),
            FetchError::Http { url, status: Some(code), message } => write!(f, "{url}: HTTP {code}: {message}"),
            FetchError::Http { url, status: None, message } => write!(f, "{url}: {message}"),
            FetchError::Io(e) => write!(f, "file error: {e}"),
            FetchError::EmptyKey { source } => write!(f, "the API key from {source} is empty"),
        }
    }
}

impl std::error::Error for FetchError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            FetchError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for FetchError {
    fn from(e: io::Error) -> Self {
        FetchError::Io(e)
    }
}

/// Percent-encode everything except the unreserved characters.
pub(crate) fn percent_encode(text: &str) -> String {
    let mut out = String::new();
    for b in text.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

pub(crate) fn agent_without_redirects() -> ureq::Agent {
    ureq::AgentBuilder::new()
        .redirects(0)
        .timeout_connect(std::time::Duration::from_secs(30))
        .timeout_read(std::time::Duration::from_secs(300))
        .build()
}

fn invalid_name(name: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, format!("unsafe output file name \"{name}\""))
}

/// Write `reader` to `out_dir/final_name` through `final_name.part`.
/// With `gzip`, the data is compressed on the way. `expected` is the required length of the received data in bytes
/// (before compression); on a mismatch the part file is removed and `SizeMismatch` is returned.
pub(crate) fn save(mut reader: impl Read, out_dir: &Path, final_name: &str, gzip: bool, expected: Option<u64>) -> Result<PathBuf, FetchError> {
    if final_name.is_empty() || final_name.contains(['/', '\\']) || final_name == "." || final_name == ".." {
        return Err(invalid_name(final_name).into());
    }
    fs::create_dir_all(out_dir)?;
    let final_path = out_dir.join(final_name);
    let part_path = out_dir.join(format!("{final_name}.part"));
    let result = (|| -> Result<u64, FetchError> {
        let file = BufWriter::new(File::create(&part_path)?);
        let received = if gzip {
            let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            let n = io::copy(&mut reader, &mut encoder)?;
            encoder.finish()?.flush()?;
            n
        } else {
            let mut file = file;
            let n = io::copy(&mut reader, &mut file)?;
            file.flush()?;
            n
        };
        if let Some(expected) = expected {
            if received != expected {
                return Err(FetchError::SizeMismatch { file: final_name.to_string(), expected, actual: received });
            }
        }
        fs::rename(&part_path, &final_path)?;
        Ok(received)
    })();
    match result {
        Ok(_) => Ok(final_path),
        Err(e) => {
            let _ = fs::remove_file(&part_path);
            Err(e)
        }
    }
}
