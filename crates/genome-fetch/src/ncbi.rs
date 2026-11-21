//! NCBI downloads. New code (not from qmutfind). The catalog decides the URLs, so http and https both work.

use crate::{save, FetchError, Fetched};
use genome_catalog::{FileRef, Genome};
use std::path::{Path, PathBuf};

fn download_one(agent: &ureq::Agent, f: &FileRef, out_dir: &Path) -> Result<PathBuf, FetchError> {
    let http_err = |status, message: String| FetchError::Http {
        url: f.url.clone(),
        status,
        message,
    };
    let response = match agent.get(&f.url).call() {
        Ok(r) => r,
        Err(ureq::Error::Status(code, _)) => {
            return Err(http_err(Some(code), "download failed".into()))
        }
        Err(e) => return Err(http_err(None, e.to_string())),
    };
    if response.status() != 200 {
        return Err(http_err(
            Some(response.status()),
            "unexpected status".into(),
        ));
    }
    // The catalog names are already `.gz`; the bytes are stored as received.
    save(
        response.into_reader(),
        out_dir,
        &f.name,
        false,
        Some(f.bytes),
    )
}

/// Download the FASTA and, when the catalog has one, the GFF3 of `g` into `out_dir`.
/// Files keep the catalog names. Each file must have exactly the catalog byte length.
/// A failed or short download leaves no file under the final name.
pub fn fetch_ncbi(g: &Genome, out_dir: &Path) -> Result<Fetched, FetchError> {
    // No key is involved, so ureq may follow redirects itself.
    let agent = ureq::AgentBuilder::new()
        .timeout_connect(std::time::Duration::from_secs(30))
        .timeout_read(std::time::Duration::from_secs(300))
        .build();
    let fasta = download_one(&agent, &g.fasta, out_dir)?;
    let gff = match &g.gff {
        Some(f) => Some(download_one(&agent, f, out_dir)?),
        None => None,
    };
    Ok(Fetched { fasta, gff })
}
