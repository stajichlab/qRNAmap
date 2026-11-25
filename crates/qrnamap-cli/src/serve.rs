//! Web UI server for `qrnamap ui`.
//! Serves the built `dist-site/` assets locally.

use std::path::Path;
use std::thread;
use std::time::Duration;
use tiny_http::{Header, Response, Server, StatusCode};

/// The built web assets.
///
/// `dist-site/` is a build artifact (`scripts/build-site.sh`), so it is missing in a fresh
/// checkout. `allow_missing` keeps this crate compiling without it; the assets are embedded when
/// the site was built before the crate was compiled.
///
/// The reads of the test sets (`testsets/reads/`) are left out: they are 216 MB, and the page
/// falls back to the ENA stream for a test set whose reads are not there. The site on disk still
/// has them, so `qrnamap ui` in a checkout serves them from there.
#[derive(rust_embed::RustEmbed)]
#[folder = "../../dist-site/"]
#[allow_missing = true]
#[exclude = "testsets/reads/**"]
struct DistAssets;

/// The file `path` of the site next to the working directory, read from disk. This is the site
/// that `scripts/build-site.sh` wrote, so a site built after the crate was compiled is served
/// without a rebuild. `None` when there is no such file.
fn local_file(path: &str) -> Option<Vec<u8>> {
    std::fs::read(Path::new("dist-site").join(path)).ok()
}

/// True when the web assets can be served: embedded in the binary, or in `dist-site/` next to the
/// working directory.
fn assets_present() -> bool {
    DistAssets::iter().next().is_some() || local_file("analyze/index.html").is_some()
}

pub fn run_ui(port: u16, host: &str, open_browser: bool) -> Result<(), Box<dyn std::error::Error>> {
    if !assets_present() {
        return Err(
            "the web assets are not built: run `pixi run build-wasm-web` and \
             `bash scripts/build-site.sh` from the repository root, then rebuild qrnamap"
                .into(),
        );
    }

    let addr = format!("{host}:{port}");
    let server = Server::http(&addr)
        .map_err(|e| format!("Failed to bind local server on {addr}: {e}"))?;

    let url = format!("http://{addr}/analyze/");
    eprintln!("===========================================================");
    eprintln!("  qRNAmap local server running at: {url}");
    eprintln!("  Serving embedded web application.");
    eprintln!("  Press Ctrl+C to stop.");
    eprintln!("===========================================================");

    if open_browser {
        let open_url = url.clone();
        thread::spawn(move || {
            thread::sleep(Duration::from_millis(250));
            let _ = open::that(&open_url);
        });
    }

    for request in server.incoming_requests() {
        let raw_url = request.url().split('?').next().unwrap_or("/");
        let clean_path = if raw_url == "/" || raw_url.is_empty() {
            "/analyze/"
        } else {
            raw_url
        };

        let trimmed = clean_path.trim_start_matches('/');
        let target_path = if trimmed.is_empty() {
            "analyze/index.html".to_string()
        } else if trimmed.ends_with('/') {
            format!("{trimmed}index.html")
        } else {
            trimmed.to_string()
        };

        // The embedded assets come first: they are what the binary ships. The site next to the
        // working directory is the fallback for a binary built without a site.
        let file = DistAssets::get(&target_path)
            .or_else(|| {
                if !target_path.contains('.') {
                    DistAssets::get(&format!("{target_path}/index.html"))
                } else {
                    None
                }
            })
            .map(|content| content.data.into_owned())
            .or_else(|| local_file(&target_path))
            .or_else(|| {
                if !target_path.contains('.') {
                    local_file(&format!("{target_path}/index.html"))
                } else {
                    None
                }
            });

        match file {
            Some(content) => {
                let mime = mime_type(&target_path);
                let header = Header::from_bytes(&b"Content-Type"[..], mime.as_bytes())
                    .unwrap_or_else(|_| Header::from_bytes(&b"Content-Type"[..], &b"application/octet-stream"[..]).unwrap());
                let mut response = Response::from_data(content);
                response.add_header(header);
                let _ = request.respond(response);
            }
            None => {
                let response = Response::from_string("404 Not Found")
                    .with_status_code(StatusCode(404));
                let _ = request.respond(response);
            }
        }
    }

    Ok(())
}

fn mime_type(path: &str) -> &'static str {
    let ext = Path::new(path)
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    match ext {
        "html" => "text/html; charset=utf-8",
        "js" | "mjs" => "application/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "json" => "application/json",
        "wasm" => "application/wasm",
        "svg" => "image/svg+xml",
        "png" => "image/png",
        "jpg" | "jpeg" => "image/jpeg",
        "gz" => "application/gzip",
        "tsv" => "text/tab-separated-values; charset=utf-8",
        "txt" => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}
