//! Embedded web UI server for `qrnamap ui`.
//! Serves the pre-built `dist-site/` assets locally.

use std::path::Path;
use std::thread;
use std::time::Duration;
use tiny_http::{Header, Response, Server, StatusCode};

#[derive(rust_embed::RustEmbed)]
#[folder = "default-dist/"]
struct DistAssets;

pub fn run_ui(port: u16, host: &str, open_browser: bool) -> Result<(), Box<dyn std::error::Error>> {
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

        let local_file = Path::new("dist-site").join(&target_path);
        if local_file.is_file() {
            if let Ok(bytes) = std::fs::read(&local_file) {
                let mime = mime_type(&target_path);
                let header = Header::from_bytes(&b"Content-Type"[..], mime.as_bytes())
                    .unwrap_or_else(|_| Header::from_bytes(&b"Content-Type"[..], &b"application/octet-stream"[..]).unwrap());
                let mut response = Response::from_data(bytes);
                response.add_header(header);
                let _ = request.respond(response);
                continue;
            }
        }

        let file = DistAssets::get(&target_path).or_else(|| {
            if !target_path.contains('.') {
                DistAssets::get(&format!("{target_path}/index.html"))
            } else {
                None
            }
        });

        match file {
            Some(content) => {
                let mime = mime_type(&target_path);
                let header = Header::from_bytes(&b"Content-Type"[..], mime.as_bytes())
                    .unwrap_or_else(|_| Header::from_bytes(&b"Content-Type"[..], &b"application/octet-stream"[..]).unwrap());
                let mut response = Response::from_data(content.data.into_owned());
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
