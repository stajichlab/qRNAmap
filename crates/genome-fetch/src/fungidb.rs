//! FungiDB downloads. Moved and adapted from qmutfind `qmutfind-cli/src/fetch.rs`.
//!
//! The organism list is public. The file downloads need an API key. The key goes only to the host of the base URL,
//! never over plain http when the site is https, and it is never printed or put in an error.

use crate::{
    agent_without_redirects, percent_encode, save, FetchError, Fetched, KeySource, KEY_VARIABLE,
};
use std::fs;
use std::path::Path;

const MAX_REDIRECTS: usize = 5;

/// The host part of an http(s) URL, lower case, with the port when the URL has one.
/// Two hosts are the same only when this whole text is equal ("localhost:80" and "127.0.0.1:80" differ).
/// `None` when the URL has another form.
fn host_of(url: &str) -> Option<String> {
    let rest = url
        .strip_prefix("https://")
        .or_else(|| url.strip_prefix("http://"))?;
    Some(rest.split(['/', '?', '#']).next()?.to_ascii_lowercase())
}

/// The decision function for the key: send it to `url` only when the host equals the host of `base`, and not when
/// `base` is https but `url` is plain http.
fn key_allowed(base: &str, url: &str) -> bool {
    let same_host = matches!((host_of(base), host_of(url)), (Some(a), Some(b)) if a == b);
    let downgrade = base.starts_with("https://") && url.starts_with("http://");
    same_host && !downgrade
}

struct Organism {
    name: String,
    fasta: String,
    gff: String,
}

fn http_err(url: &str, status: Option<u16>, message: impl Into<String>) -> FetchError {
    FetchError::Http {
        url: url.to_string(),
        status,
        message: message.into(),
    }
}

/// Ask the site for all organisms with their file URLs. This request needs no key.
fn list_organisms(agent: &ureq::Agent, base: &str) -> Result<Vec<Organism>, FetchError> {
    let config = r#"{"attributes":["primary_key","URLGenomeFasta","URLgff"],"tables":[],"attributeFormat":"text"}"#;
    let url = format!(
        "{base}/fungidb/service/record-types/organism/searches/GenomeDataTypes/reports/standard?reportConfig={}",
        percent_encode(config)
    );
    let body = match agent.get(&url).call() {
        Ok(r) => r.into_string().map_err(FetchError::Io)?,
        Err(ureq::Error::Status(code, _)) => {
            return Err(http_err(&url, Some(code), "cannot list organisms"))
        }
        Err(e) => return Err(http_err(&url, None, format!("cannot list organisms: {e}"))),
    };
    let json: serde_json::Value = serde_json::from_str(&body)
        .map_err(|e| http_err(&url, None, format!("the organism list is not JSON: {e}")))?;
    let records = json["records"]
        .as_array()
        .ok_or_else(|| http_err(&url, None, "the organism list has no records"))?;
    Ok(records
        .iter()
        .filter_map(|r| {
            let a = &r["attributes"];
            Some(Organism {
                name: a["primary_key"].as_str()?.to_string(),
                fasta: a["URLGenomeFasta"].as_str()?.to_string(),
                gff: a["URLgff"].as_str()?.to_string(),
            })
        })
        .collect())
}

/// Pick the organism whose name equals the request (ignoring case). Otherwise list the names that contain it.
fn choose<'a>(
    organisms: &'a [Organism],
    wanted: &str,
    base: &str,
) -> Result<&'a Organism, FetchError> {
    let w = wanted
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_lowercase();
    if let Some(o) = organisms.iter().find(|o| o.name.to_lowercase() == w) {
        return Ok(o);
    }
    let close: Vec<&str> = organisms
        .iter()
        .filter(|o| o.name.to_lowercase().contains(&w))
        .map(|o| o.name.as_str())
        .collect();
    let message = if close.is_empty() {
        format!("no organism named \"{wanted}\" at this site")
    } else {
        format!(
            "no organism named exactly \"{wanted}\". Close matches:\n  {}",
            close.join("\n  ")
        )
    };
    Err(http_err(base, None, message))
}

/// One GET with manual redirects. The key is sent only while `key_allowed` is true for the current URL.
// `ureq::Error` is large (clippy::result_large_err). It is kept unboxed on purpose: it is ureq's own
// type, the one caller matches on its `Status` and `Transport` variants, and the function runs once
// per downloaded file, so the size of the `Err` value costs nothing measurable.
#[allow(clippy::result_large_err)]
fn get(
    agent: &ureq::Agent,
    url: &str,
    auth: Option<&str>,
    base: &str,
) -> Result<ureq::Response, ureq::Error> {
    let mut url = url.to_string();
    for _ in 0..=MAX_REDIRECTS {
        let mut request = agent.get(&url);
        if let Some(value) = auth {
            if key_allowed(base, &url) {
                request = request.set("Authorization", value);
            }
        }
        let response = request.call()?;
        if matches!(response.status(), 301 | 302 | 303 | 307 | 308) {
            match response.header("Location") {
                Some(next) if next.starts_with('/') => {
                    // A relative Location stays on the current host.
                    let scheme = if url.starts_with("https://") {
                        "https://"
                    } else {
                        "http://"
                    };
                    url = format!("{scheme}{}{next}", host_of(&url).unwrap_or_default());
                }
                Some(next) => url = next.to_string(),
                None => return Ok(response),
            }
        } else {
            return Ok(response);
        }
    }
    Err(ureq::Error::Status(
        310,
        ureq::Response::new(310, "Too many redirects", "").unwrap(),
    ))
}

fn download(
    agent: &ureq::Agent,
    url: &str,
    key: Option<&str>,
    key_source: &str,
    base: &str,
    out_dir: &Path,
) -> Result<std::path::PathBuf, FetchError> {
    // The organism list gives `/common/downloads/` URLs. They answer 404 on the real site (2026-11-12). The same
    // file is served, with the key, under `/a/service/raw-files/`.
    let url = url.replace("/common/downloads/", "/a/service/raw-files/");
    // Never send the key over plain http when the site is https.
    let url = if base.starts_with("https://") && url.starts_with("http://") {
        format!("https://{}", &url["http://".len()..])
    } else {
        url
    };
    let name = url
        .split(['?', '#'])
        .next()
        .and_then(|u| u.rsplit('/').next())
        .filter(|n| !n.is_empty())
        .ok_or_else(|| http_err(&url, None, "cannot find a file name in the URL"))?;
    let already_gz = name.ends_with(".gz");
    let final_name = if already_gz {
        name.to_string()
    } else {
        format!("{name}.gz")
    };

    let bearer = key.map(|k| format!("Bearer {k}"));
    let response = match get(agent, &url, bearer.as_deref(), base) {
        Ok(r) => r,
        Err(ureq::Error::Status(401, _)) if key.is_none() => {
            return Err(FetchError::NeedsKey { url })
        }
        Err(ureq::Error::Status(401, _)) => {
            return Err(http_err(
                &url,
                Some(401),
                format!("the server rejected the key from {key_source}"),
            ))
        }
        Err(ureq::Error::Status(code, _)) => {
            return Err(http_err(&url, Some(code), "download failed"))
        }
        // ureq's text for a transport error can hold the whole request header line (kind BadHeader gives
        // "invalid header 'Authorization: Bearer <key>'"). So with a key, use a fixed text and only the error kind.
        Err(e) if key.is_some() => {
            return Err(http_err(
                &url,
                None,
                format!(
                "network error ({:?}) on a request that carries the API key; details are withheld",
                e.kind()
            ),
            ))
        }
        Err(e) => return Err(http_err(&url, None, e.to_string())),
    };
    if response.status() != 200 {
        return Err(http_err(&url, Some(response.status()), "unexpected status"));
    }
    save(
        response.into_reader(),
        out_dir,
        &final_name,
        !already_gz,
        None,
    )
}

/// The value of a line like `EUPATHDB_KEY=abc`, `export EUPATHDB_KEY="abc"` or `EUPATHDB_KEY='abc'`.
/// `None` for any other line. Only the name `EUPATHDB_KEY` counts, so a key that ends in `=` is not cut.
fn shell_assignment_value(line: &str) -> Option<&str> {
    let rest = line
        .strip_prefix("export ")
        .map(str::trim_start)
        .unwrap_or(line);
    let value = rest.strip_prefix(KEY_VARIABLE)?.strip_prefix('=')?.trim();
    let value = value
        .strip_prefix('"')
        .and_then(|v| v.strip_suffix('"'))
        .or_else(|| value.strip_prefix('\'').and_then(|v| v.strip_suffix('\'')))
        .unwrap_or(value);
    Some(value)
}

/// Parse key file text. Each line is `KEY` (any site) or `HOST KEY` (that site only). A line for the host wins over a
/// line without a host. Blank lines and lines that start with `#` are ignored.
fn key_from_text(text: &str, host: &str) -> Option<String> {
    let mut default = None;
    for line in text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        // `EUPATHDB_KEY=""` (empty after unquoting) is an empty key: this line gives no key (R29-14).
        let line = match shell_assignment_value(line) {
            Some("") => continue,
            Some(v) => v,
            None => line,
        };
        let mut words = line.split_whitespace();
        match (words.next(), words.next(), words.next()) {
            (Some(key), None, _) => default = default.or(Some(key)),
            (Some(h), Some(key), None) if h.eq_ignore_ascii_case(host) => {
                return Some(key.to_string())
            }
            _ => {}
        }
    }
    default.map(str::to_string)
}

/// Read the key for `host` from a key file. Returns the key and whether other users can read the file.
fn key_from_file(path: &Path, host: &str) -> Result<(String, bool), FetchError> {
    use std::os::unix::fs::PermissionsExt;
    let meta = fs::metadata(path)?;
    let open = meta.permissions().mode() & 0o077 != 0;
    let text = fs::read_to_string(path)?;
    match key_from_text(&text, host) {
        Some(k) => Ok((k, open)),
        None => Err(FetchError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("no key for {host} in key file {}", path.display()),
        ))),
    }
}

/// Refuse a key with any byte outside 0x21-0x7E (space, control characters such as a newline, non-ASCII).
/// Such a key cannot be a valid header value, and an HTTP library may print the header line in its error.
/// The error names the source and never holds the key or any part of it. It is `FetchError::Io` (InvalidData).
/// An empty key is `FetchError::EmptyKey` (R29-14).
fn validate_key(key: &str, described_source: &str) -> Result<(), FetchError> {
    if key.is_empty() {
        return Err(FetchError::EmptyKey {
            source: described_source.to_string(),
        });
    }
    if key.bytes().all(|b| (0x21..=0x7E).contains(&b)) {
        Ok(())
    } else {
        Err(FetchError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("the key from {described_source} holds characters that are not allowed (only printable ASCII without spaces is accepted)"),
        )))
    }
}

/// Resolve the key. Returns the key (if any) and a description of its source for messages.
fn resolve_key(
    source: &KeySource,
    base_host: &str,
) -> Result<(Option<String>, String), FetchError> {
    match source {
        KeySource::File(path) => {
            let (key, open) = key_from_file(path, base_host)?;
            let described = format!("key file {}", path.display());
            validate_key(&key, &described)?;
            if open {
                eprintln!("genome-fetch: warning: key file {} is readable by other users. Run: chmod 600 {}", path.display(), path.display());
            }
            Ok((Some(key), described))
        }
        KeySource::Env => {
            let key = std::env::var(KEY_VARIABLE)
                .ok()
                .map(|k| k.trim().to_string())
                .filter(|k| !k.is_empty());
            let described = format!("the environment variable {KEY_VARIABLE}");
            if let Some(k) = &key {
                validate_key(k, &described)?;
            }
            Ok((key, described))
        }
        KeySource::None => Ok((None, "no key".to_string())),
    }
}

/// Download the FASTA and GFF of one FungiDB organism into `out_dir`. Files that do not end in `.gz` are compressed
/// with gzip and get a `.gz` name.
pub fn fetch_fungidb(
    organism: &str,
    key: &KeySource,
    base_url: &str,
    out_dir: &Path,
) -> Result<Fetched, FetchError> {
    let base = base_url.trim_end_matches('/');
    let base_host = host_of(base).ok_or_else(|| {
        http_err(
            base,
            None,
            "the base URL must start with http:// or https://",
        )
    })?;
    let (key, key_source) = resolve_key(key, &base_host)?;
    let agent = agent_without_redirects();
    let organisms = list_organisms(&agent, base)?;
    let o = choose(&organisms, organism, base)?;
    let fasta = download(&agent, &o.fasta, key.as_deref(), &key_source, base, out_dir)?;
    let gff = download(&agent, &o.gff, key.as_deref(), &key_source, base, out_dir)?;
    Ok(Fetched {
        fasta,
        gff: Some(gff),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn never_sends_key_over_http_when_site_is_https() {
        let base = "https://fungidb.org";
        assert!(key_allowed(base, "https://fungidb.org/a/x.gz"));
        assert!(
            !key_allowed(base, "http://fungidb.org/a/x.gz"),
            "downgrade to http"
        );
        assert!(!key_allowed(base, "https://other.org/a/x.gz"), "other host");
        assert!(
            !key_allowed(base, "https://fungidb.org:8443/a/x.gz"),
            "other port is another host"
        );
        // An http site (a local test server) may use http.
        assert!(key_allowed("http://127.0.0.1:9", "http://127.0.0.1:9/x"));
        assert!(!key_allowed("http://127.0.0.1:9", "http://localhost:9/x"));
        assert!(!key_allowed(base, "ftp://fungidb.org/x"));
    }

    #[test]
    fn key_file_host_line_beats_default_line() {
        let text = "# comment\nDEFAULTKEY\nfungidb.org HOSTKEY\nother.org OTHERKEY\n";
        assert_eq!(
            key_from_text(text, "fungidb.org").as_deref(),
            Some("HOSTKEY")
        );
        assert_eq!(
            key_from_text(text, "FungiDB.org").as_deref(),
            Some("HOSTKEY")
        );
        assert_eq!(
            key_from_text(text, "third.org").as_deref(),
            Some("DEFAULTKEY")
        );
        assert_eq!(
            key_from_text("export EUPATHDB_KEY=\"abc==\"\n", "x").as_deref(),
            Some("abc==")
        );
        assert_eq!(key_from_text("# only a comment\n", "x"), None);
    }

    #[test]
    fn validate_key_rejects_unsafe_bytes_without_leak() {
        let source = "the environment variable EUPATHDB_KEY";
        let cases = [
            "SENTINEL-A\nSENTINEL-B",
            "SENTINEL-A\rSENTINEL-B",
            "SENTINEL-A\tSENTINEL-B",
            "SENTINEL-A SENTINEL-B",
            "a b",
            "SENTINEL-A\u{7f}SENTINEL-B",
            "SENTINEL-\u{e9}",
        ];
        for key in cases {
            let err = validate_key(key, source).expect_err(key);
            let (d, g) = (err.to_string(), format!("{err:?}"));
            assert!(d.contains(source) && d.contains("not allowed"), "{d}");
            for part in ["SENTINEL-A", "SENTINEL-B", "SENTINEL-\u{e9}", "SENTINEL"] {
                assert!(
                    !d.contains(part) && !g.contains(part),
                    "leak of {part} for {key:?}: {d} / {g}"
                );
            }
        }
        // "a b" holds no sentinel; check its parts are absent too.
        let err = validate_key("a b", source).unwrap_err();
        assert!(!err.to_string().contains("a b") && !format!("{err:?}").contains("a b"));
    }

    #[test]
    fn validate_key_accepts_printable_ascii() {
        assert!(validate_key("AbC-123_x.y/z=+", "x").is_ok());
        assert!(
            validate_key("!~", "x").is_ok(),
            "0x21 and 0x7E are the limits"
        );
    }

    #[test]
    fn validate_key_rejects_empty_key_with_named_error() {
        // R29-14: an empty key is a named error, not Io.
        match validate_key("", "key file /k") {
            Err(FetchError::EmptyKey { source }) => assert_eq!(source, "key file /k"),
            other => panic!("expected EmptyKey, got {other:?}"),
        }
        let e = validate_key("", "key file /k").unwrap_err().to_string();
        assert!(e.contains("empty") && e.contains("key file /k"), "{e}");
    }

    #[test]
    fn empty_quoted_assignment_is_no_key() {
        // R29-14: `EUPATHDB_KEY=""` is an empty key, so the line gives no key.
        for text in [
            "EUPATHDB_KEY=\"\"\n",
            "export EUPATHDB_KEY=''\n",
            "EUPATHDB_KEY=\n",
        ] {
            assert_eq!(key_from_text(text, "fungidb.org"), None, "{text:?}");
        }
        assert_eq!(
            key_from_text("EUPATHDB_KEY=\"\"\nREALKEY\n", "x").as_deref(),
            Some("REALKEY")
        );
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("k");
        fs::write(&p, "EUPATHDB_KEY=\"\"\n").unwrap();
        let e = key_from_file(&p, "fungidb.org").unwrap_err().to_string();
        assert!(e.contains("no key"), "{e}");
    }

    #[test]
    fn key_file_world_readable_is_flagged() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("k");
        fs::write(&p, "KEY1\n").unwrap();
        fs::set_permissions(&p, fs::Permissions::from_mode(0o644)).unwrap();
        assert_eq!(key_from_file(&p, "h").unwrap(), ("KEY1".to_string(), true));
        fs::set_permissions(&p, fs::Permissions::from_mode(0o600)).unwrap();
        assert_eq!(key_from_file(&p, "h").unwrap(), ("KEY1".to_string(), false));
        // The warning text itself goes to stderr in resolve_key; its capture is not tested here.
    }
}
