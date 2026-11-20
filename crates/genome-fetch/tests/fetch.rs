mod common;

use common::{Response, Server};
use genome_catalog::{FileRef, Genome};
use genome_fetch::{fetch_fungidb, fetch_ncbi, FetchError, KeySource};
use std::fs;

const FASTA: &[u8] = b"fasta-gz-bytes-0123456789";
const GFF: &[u8] = b"gff-gz-bytes";

fn genome(base: &str, fasta_bytes: u64, with_gff: bool) -> Genome {
    Genome {
        accession: "GCF_000000001.1".into(),
        source: "ncbi".into(),
        label: "test".into(),
        organism: "Testus fungus".into(),
        strain: None,
        level: "Complete".into(),
        assembly: "ASM1".into(),
        genome_size: 1000,
        genes: None,
        fasta: FileRef { name: "t_genomic.fna.gz".into(), url: format!("{base}/t_genomic.fna.gz"), bytes: fasta_bytes },
        gff: with_gff.then(|| FileRef { name: "t_genomic.gff.gz".into(), url: format!("{base}/t_genomic.gff.gz"), bytes: GFF.len() as u64 }),
    }
}

fn ncbi_server() -> Server {
    Server::start(|req, _| match req.path.as_str() {
        "/t_genomic.fna.gz" => Response::ok(FASTA),
        "/t_genomic.gff.gz" => Response::ok(GFF),
        _ => Response::status(404),
    })
}

#[test]
fn ncbi_download_writes_files() {
    let s = ncbi_server();
    let dir = tempfile::tempdir().unwrap();
    let g = genome(&s.base, FASTA.len() as u64, true);
    let out = fetch_ncbi(&g, dir.path()).unwrap();
    assert_eq!(out.fasta, dir.path().join("t_genomic.fna.gz"));
    assert_eq!(fs::read(&out.fasta).unwrap(), FASTA);
    assert_eq!(fs::read(out.gff.unwrap()).unwrap(), GFF);
    let names: Vec<_> = fs::read_dir(dir.path()).unwrap().map(|e| e.unwrap().file_name()).collect();
    assert_eq!(names.len(), 2, "no temporary files left: {names:?}");
}

#[test]
fn ncbi_download_checks_size() {
    let s = ncbi_server();
    let dir = tempfile::tempdir().unwrap();
    let g = genome(&s.base, FASTA.len() as u64 + 10, false);
    let err = fetch_ncbi(&g, dir.path()).unwrap_err();
    match err {
        FetchError::SizeMismatch { expected, actual, .. } => {
            assert_eq!(expected, FASTA.len() as u64 + 10);
            assert_eq!(actual, FASTA.len() as u64);
        }
        other => panic!("wrong error: {other}"),
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0, "no file may remain");
}

/// FungiDB organism list JSON pointing at the given FASTA and GFF URLs.
fn organism_json(fasta: &str, gff: &str) -> Vec<u8> {
    format!(r#"{{"records":[{{"attributes":{{"primary_key":"Testus fungus","URLGenomeFasta":"{fasta}","URLgff":"{gff}"}}}}]}}"#).into_bytes()
}

fn write_key_file(dir: &std::path::Path, text: &str) -> KeySource {
    use std::os::unix::fs::PermissionsExt;
    let p = dir.join("key.txt");
    fs::write(&p, text).unwrap();
    fs::set_permissions(&p, fs::Permissions::from_mode(0o600)).unwrap();
    KeySource::File(p)
}

#[test]
fn fungidb_key_sent_only_to_key_host() {
    // Host comparison: the host is the text between "://" and the next "/", with the port, in lower case.
    // "127.0.0.1:P1" (key host) and "localhost:P2" differ, so the key must not go to the second server.
    let second = Server::start(|_, _| Response::ok(b"data-gz".as_slice()));
    let second_url = format!("http://localhost:{}/file.fa.gz", second.port);
    let first = Server::start(move |req, base| {
        if req.path.starts_with("/fungidb/service/") {
            Response::ok(organism_json(&format!("{base}/a/service/raw-files/x.fa.gz"), &format!("{base}/a/service/raw-files/x.gff.gz")))
        } else if req.path == "/a/service/raw-files/x.fa.gz" {
            Response::redirect(&second_url)
        } else {
            Response::ok(b"gff-gz".as_slice())
        }
    });
    let dir = tempfile::tempdir().unwrap();
    let key = write_key_file(dir.path(), "SENTINEL-KEY-1\n");
    let out = dir.path().join("out");
    fetch_fungidb("Testus fungus", &key, &first.base, &out).unwrap();
    let on_first: Vec<_> = first.requests().into_iter().filter(|r| r.path.contains("x.fa.gz")).collect();
    assert_eq!(on_first.len(), 1);
    assert_eq!(on_first[0].header("authorization"), Some("Bearer SENTINEL-KEY-1"));
    let seen = second.requests();
    assert_eq!(seen.len(), 1);
    assert!(seen[0].header("authorization").is_none(), "key leaked to second host");
}

#[test]
fn fungidb_401_without_key_names_variable() {
    let s = Server::start(|req, base| {
        if req.path.starts_with("/fungidb/service/") {
            Response::ok(organism_json(&format!("{base}/a/service/raw-files/x.fa.gz"), &format!("{base}/a/service/raw-files/x.gff.gz")))
        } else {
            Response::status(401)
        }
    });
    let dir = tempfile::tempdir().unwrap();
    let err = fetch_fungidb("Testus fungus", &KeySource::None, &s.base, dir.path()).unwrap_err();
    assert!(matches!(err, FetchError::NeedsKey { .. }));
    let text = err.to_string();
    assert!(text.contains("EUPATHDB_KEY"), "{text}");
    assert!(text.contains("FungiDB account"), "{text}");
}

#[test]
fn key_never_in_error_text() {
    let sentinel = "SENTINEL-SECRET-KEY-XYZ";
    let s = Server::start(|req, base| {
        if req.path.starts_with("/fungidb/service/") {
            Response::ok(organism_json(&format!("{base}/a/service/raw-files/x.fa.gz"), &format!("{base}/a/service/raw-files/x.gff.gz")))
        } else {
            Response::status(401)
        }
    });
    let dir = tempfile::tempdir().unwrap();
    let key = write_key_file(dir.path(), &format!("{sentinel}\n"));
    let err = fetch_fungidb("Testus fungus", &key, &s.base, &dir.path().join("o")).unwrap_err();
    let (d, g) = (err.to_string(), format!("{err:?}"));
    assert!(d.contains("rejected"), "{d}");
    assert!(!d.contains(sentinel) && !g.contains(sentinel), "{d} {g}");
    // Another error path: an organism that is not in the list.
    let err = fetch_fungidb("No such organism", &key, &s.base, &dir.path().join("o")).unwrap_err();
    let (d, g) = (err.to_string(), format!("{err:?}"));
    assert!(!d.contains(sentinel) && !g.contains(sentinel));
    assert!(!format!("{key:?}").contains(sentinel));
}

fn list_server(file_response: fn() -> Response) -> Server {
    Server::start(move |req, base| {
        if req.path.starts_with("/fungidb/service/") {
            Response::ok(organism_json(&format!("{base}/a/service/raw-files/x.fa.gz"), &format!("{base}/a/service/raw-files/x.gff.gz")))
        } else {
            file_response()
        }
    })
}

fn assert_no_leak(err: &FetchError, parts: &[&str]) {
    let (d, g) = (err.to_string(), format!("{err:?}"));
    for p in parts {
        assert!(!d.contains(p) && !g.contains(p), "key part {p} leaked: {d} / {g}");
    }
}

#[test]
fn key_file_second_line_is_ignored_and_401_error_has_no_leak() {
    // The file parser takes the first one-word line as the key, so the key is `SENTINEL-A` (valid) and the
    // server answers 401. This checks only that the "key rejected" error does not hold either word.
    let s = list_server(|| Response::status(401));
    let dir = tempfile::tempdir().unwrap();
    let key = write_key_file(dir.path(), "SENTINEL-A\nSENTINEL-B\n");
    let err = fetch_fungidb("Testus fungus", &key, &s.base, &dir.path().join("o")).unwrap_err();
    assert!(err.to_string().contains("rejected"), "{err}");
    assert_no_leak(&err, &["SENTINEL-A", "SENTINEL-B"]);
}

#[test]
fn key_file_with_curly_quotes_is_rejected_without_leak() {
    let s = list_server(|| Response::status(401));
    let dir = tempfile::tempdir().unwrap();
    // Curly quotes are not removed by the parser, so the key holds non-ASCII bytes and must be refused before any request.
    let key = write_key_file(dir.path(), "EUPATHDB_KEY=\u{201c}SENTINEL-A\u{201d}\n");
    let before = s.requests().len();
    let err = fetch_fungidb("Testus fungus", &key, &s.base, &dir.path().join("o")).unwrap_err();
    assert_no_leak(&err, &["SENTINEL-A"]);
    assert!(err.to_string().contains("not allowed"), "{err}");
    assert_eq!(s.requests().len(), before, "no request may be sent");
}

#[test]
fn key_with_non_ascii_is_rejected_without_leak() {
    let s = list_server(|| Response::ok(b"x".as_slice()));
    let dir = tempfile::tempdir().unwrap();
    let key = write_key_file(dir.path(), "SENTINEL-\u{e9}-C\n");
    let err = fetch_fungidb("Testus fungus", &key, &s.base, &dir.path().join("o")).unwrap_err();
    assert_no_leak(&err, &["SENTINEL", "-C"]);
    let text = err.to_string();
    assert!(text.contains("key.txt") && text.contains("not allowed"), "{text}");
    assert_eq!(s.requests().len(), 0, "no request may be sent");
}

#[test]
fn transport_error_on_authenticated_request_has_no_key() {
    let s = list_server(Response::close);
    let dir = tempfile::tempdir().unwrap();
    let key = write_key_file(dir.path(), "SENTINEL-TRANSPORT-KEY\n");
    let err = fetch_fungidb("Testus fungus", &key, &s.base, &dir.path().join("o")).unwrap_err();
    assert!(matches!(err, FetchError::Http { status: None, .. }), "{err:?}");
    assert_no_leak(&err, &["SENTINEL-TRANSPORT-KEY"]);
}
