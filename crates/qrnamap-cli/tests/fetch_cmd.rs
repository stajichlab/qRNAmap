//! `qrnamap fetch` through the binary, against a local HTTP server (no network).
//!
//! NCBI tests pass a catalog file with the hidden test-only option `--catalog FILE`, so the
//! genome URLs point at the local server. No API key ever reaches the output.

mod common;

use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use common::fetch_server::{Response, Server};

const FASTA: &[u8] = b"fasta-gz-bytes-0123456789";
const GFF: &[u8] = b"gff-gz-bytes";

fn qrnamap(args: &[&str], env: &[(&str, &str)]) -> Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_qrnamap"));
    c.args(args).env_remove("EUPATHDB_KEY");
    for (k, v) in env {
        c.env(k, v);
    }
    c.output().unwrap()
}

fn text(o: &Output) -> (String, String) {
    (
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// A catalog JSON with one genome whose files are on `base`.
fn catalog(dir: &Path, base: &str, with_gff: bool, genome_size: u64) -> std::path::PathBuf {
    let gff = if with_gff {
        format!(
            r#","gff":{{"name":"t_genomic.gff.gz","url":"{base}/t_genomic.gff.gz","bytes":{}}}"#,
            GFF.len()
        )
    } else {
        String::new()
    };
    let json = format!(
        r#"{{"schema_version":1,"generated":"2026-11-19","genomes":[{{"accession":"GCF_000000001.1","source":"RefSeq","label":"test","organism":"Testus fungus","level":"Complete","assembly":"ASM1","genome_size":{genome_size},"fasta":{{"name":"t_genomic.fna.gz","url":"{base}/t_genomic.fna.gz","bytes":{}}}{gff}}}]}}"#,
        FASTA.len()
    );
    let p = dir.join("catalog.json");
    fs::write(&p, json).unwrap();
    p
}

fn ncbi_server() -> Server {
    Server::start(|req, _| match req.path.as_str() {
        "/t_genomic.fna.gz" => Response::ok(FASTA),
        "/t_genomic.gff.gz" => Response::ok(GFF),
        _ => Response::status(404),
    })
}

#[test]
fn fetch_ncbi_writes_files_to_out_dir() {
    let srv = ncbi_server();
    let dir = tempfile::tempdir().unwrap();
    let cat = catalog(dir.path(), &srv.base, true, 1000);
    let out = dir.path().join("genome");
    let o = qrnamap(
        &[
            "fetch",
            "--source",
            "ncbi",
            "--id",
            "GCF_000000001.1",
            "--out",
            s(&out),
            "--catalog",
            s(&cat),
        ],
        &[],
    );
    let (stdout, stderr) = text(&o);
    assert!(o.status.success(), "{stdout}\n{stderr}");
    assert_eq!(fs::read(out.join("t_genomic.fna.gz")).unwrap(), FASTA);
    assert_eq!(fs::read(out.join("t_genomic.gff.gz")).unwrap(), GFF);
    assert_eq!(fs::read_dir(&out).unwrap().count(), 2, "no other files");
    assert!(
        stdout.contains("t_genomic.fna.gz") && stdout.contains("t_genomic.gff.gz"),
        "{stdout}"
    );
    assert!(!stderr.contains("warning"), "{stderr}");
}

#[test]
fn fetch_catalog_option_is_hidden_from_help() {
    let o = qrnamap(&["fetch", "--help"], &[]);
    let (stdout, _) = text(&o);
    assert!(o.status.success());
    assert!(!stdout.contains("--catalog"), "{stdout}");
    assert!(stdout.contains("--base-url"), "{stdout}");
}

#[test]
fn fetch_unknown_accession_is_named_error() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("g");
    let o = qrnamap(
        &[
            "fetch",
            "--source",
            "ncbi",
            "--id",
            "GCF_999999999.9",
            "--out",
            s(&out),
        ],
        &[],
    );
    let (_, stderr) = text(&o);
    assert!(!o.status.success());
    assert!(
        stderr.contains("GCF_999999999.9") && stderr.contains("not in the genome catalog"),
        "{stderr}"
    );
    assert!(!out.exists(), "no output directory for a failed lookup");
}

#[test]
fn fetch_no_gff_genome_warns() {
    let srv = ncbi_server();
    let dir = tempfile::tempdir().unwrap();
    let cat = catalog(dir.path(), &srv.base, false, 60_000_000);
    let out = dir.path().join("genome");
    let o = qrnamap(
        &[
            "fetch",
            "--source",
            "ncbi",
            "--id",
            "GCF_000000001.1",
            "--out",
            s(&out),
            "--catalog",
            s(&cat),
        ],
        &[],
    );
    let (stdout, stderr) = text(&o);
    assert!(o.status.success(), "{stdout}\n{stderr}");
    assert_eq!(fs::read(out.join("t_genomic.fna.gz")).unwrap(), FASTA);
    assert!(!out.join("t_genomic.gff.gz").exists());
    assert!(
        stderr.contains("warning") && stderr.contains("cannot be used for RNA-seq counting"),
        "{stderr}"
    );
    // The catalog's size warning (genome above 50 Mb) is printed too.
    assert!(stderr.contains("60.0 Mb"), "{stderr}");
}

/// FungiDB-like server: the organism list is public; files need `Bearer <want>`.
fn fungidb_server(want: &'static str) -> Server {
    Server::start(move |req, base| {
        if req.path.contains("/searches/GenomeDataTypes/") {
            let body = format!(
                r#"{{"records":[{{"attributes":{{"primary_key":"Testus fungus","URLGenomeFasta":"{base}/a/t.fasta","URLgff":"{base}/a/t.gff"}}}}]}}"#
            );
            return Response::ok(body.into_bytes());
        }
        if req.header("authorization") != Some(&format!("Bearer {want}")) {
            return Response::status(401);
        }
        match req.path.as_str() {
            "/a/t.fasta" => Response::ok(&b">c\nACGT\n"[..]),
            "/a/t.gff" => Response::ok(&b"##gff-version 3\n"[..]),
            _ => Response::status(404),
        }
    })
}

fn key_file(dir: &Path, text: &str) -> std::path::PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let p = dir.join("key.txt");
    fs::write(&p, text).unwrap();
    fs::set_permissions(&p, fs::Permissions::from_mode(0o600)).unwrap();
    p
}

#[test]
fn fetch_key_file_wins_over_env() {
    let srv = fungidb_server("FILEKEY-SENTINEL");
    let dir = tempfile::tempdir().unwrap();
    let kf = key_file(dir.path(), "FILEKEY-SENTINEL\n");
    let out = dir.path().join("fdb");
    let o = qrnamap(
        &[
            "fetch",
            "--source",
            "fungidb",
            "--id",
            "Testus fungus",
            "--out",
            s(&out),
            "--eupathdb-key-file",
            s(&kf),
            "--base-url",
            &srv.base,
        ],
        &[("EUPATHDB_KEY", "ENVKEY-SENTINEL")],
    );
    let (stdout, stderr) = text(&o);
    assert!(o.status.success(), "{stdout}\n{stderr}");
    assert!(out.join("t.fasta.gz").is_file() && out.join("t.gff.gz").is_file());
    let auth: Vec<String> = srv
        .requests()
        .iter()
        .filter_map(|r| r.header("authorization").map(str::to_string))
        .collect();
    assert_eq!(auth.len(), 2, "{auth:?}");
    assert!(
        auth.iter().all(|a| a == "Bearer FILEKEY-SENTINEL"),
        "{auth:?}"
    );
    for k in ["FILEKEY", "ENVKEY"] {
        assert!(
            !stdout.contains(k) && !stderr.contains(k),
            "key leaked: {stdout}\n{stderr}"
        );
    }
}

#[test]
fn fetch_env_key_is_used_without_key_file() {
    let srv = fungidb_server("ENVKEY-SENTINEL");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("fdb");
    let o = qrnamap(
        &[
            "fetch",
            "--source",
            "fungidb",
            "--id",
            "Testus fungus",
            "--out",
            s(&out),
            "--base-url",
            &srv.base,
        ],
        &[("EUPATHDB_KEY", "ENVKEY-SENTINEL")],
    );
    let (stdout, stderr) = text(&o);
    assert!(o.status.success(), "{stdout}\n{stderr}");
    assert!(!stdout.contains("ENVKEY") && !stderr.contains("ENVKEY"));
}

#[test]
fn fetch_empty_key_is_named_error() {
    let srv = fungidb_server("X");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("fdb");
    let base = srv.base.clone();
    let run = |extra: &[&str], env: &[(&str, &str)]| {
        let mut a = vec![
            "fetch",
            "--source",
            "fungidb",
            "--id",
            "Testus fungus",
            "--out",
            s(&out),
            "--base-url",
            &base,
        ];
        a.extend_from_slice(extra);
        qrnamap(&a, env)
    };
    // Empty environment variable, no key file.
    let o = run(&[], &[("EUPATHDB_KEY", "  ")]);
    let (_, e1) = text(&o);
    assert!(!o.status.success());
    assert!(e1.contains("empty") && e1.contains("EUPATHDB_KEY"), "{e1}");
    // Key file with an empty assignment, and an empty key file.
    for (i, body) in ["EUPATHDB_KEY=\"\"\n", "# no key here\n\n"]
        .iter()
        .enumerate()
    {
        let sub = dir.path().join(format!("k{i}"));
        fs::create_dir_all(&sub).unwrap();
        let kf = key_file(&sub, body);
        let o = run(&["--eupathdb-key-file", s(&kf)], &[]);
        let (_, e) = text(&o);
        assert!(!o.status.success(), "{body:?}");
        assert!(
            e.contains("empty") && e.contains("key file"),
            "{body:?}: {e}"
        );
    }
    assert!(
        srv.requests().is_empty(),
        "no request may be made with an empty key"
    );
    assert!(!out.exists());
}

#[test]
fn fetch_base_url_is_fungidb_only() {
    let dir = tempfile::tempdir().unwrap();
    let o = qrnamap(
        &[
            "fetch",
            "--source",
            "ncbi",
            "--id",
            "X",
            "--out",
            s(dir.path()),
            "--base-url",
            "http://127.0.0.1:9",
        ],
        &[],
    );
    let (_, e) = text(&o);
    assert!(!o.status.success());
    assert!(e.contains("--base-url"), "{e}");
}
