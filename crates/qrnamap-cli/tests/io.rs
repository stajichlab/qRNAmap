use flate2::write::GzEncoder;
use flate2::Compression;
use qrnamap_cli::io::{open_input, InputError};
use std::io::{Read, Write};
use std::path::Path;

fn gz(data: &[u8]) -> Vec<u8> {
    let mut e = GzEncoder::new(Vec::new(), Compression::default());
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

fn read_all(path: &Path) -> String {
    let mut s = String::new();
    open_input(path).unwrap().read_to_string(&mut s).unwrap();
    s
}

fn write_tmp(name: &str, bytes: &[u8]) -> (tempfile::TempDir, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join(name);
    std::fs::write(&p, bytes).unwrap();
    (dir, p)
}

#[test]
fn reads_plain() {
    let (_d, p) = write_tmp("a.txt", b"hello\nworld\n");
    assert_eq!(read_all(&p), "hello\nworld\n");
}

#[test]
fn reads_gzip_single_member() {
    let (_d, p) = write_tmp("a.txt.gz", &gz(b"one\ntwo\n"));
    assert_eq!(read_all(&p), "one\ntwo\n");
}

#[test]
fn reads_gzip_multi_member() {
    let mut b = gz(b"first\n");
    b.extend(gz(b"second\n"));
    let (_d, p) = write_tmp("m.gz", &b);
    assert_eq!(read_all(&p), "first\nsecond\n");
}

#[test]
fn reads_zstd() {
    let b = zstd::encode_all(&b"zstd data\nline2\n"[..], 3).unwrap();
    let (_d, p) = write_tmp("a.txt.zst", &b);
    assert_eq!(read_all(&p), "zstd data\nline2\n");
}

#[test]
fn extension_wrong_magic_wins() {
    let (_d, p) = write_tmp("a.txt", &gz(b"really gzip\n"));
    assert_eq!(read_all(&p), "really gzip\n");
    let b = zstd::encode_all(&b"really zstd\n"[..], 3).unwrap();
    let (_d2, p2) = write_tmp("b.gz", &b);
    assert_eq!(read_all(&p2), "really zstd\n");
}

#[test]
fn empty_file_reads_empty() {
    let (_d, p) = write_tmp("e.gz", b"");
    assert_eq!(read_all(&p), "");
}

#[test]
fn missing_file_is_io_error() {
    let r = open_input(Path::new("/nonexistent/qrnamap/none.txt"));
    assert!(matches!(r, Err(InputError::Io(_))));
}

fn assert_truncated(path: &Path) {
    let mut r = open_input(path).unwrap();
    let mut out = Vec::new();
    let err = r.read_to_end(&mut out).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    let inner = err
        .get_ref()
        .and_then(|e| e.downcast_ref::<InputError>())
        .expect("inner InputError");
    match inner {
        InputError::Truncated { path: p } => assert_eq!(p, path),
        other => panic!("wrong variant: {other:?}"),
    }
    let text = inner.to_string();
    assert!(text.contains("truncated"));
    assert!(text.contains(&path.display().to_string()));
}

#[test]
fn truncated_gzip_is_named_error() {
    let data: Vec<u8> = (0..20000u32).flat_map(|i| i.to_string().into_bytes()).collect();
    let full = gz(&data);
    let (_d, p) = write_tmp("t.gz", &full[..full.len() / 2]);
    assert_truncated(&p);
}

#[test]
fn truncated_zstd_is_named_error() {
    let data: Vec<u8> = (0..20000u32).flat_map(|i| i.to_string().into_bytes()).collect();
    let full = zstd::encode_all(&data[..], 3).unwrap();
    let (_d, p) = write_tmp("t.zst", &full[..full.len() / 2]);
    assert_truncated(&p);
}
