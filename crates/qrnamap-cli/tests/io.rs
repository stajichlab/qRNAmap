use flate2::write::GzEncoder;
use flate2::Compression;
use qrnamap_cli::io::{open_input, open_reader, InputError};
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
    // A 0-byte file without a compressed name is an empty plain file. (A 0-byte `.gz` or
    // `.zst` file is the named error `Empty`, final review R29-8.)
    let (_d, p) = write_tmp("e.txt", b"");
    assert_eq!(read_all(&p), "");
}

#[test]
fn empty_compressed_name_is_named_error() {
    // R29-8: a 0-byte file named .gz or .zst is not valid compressed data.
    for name in ["e.gz", "e.fq.zst"] {
        let (_d, p) = write_tmp(name, b"");
        match open_input(&p) {
            Err(InputError::Empty { path }) => assert_eq!(path, p),
            Err(other) => panic!("{name}: expected Empty, got {other:?}"),
            Ok(_) => panic!("{name}: expected Empty, got a reader"),
        }
        let e = open_input(&p).err().unwrap().to_string();
        assert!(e.contains("empty") && e.contains(name), "{e}");
    }
}

/// Read to the first error, then read again: the second read must give the same named error,
/// not a clean end of file (R29-8).
fn assert_error_is_sticky(path: &Path, want: fn(&InputError) -> bool) {
    let mut r = open_input(path).unwrap();
    let mut out = Vec::new();
    let first = r.read_to_end(&mut out).unwrap_err();
    let inner = first
        .get_ref()
        .and_then(|e| e.downcast_ref::<InputError>())
        .expect("inner InputError");
    assert!(want(inner), "first error: {inner:?}");
    let mut buf = [0u8; 64];
    let second = r.read(&mut buf);
    match second {
        Err(e) => {
            let inner = e
                .get_ref()
                .and_then(|e| e.downcast_ref::<InputError>())
                .expect("inner InputError");
            assert!(want(inner), "second error: {inner:?}");
        }
        Ok(n) => panic!(
            "a read after the error returned Ok({n}) for {}",
            path.display()
        ),
    }
}

#[test]
fn gzip_error_is_returned_again_on_later_reads() {
    let data: Vec<u8> = (0..20000u32)
        .flat_map(|i| i.to_string().into_bytes())
        .collect();
    let full = gz(&data);
    let (_d, p) = write_tmp("t.gz", &full[..full.len() / 2]);
    assert_error_is_sticky(&p, |e| matches!(e, InputError::Truncated { .. }));
    let mut bad = full.clone();
    let n = bad.len();
    bad[n - 8] ^= 0xff; // CRC32 field
    let (_d2, p2) = write_tmp("c.gz", &bad);
    assert_error_is_sticky(&p2, |e| matches!(e, InputError::Corrupt { .. }));
}

#[test]
fn zstd_error_is_returned_again_on_later_reads() {
    let data: Vec<u8> = (0..20000u32)
        .flat_map(|i| i.to_string().into_bytes())
        .collect();
    let full = zstd::encode_all(&data[..], 3).unwrap();
    let (_d, p) = write_tmp("t.zst", &full[..full.len() / 2]);
    assert_error_is_sticky(&p, |e| matches!(e, InputError::Truncated { .. }));
}

#[test]
fn missing_file_is_io_error() {
    let r = open_input(Path::new("/nonexistent/qrnamap/none.txt"));
    assert!(matches!(r, Err(InputError::Io { .. })));
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
    let data: Vec<u8> = (0..20000u32)
        .flat_map(|i| i.to_string().into_bytes())
        .collect();
    let full = gz(&data);
    let (_d, p) = write_tmp("t.gz", &full[..full.len() / 2]);
    assert_truncated(&p);
}

#[test]
fn truncated_zstd_is_named_error() {
    let data: Vec<u8> = (0..20000u32)
        .flat_map(|i| i.to_string().into_bytes())
        .collect();
    let full = zstd::encode_all(&data[..], 3).unwrap();
    let (_d, p) = write_tmp("t.zst", &full[..full.len() / 2]);
    assert_truncated(&p);
}

fn read_until_end_or_error(path: &Path) -> (Vec<u8>, std::io::Error) {
    let mut r = open_input(path).unwrap();
    let mut out = Vec::new();
    let err = r.read_to_end(&mut out).unwrap_err();
    (out, err)
}

fn assert_corrupt(err: &std::io::Error, path: &Path) -> String {
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    let inner = err
        .get_ref()
        .and_then(|e| e.downcast_ref::<InputError>())
        .expect("inner InputError");
    match inner {
        InputError::Corrupt { path: p, detail } => {
            assert_eq!(p, path);
            let text = inner.to_string();
            assert!(text.contains("corrupt"));
            assert!(text.contains(&path.display().to_string()));
            detail.clone()
        }
        other => panic!("wrong variant: {other:?}"),
    }
}

#[test]
fn corrupt_gzip_checksum_is_named_corrupt_error() {
    let mut b = gz(b"some text to corrupt\n");
    let n = b.len();
    b[n - 5] ^= 0xff; // CRC32 field
    let (_d, p) = write_tmp("c.gz", &b);
    let (_out, err) = read_until_end_or_error(&p);
    let detail = assert_corrupt(&err, &p);
    assert!(detail.contains("checksum"), "detail: {detail}");
}

#[test]
fn corrupt_zstd_is_named_corrupt_error() {
    let mut enc = zstd::stream::Encoder::new(Vec::new(), 3).unwrap();
    enc.include_checksum(true).unwrap();
    enc.write_all(b"zstd text to corrupt\n").unwrap();
    let mut b = enc.finish().unwrap();
    let n = b.len();
    b[n - 1] ^= 0xff; // checksum field
    let (_d, p) = write_tmp("c.zst", &b);
    let (_out, err) = read_until_end_or_error(&p);
    let detail = assert_corrupt(&err, &p);
    assert!(detail.contains("checksum"), "detail: {detail}");
}

#[test]
fn gzip_trailing_junk_is_corrupt_not_truncated() {
    let mut b = gz(b"first\n");
    b.extend(gz(b"second\n"));
    b.extend(b"junk\n");
    let (_d, p) = write_tmp("j.gz", &b);
    let (out, err) = read_until_end_or_error(&p);
    assert_eq!(out, b"first\nsecond\n");
    let detail = assert_corrupt(&err, &p);
    assert_eq!(detail, "trailing bytes after last member");
}

#[test]
fn gzip_trailing_zero_padding_is_corrupt() {
    let mut b = gz(b"data\n");
    b.extend(vec![0u8; 16]);
    let (_d, p) = write_tmp("z.gz", &b);
    let (out, err) = read_until_end_or_error(&p);
    assert_eq!(out, b"data\n");
    let detail = assert_corrupt(&err, &p);
    assert_eq!(detail, "trailing bytes after last member");
}

struct OneByte(std::io::Cursor<Vec<u8>>);
impl Read for OneByte {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = buf.len().min(1);
        self.0.read(&mut buf[..n])
    }
}

#[test]
fn short_first_read_still_detects_gzip() {
    let b = gz(b"slow pipe\nline2\n");
    let mut r = open_reader(OneByte(std::io::Cursor::new(b)), Path::new("pipe.gz")).unwrap();
    let mut s = String::new();
    r.read_to_string(&mut s).unwrap();
    assert_eq!(s, "slow pipe\nline2\n");
}

#[test]
fn short_first_read_still_detects_zstd() {
    let b = zstd::encode_all(&b"slow zstd\n"[..], 3).unwrap();
    let mut r = open_reader(OneByte(std::io::Cursor::new(b)), Path::new("pipe")).unwrap();
    let mut s = String::new();
    r.read_to_string(&mut s).unwrap();
    assert_eq!(s, "slow zstd\n");
}

#[test]
fn open_missing_file_error_names_path() {
    let path = Path::new("/nonexistent/qrnamap/none.txt");
    let err = match open_input(path) {
        Err(e) => e,
        Ok(_) => panic!("expected error"),
    };
    assert!(err.to_string().contains("/nonexistent/qrnamap/none.txt"));
}
