//! `qrnamap run` and the reads-per-sample check.

mod common;

use std::fs;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use qrnamap_cli::io::open_reader;
use qrnamap_cli::run::{check_reads_per_sample, count_pairs, count_records};
use qrnamap_core::SampleSheet;

fn fastq(n: usize) -> Vec<u8> {
    let mut v = Vec::new();
    for i in 0..n {
        writeln!(v, "@r{i} extra\nACGTACGTAC\n+\nIIIIIIIIII").unwrap();
    }
    v
}

fn gz(data: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(data).unwrap();
    e.finish().unwrap()
}

#[test]
fn count_pairs_plain_and_gzip() {
    let d = tempfile::tempdir().unwrap();
    let plain = d.path().join("a.fq");
    fs::write(&plain, fastq(123)).unwrap();
    // Multi-member gzip: two members of 100 and 23 records.
    let multi = d.path().join("a.fq.gz");
    let mut bytes = gz(&fastq(100));
    bytes.extend(gz(&fastq(23)));
    fs::write(&multi, bytes).unwrap();
    let zst = d.path().join("a.fq.zst");
    fs::write(&zst, zstd::encode_all(&fastq(123)[..], 3).unwrap()).unwrap();
    for p in [&plain, &multi, &zst] {
        assert_eq!(count_pairs(p, None).unwrap(), 123, "{}", p.display());
        assert_eq!(count_pairs(p, Some(1000)).unwrap(), 123);
        assert_eq!(count_pairs(p, Some(123)).unwrap(), 123, "equal to the limit is not over it");
        assert_eq!(count_pairs(p, Some(50)).unwrap(), 51, "stops at limit + 1");
    }
    // A truncated record is an error, not a short count.
    let bad = d.path().join("bad.fq");
    fs::write(&bad, b"@r1\nACGT\n+\n").unwrap();
    assert!(count_pairs(&bad, None).is_err());
    // A truncated gzip stream is a named error.
    let full = gz(&fastq(500));
    let cut = d.path().join("cut.fq.gz");
    fs::write(&cut, &full[..full.len() / 2]).unwrap();
    let e = count_pairs(&cut, None).unwrap_err();
    assert!(e.to_string().contains("truncated"), "{e}");
}

/// A reader that records how many bytes were taken from it.
struct Counting {
    inner: std::io::Cursor<Vec<u8>>,
    taken: Arc<AtomicUsize>,
}

impl Read for Counting {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.taken.fetch_add(n, Ordering::SeqCst);
        Ok(n)
    }
}

#[test]
fn count_pairs_stops_after_limit() {
    let data = fastq(200_000);
    let total = data.len();
    let taken = Arc::new(AtomicUsize::new(0));
    let r = open_reader(Counting { inner: std::io::Cursor::new(data), taken: taken.clone() }, Path::new("mem.fq")).unwrap();
    assert_eq!(count_records(r, Path::new("mem.fq"), Some(10)).unwrap(), 11);
    let used = taken.load(Ordering::SeqCst);
    eprintln!("read {used} of {total} bytes");
    assert!(used < 64 * 1024, "read {used} of {total} bytes; it must stop after the limit");
}

fn sheet_with(d: &Path, rows: &[(&str, usize)]) -> SampleSheet {
    let mut t = String::from("sample\tfastq1\tfastq2\n");
    for (name, n) in rows {
        let r1 = d.join(format!("{name}_1.fq"));
        let r2 = d.join(format!("{name}_2.fq"));
        fs::write(&r1, fastq(*n)).unwrap();
        fs::write(&r2, fastq(*n)).unwrap();
        t.push_str(&format!("{name}\t{}\t{}\n", r1.display(), r2.display()));
    }
    SampleSheet::parse(&t).unwrap()
}

#[test]
fn exceeds_reads_per_sample_gives_warning_not_error() {
    let d = tempfile::tempdir().unwrap();
    let sheet = sheet_with(d.path(), &[("A", 30), ("B", 5), ("C", 6)]);
    let w = check_reads_per_sample(&sheet, 5).unwrap();
    let names: Vec<&str> = w.iter().map(|x| x.sample.as_str()).collect();
    assert_eq!(names, ["A", "C"]);
    let text = w[0].to_string();
    assert!(text.contains("'A'") && text.contains("more than 5") && text.contains("reads_per_sample"), "{text}");
    assert!(check_reads_per_sample(&sheet, 30).unwrap().is_empty());
}

#[test]
fn reads_per_sample_sums_lanes() {
    let d = tempfile::tempdir().unwrap();
    let (a, b) = (d.path().join("l1.fq"), d.path().join("l2.fq"));
    fs::write(&a, fastq(4)).unwrap();
    fs::write(&b, fastq(4)).unwrap();
    let t = format!("sample\tfastq1\nS\t{},{}\n", a.display(), b.display());
    let sheet = SampleSheet::parse(&t).unwrap();
    assert!(check_reads_per_sample(&sheet, 8).unwrap().is_empty());
    assert_eq!(check_reads_per_sample(&sheet, 7).unwrap().len(), 1);
}
