//! FASTQ parsing over a byte reader, and a chunked parser for browser input. No file access and no threads.

use std::io::{self, BufRead};

/// Iterator over FASTQ records as (name, sequence). The name is the first word of the header.
pub struct FastqReader<R: BufRead> {
    lines: io::Lines<R>,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader { lines: reader.lines() }
    }

    fn next_line(&mut self, what: &str) -> io::Result<String> {
        match self.lines.next() {
            Some(line) => line,
            None => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("truncated FASTQ record: missing {what}"),
            )),
        }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = io::Result<(String, Vec<u8>)>;

    fn next(&mut self) -> Option<Self::Item> {
        let header = match self.lines.next()? {
            Ok(h) => h,
            Err(e) => return Some(Err(e)),
        };
        if header.trim().is_empty() {
            // Blank lines after the last record end the input (R29-12). A blank line followed by
            // more text is an error.
            for line in self.lines.by_ref() {
                match line {
                    Ok(l) if l.trim().is_empty() => {}
                    Ok(_) => {
                        return Some(Err(io::Error::new(io::ErrorKind::InvalidData, "blank line between FASTQ records")))
                    }
                    Err(e) => return Some(Err(e)),
                }
            }
            return None;
        }
        let record = (|| {
            let name = header
                .strip_prefix('@')
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "FASTQ header does not start with @"))?
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let seq = self.next_line("sequence")?;
            self.next_line("+ line")?;
            self.next_line("quality")?;
            Ok((name, seq.trim().bytes().map(|b| b.to_ascii_uppercase()).collect()))
        })();
        Some(record)
    }
}

/// Parse FASTQ from byte chunks that can split a record anywhere. Call `push` for each
/// chunk, then `finish` once. A truncated last record is an error.
pub struct FastqChunker {
    buffer: Vec<u8>,
    /// The most bytes of an unfinished record that are kept.
    limit: usize,
    /// Bytes at the start of `buffer` that were already searched for line breaks, and the line breaks found in them. A push searches
    /// only the new bytes, so a long record that arrives in many small chunks is not read again at every chunk.
    scanned: usize,
    newlines: usize,
    /// Bytes examined by the newline scan since the chunker was made. Each byte should be
    /// examined once; the tests check this instead of wall-clock time (R29-15).
    scan_work: u64,
}

/// The default limit for one unfinished record: 64 MiB. A Nanopore read is some kilobases to a few megabases.
pub const DEFAULT_MAX_RECORD_BYTES: usize = 64 << 20;

impl Default for FastqChunker {
    fn default() -> Self {
        FastqChunker::with_limit(DEFAULT_MAX_RECORD_BYTES)
    }
}

impl FastqChunker {
    pub fn new() -> Self {
        FastqChunker::default()
    }

    /// A chunker that gives an error when an unfinished record grows past `limit` bytes. Input with no line breaks, or text that is not
    /// FASTQ, would otherwise fill the memory of the tab before the end of the file shows the problem.
    pub fn with_limit(limit: usize) -> Self {
        FastqChunker { buffer: Vec::new(), limit, scanned: 0, newlines: 0, scan_work: 0 }
    }

    /// Add a chunk. Returns the records that are complete now.
    pub fn push(&mut self, chunk: &[u8]) -> io::Result<Vec<(String, Vec<u8>)>> {
        self.buffer.extend_from_slice(chunk);
        // A record ends at every fourth newline. Only the bytes added since the last push are searched.
        let mut consumed = 0;
        for (i, &b) in self.buffer.iter().enumerate().skip(self.scanned) {
            self.scan_work += 1;
            if b == b'\n' {
                self.newlines += 1;
                if self.newlines.is_multiple_of(4) {
                    consumed = i + 1;
                }
            }
        }
        self.scanned = self.buffer.len();
        let records = FastqReader::new(&self.buffer[..consumed]).collect::<io::Result<Vec<_>>>()?;
        self.buffer.drain(..consumed);
        // The complete records are gone; the rest holds the line breaks after the last record end (0 to 3 of them).
        self.scanned = self.buffer.len();
        self.newlines %= 4;
        if self.buffer.len() > self.limit {
            let limit = self.limit;
            self.buffer = Vec::new();
            self.scanned = 0;
            self.newlines = 0;
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("a FASTQ record is larger than {limit} bytes: the input may not be FASTQ, or has no line breaks"),
            ));
        }
        Ok(records)
    }

    /// Finish the input. Returns a last record that had no trailing newline, or an error when
    /// the input stops inside a record.
    pub fn finish(&mut self) -> io::Result<Vec<(String, Vec<u8>)>> {
        if self.buffer.iter().all(|b| b.is_ascii_whitespace()) {
            self.buffer.clear();
            self.scanned = 0;
            self.newlines = 0;
            return Ok(Vec::new());
        }
        self.buffer.push(b'\n');
        let records = FastqReader::new(&self.buffer[..]).collect::<io::Result<Vec<_>>>()?;
        self.buffer.clear();
        self.scanned = 0;
        self.newlines = 0;
        Ok(records)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_long_record_arriving_in_small_chunks_is_scanned_once() {
        // 6 MB in one record, in chunks of 2 KiB: 3,000 pushes. A scan from the start of the
        // buffer on every push would examine about 9 GB; a scan of the new bytes examines each
        // byte once.
        let record = format!("@long\n{}\n+\n{}\n", "A".repeat(3_000_000), "I".repeat(3_000_000));
        let mut chunker = FastqChunker::new();
        let mut n = 0;
        for chunk in record.as_bytes().chunks(2048) {
            n += chunker.push(chunk).unwrap().len();
        }
        n += chunker.finish().unwrap().len();
        assert_eq!(n, 1);
        assert_eq!(chunker.scan_work, record.len() as u64, "each byte is scanned once");
    }

    #[test]
    fn a_record_with_no_line_breaks_reaches_the_limit_after_one_scan_of_each_byte() {
        let limit = 4 << 20;
        let mut chunker = FastqChunker::with_limit(limit);
        let piece = vec![b'A'; 2048];
        let mut pushed = 0u64;
        let mut failed = false;
        for _ in 0..4000 {
            pushed += piece.len() as u64;
            if chunker.push(&piece).is_err() {
                failed = true;
                break;
            }
        }
        assert!(failed, "the limit was not reached");
        assert!(pushed <= limit as u64 + 2048, "the error came late: {pushed} bytes pushed");
        assert_eq!(chunker.scan_work, pushed, "each byte is scanned once");
    }
}
