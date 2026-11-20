//! FASTA parsing over a byte reader. No file access and no threads.

use std::io::{self, BufRead};

/// Read all FASTA records. The name is the first word of the header. Lines of a record are joined.
/// Bases keep their case (differs from qmutfind, which uppercased them: qRNAmap keeps soft-masked bases as given).
pub fn read_fasta<R: BufRead>(reader: R) -> io::Result<Vec<(String, Vec<u8>)>> {
    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            contigs.push((name, Vec::new()));
        } else if let Some((_, seq)) = contigs.last_mut() {
            seq.extend(line.trim().bytes());
        }
    }
    Ok(contigs)
}
