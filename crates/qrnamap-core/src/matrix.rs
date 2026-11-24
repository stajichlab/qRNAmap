//! The gene-by-sample count matrix step shared by the CLI and the browser page.
//! No file or network I/O.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;

use crate::counts::{assemble, unlisted_genes, CountMatrix};

/// At most this many unlisted gene names are shown in a warning.
const MAX_LISTED: usize = 20;

/// Errors from the matrix step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixError {
    /// A gene label has a tab, newline or carriage return, which would break the TSV.
    BadGeneLabel(String),
}

impl fmt::Display for MatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixError::BadGeneLabel(l) => write!(
                f,
                "gene label {l:?} has a tab, newline or carriage return and cannot go in a table"
            ),
        }
    }
}

impl std::error::Error for MatrixError {}

/// A `tx2gene` table with a row that cannot be read. `line` is 1-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tx2GeneError {
    pub line: usize,
    pub detail: String,
}

impl fmt::Display for Tx2GeneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "tx2gene line {}: {}", self.line, self.detail)
    }
}

impl std::error::Error for Tx2GeneError {}

fn check_label(l: &str) -> Result<(), MatrixError> {
    if l.contains(['\t', '\n', '\r']) {
        return Err(MatrixError::BadGeneLabel(l.to_string()));
    }
    Ok(())
}

/// Unique gene ids of a `transcript\tgene` table, in order of first appearance. Blank lines and
/// CRLF are accepted. A row without exactly two columns, or with an empty gene id, is an error.
pub fn genes_from_tx2gene(text: &str) -> Result<Vec<String>, Tx2GeneError> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for (i, l) in text.lines().enumerate() {
        let l = l.trim_end_matches('\r');
        if l.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        let [_, g] = f.as_slice() else {
            return Err(Tx2GeneError {
                line: i + 1,
                detail: "unexpected number of columns".into(),
            });
        };
        if g.is_empty() {
            return Err(Tx2GeneError {
                line: i + 1,
                detail: "empty gene id".into(),
            });
        }
        if seen.insert(g.to_string()) {
            out.push(g.to_string());
        }
    }
    Ok(out)
}

/// Read a `transcript\tgene` table (no header) into a map. Blank lines and CRLF are accepted. A row
/// without exactly two non-empty columns, or a transcript id that occurs twice, is an error.
pub fn parse_tx2gene(text: &str) -> Result<HashMap<String, String>, Tx2GeneError> {
    let mut m = HashMap::new();
    for (i, l) in text.lines().enumerate() {
        if l.trim().is_empty() {
            continue;
        }
        match l
            .trim_end_matches('\r')
            .split('\t')
            .collect::<Vec<_>>()
            .as_slice()
        {
            [t, g] if !t.is_empty() && !g.is_empty() => {
                if m.insert(t.to_string(), g.to_string()).is_some() {
                    return Err(Tx2GeneError {
                        line: i + 1,
                        detail: format!("transcript '{t}' occurs twice"),
                    });
                }
            }
            _ => {
                return Err(Tx2GeneError {
                    line: i + 1,
                    detail: "expected 'transcript<TAB>gene'".into(),
                })
            }
        }
    }
    Ok(m)
}

/// The warning for a sample whose counts are all zero or empty, else `None`.
pub fn zero_warning(sample: &str, counts: &BTreeMap<String, u64>) -> Option<String> {
    if counts.values().all(|&c| c == 0) {
        Some(format!(
            "sample '{sample}': all-zero sample (no counts in any gene)"
        ))
    } else {
        None
    }
}

/// Build the gene-by-sample matrix and the warnings for the user.
///
/// Rows: `genes` when given, else the sorted union of the sample genes. Any gene label with a
/// tab, newline or CR is `BadGeneLabel`. Warnings: genes in a sample but not in `genes`
/// (dropped from the matrix), and each sample whose counts are all zero or empty.
pub fn count_matrix(
    per_sample: &[(String, BTreeMap<String, u64>)],
    genes: Option<&[String]>,
) -> Result<(CountMatrix, Vec<String>), MatrixError> {
    for (_, m) in per_sample {
        m.keys().try_for_each(|g| check_label(g))?;
    }
    let mut warnings = Vec::new();
    let rows: Vec<String> = match genes {
        Some(g) => {
            g.iter().try_for_each(|x| check_label(x))?;
            let extra = unlisted_genes(per_sample, g);
            if !extra.is_empty() {
                let shown: Vec<&str> = extra.iter().take(MAX_LISTED).map(String::as_str).collect();
                let more = if extra.len() > MAX_LISTED {
                    ", ..."
                } else {
                    ""
                };
                warnings.push(format!(
                    "{} gene(s) in the sample tables are not in the gene list and are left out of the matrix: {}{more}",
                    extra.len(),
                    shown.join(", ")
                ));
            }
            g.to_vec()
        }
        None => per_sample
            .iter()
            .flat_map(|(_, m)| m.keys().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect(),
    };
    for (s, m) in per_sample {
        if let Some(w) = zero_warning(s, m) {
            warnings.push(w);
        }
    }
    Ok((assemble(per_sample, &rows), warnings))
}
