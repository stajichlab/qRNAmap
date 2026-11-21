//! Count tables: kallisto abundance and featureCounts parsers, and the gene-by-sample matrix.
//! The functions take text, not paths. They do no file I/O.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;

/// Errors from the count parsers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CountError {
    /// A transcript id in the abundance table is not in `tx2gene`.
    UnknownTranscript(String),
    /// A row cannot be read. `line` is 1-based over the whole input, comment lines included.
    BadRow { line: usize },
    /// A required header column is missing from the abundance table.
    MissingColumn(String),
    /// A featureCounts table has no header line.
    MissingHeader,
    /// A gene id occurs in two rows of a featureCounts table.
    DuplicateGene(String),
}

impl fmt::Display for CountError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CountError::UnknownTranscript(id) => {
                write!(f, "transcript '{id}' is not in the tx2gene table")
            }
            CountError::BadRow { line } => write!(f, "bad count row at line {line}"),
            CountError::MissingColumn(c) => write!(f, "missing column '{c}' in the header"),
            CountError::MissingHeader => write!(f, "no header line found"),
            CountError::DuplicateGene(g) => write!(f, "gene '{g}' occurs in more than one row"),
        }
    }
}

impl std::error::Error for CountError {}

/// Gene counts from a kallisto `abundance.tsv`.
///
/// Columns are found by header name (`target_id`, `est_counts`). The count of a gene is the sum of
/// `est_counts` over its transcripts, rounded once per gene, half away from zero. A transcript id that
/// is not in `tx2gene` is `UnknownTranscript`. A negative or non-finite `est_counts`, or a row with too
/// few columns, is `BadRow`. Blank lines and CRLF are accepted.
///
/// A gene that is in `tx2gene` but has no row in the table gets no entry. `assemble` gives it 0.
/// A gene whose rows all have `est_counts` 0 is in the result with 0.
pub fn gene_counts_from_abundance(
    abundance_tsv: &str,
    tx2gene: &HashMap<String, String>,
) -> Result<BTreeMap<String, u64>, CountError> {
    sum_by_gene(abundance_tsv, tx2gene, "target_id", "est_counts")
}

/// Gene counts from a salmon `quant.sf`.
///
/// Same rules as `gene_counts_from_abundance`, with columns `Name` and `NumReads`.
pub fn gene_counts_from_salmon(
    quant_sf: &str,
    tx2gene: &HashMap<String, String>,
) -> Result<BTreeMap<String, u64>, CountError> {
    sum_by_gene(quant_sf, tx2gene, "Name", "NumReads")
}

/// Sum a count column over the transcripts of each gene and round once per gene.
fn sum_by_gene(
    abundance_tsv: &str,
    tx2gene: &HashMap<String, String>,
    id_name: &str,
    count_name: &str,
) -> Result<BTreeMap<String, u64>, CountError> {
    let mut lines = abundance_tsv
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l))
        .filter(|(_, l)| !l.trim().is_empty());
    let (_, header) = lines
        .next()
        .ok_or_else(|| CountError::MissingColumn(id_name.into()))?;
    let cols: Vec<&str> = header.split('\t').collect();
    let find = |name: &str| {
        cols.iter()
            .position(|c| *c == name)
            .ok_or_else(|| CountError::MissingColumn(name.to_string()))
    };
    let id_col = find(id_name)?;
    let est_col = find(count_name)?;

    let mut sums: BTreeMap<&str, f64> = BTreeMap::new();
    for (line, row) in lines {
        let f: Vec<&str> = row.split('\t').collect();
        let (Some(id), Some(est)) = (f.get(id_col), f.get(est_col)) else {
            return Err(CountError::BadRow { line });
        };
        let est: f64 = est
            .trim()
            .parse()
            .map_err(|_| CountError::BadRow { line })?;
        if !est.is_finite() || est < 0.0 {
            return Err(CountError::BadRow { line });
        }
        let gene = tx2gene
            .get(*id)
            .ok_or_else(|| CountError::UnknownTranscript(id.to_string()))?;
        *sums.entry(gene.as_str()).or_insert(0.0) += est;
    }
    Ok(sums
        .into_iter()
        .map(|(g, s)| (g.to_string(), s.round() as u64))
        .collect())
}

/// Gene counts from a featureCounts table.
///
/// Lines that start with `#` are skipped (the `#` line is optional). The next line is the header; its
/// first cell must be exactly `Geneid`, else `MissingHeader` (also for empty or comment-only input).
/// Every later non-blank row must have the same number of tab-separated fields as the header, else
/// `BadRow`. Each row gives the gene id (first column) and the count (last column, a non-negative integer). The Chr, Start,
/// End and Strand columns hold semicolon-separated lists with `-F SAF`. They are not parsed.
/// A bad count is `BadRow`; line numbers are 1-based over the whole input. A repeated gene id is
/// `DuplicateGene`. Blank lines and CRLF are accepted.
pub fn parse_featurecounts(tsv: &str) -> Result<BTreeMap<String, u64>, CountError> {
    // '#' lines are skipped only before the header; after it a line starting with '#' is a
    // gene row (R29-7).
    let mut lines = tsv
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l))
        .filter(|(_, l)| !l.trim().is_empty())
        .skip_while(|(_, l)| l.starts_with('#'));
    let (_, header) = lines.next().ok_or(CountError::MissingHeader)?;
    let header_fields: Vec<&str> = header.split('\t').collect();
    if header_fields[0] != "Geneid" {
        return Err(CountError::MissingHeader);
    }
    let width = header_fields.len();
    let mut out = BTreeMap::new();
    for (line, row) in lines {
        let f: Vec<&str> = row.split('\t').collect();
        if f.len() != width || f.len() < 2 || f[0].is_empty() {
            return Err(CountError::BadRow { line });
        }
        let last = f[f.len() - 1].trim();
        if last.is_empty() || !last.bytes().all(|b| b.is_ascii_digit()) {
            return Err(CountError::BadRow { line });
        }
        let n: u64 = last.parse().map_err(|_| CountError::BadRow { line })?;
        if out.insert(f[0].to_string(), n).is_some() {
            return Err(CountError::DuplicateGene(f[0].to_string()));
        }
    }
    Ok(out)
}

/// Category counts from a featureCounts `.summary` file (one input file).
///
/// The first non-blank line must start with `Status\t`, else `MissingHeader`. Every later
/// non-blank line is `<category>\t<count>`; anything else is `BadRow` (1-based line). CRLF is
/// accepted. A repeated category is `DuplicateGene` (the field holds the category).
pub fn parse_featurecounts_summary(text: &str) -> Result<BTreeMap<String, u64>, CountError> {
    let mut lines = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.strip_suffix('\r').unwrap_or(l)))
        .filter(|(_, l)| !l.trim().is_empty());
    let (_, header) = lines.next().ok_or(CountError::MissingHeader)?;
    if !header.starts_with("Status\t") {
        return Err(CountError::MissingHeader);
    }
    let mut out = BTreeMap::new();
    for (line, row) in lines {
        let f: Vec<&str> = row.split('\t').collect();
        if f.len() != 2 || f[0].is_empty() || !f[1].bytes().all(|b| b.is_ascii_digit()) {
            return Err(CountError::BadRow { line });
        }
        let n: u64 = f[1].parse().map_err(|_| CountError::BadRow { line })?;
        if out.insert(f[0].to_string(), n).is_some() {
            return Err(CountError::DuplicateGene(f[0].to_string()));
        }
    }
    Ok(out)
}

/// An integer gene-by-sample matrix. `counts[gene][sample]`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CountMatrix {
    pub genes: Vec<String>,
    pub samples: Vec<String>,
    pub counts: Vec<Vec<u64>>,
}

/// Build the matrix. `genes` fixes the row set and order. Sample order follows `per_sample`.
/// A gene missing from a sample is 0. Genes in a sample map that are not in `genes` are dropped;
/// `unlisted_genes` reports them. Duplicate sample names are not checked; the caller must prevent them.
pub fn assemble(per_sample: &[(String, BTreeMap<String, u64>)], genes: &[String]) -> CountMatrix {
    let counts = genes
        .iter()
        .map(|g| {
            per_sample
                .iter()
                .map(|(_, m)| m.get(g).copied().unwrap_or(0))
                .collect()
        })
        .collect();
    CountMatrix {
        genes: genes.to_vec(),
        samples: per_sample.iter().map(|(s, _)| s.clone()).collect(),
        counts,
    }
}

/// Genes that occur in any sample map but not in `genes`, sorted and without duplicates.
/// `assemble` drops these genes.
pub fn unlisted_genes(
    per_sample: &[(String, BTreeMap<String, u64>)],
    genes: &[String],
) -> Vec<String> {
    let listed: HashSet<&str> = genes.iter().map(String::as_str).collect();
    let extra: BTreeSet<&String> = per_sample
        .iter()
        .flat_map(|(_, m)| m.keys())
        .filter(|g| !listed.contains(g.as_str()))
        .collect();
    extra.into_iter().cloned().collect()
}

impl CountMatrix {
    /// Tab-separated text. Header `gene\t<sample>...`, then one row per gene. Ends with a newline.
    pub fn to_tsv(&self) -> String {
        let mut s = String::from("gene");
        for name in &self.samples {
            s.push('\t');
            s.push_str(name);
        }
        s.push('\n');
        for (g, row) in self.genes.iter().zip(&self.counts) {
            s.push_str(g);
            for c in row {
                s.push('\t');
                s.push_str(&c.to_string());
            }
            s.push('\n');
        }
        s
    }
}
