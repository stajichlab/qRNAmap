//! GFF3 gene model that keeps every isoform. No file or network I/O; reads from any `BufRead`.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, BufRead};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

/// One transcript. All intervals are 1-based inclusive and sorted ascending.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    /// GFF type of the transcript feature: mRNA, rRNA, tRNA, ncRNA, ...
    pub kind: String,
    pub exons: Vec<(u64, u64)>,
    pub cds: Vec<(u64, u64)>,
    /// Explicit five_prime_UTR features only.
    pub utr5: Vec<(u64, u64)>,
    /// Explicit three_prime_UTR features only.
    pub utr3: Vec<(u64, u64)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub label: String,
    pub contig: String,
    pub strand: Strand,
    pub transcripts: Vec<Transcript>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneModel {
    pub genes: Vec<Gene>,
}

#[derive(Debug)]
pub enum GffError {
    /// Fewer than 9 columns, or a bad start/end. `line` is 1-based.
    BadLine { line: usize },
    /// Strand is not `+` or `-` on a gene or transcript that has children.
    BadStrand { line: usize, value: String },
    DuplicateGeneLabel(String),
    Io(io::Error),
}

impl fmt::Display for GffError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GffError::BadLine { line } => write!(f, "GFF3 line {line}: expected 9 tab-separated columns with valid start and end"),
            GffError::BadStrand { line, value } => {
                write!(f, "GFF3 line {line}: strand must be '+' or '-', found '{value}'")
            }
            GffError::DuplicateGeneLabel(l) => write!(f, "two different genes have the label '{l}'"),
            GffError::Io(e) => write!(f, "I/O error reading GFF3: {e}"),
        }
    }
}

impl std::error::Error for GffError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            GffError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for GffError {
    fn from(e: io::Error) -> Self {
        GffError::Io(e)
    }
}

struct Feat {
    line: usize,
    kind: String,
    contig: String,
    strand: String,
    attrs: Vec<(String, String)>, // raw (still percent-encoded) values
}

impl Feat {
    fn raw(&self, key: &str) -> Option<&str> {
        self.attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str()).filter(|v| !v.is_empty())
    }
    fn get(&self, key: &str) -> Option<String> {
        self.raw(key).map(percent_decode)
    }
    fn parents(&self) -> Vec<String> {
        self.raw("Parent")
            .map(|p| p.split(',').filter(|s| !s.is_empty()).map(percent_decode).collect())
            .unwrap_or_default()
    }
    fn label(&self) -> Option<String> {
        self.get("locus_tag").or_else(|| self.get("Name")).or_else(|| self.get("ID"))
    }
}

#[derive(Default)]
struct Children {
    exons: Vec<(u64, u64)>,
    cds: Vec<(u64, u64)>,
    utr5: Vec<(u64, u64)>,
    utr3: Vec<(u64, u64)>,
}

fn percent_decode(s: &str) -> String {
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' && i + 2 < b.len() && b[i + 1].is_ascii_hexdigit() && b[i + 2].is_ascii_hexdigit() {
            let h = |c: u8| (c as char).to_digit(16).unwrap() as u8;
            out.push(h(b[i + 1]) * 16 + h(b[i + 2]));
            i += 3;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn finish(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    v.sort_unstable();
    v.dedup();
    v
}

fn strand_of(f: &Feat) -> Result<Strand, GffError> {
    match f.strand.as_str() {
        "+" => Ok(Strand::Plus),
        "-" => Ok(Strand::Minus),
        other => Err(GffError::BadStrand { line: f.line, value: other.to_string() }),
    }
}

fn is_child_kind(k: &str) -> bool {
    matches!(k, "exon" | "CDS" | "five_prime_UTR" | "three_prime_UTR")
}

pub fn parse_gff3<R: BufRead>(reader: R) -> Result<GeneModel, GffError> {
    // Pass 1: read features. Child features go straight into a map keyed by parent ID.
    let mut feats: Vec<Feat> = Vec::new();
    let mut children: HashMap<String, Children> = HashMap::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        let lineno = idx + 1;
        let line = line.trim_end_matches(['\r', '\n']);
        if line.starts_with("##FASTA") {
            break;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(GffError::BadLine { line: lineno });
        }
        let start: u64 = cols[3].trim().parse().map_err(|_| GffError::BadLine { line: lineno })?;
        let end: u64 = cols[4].trim().parse().map_err(|_| GffError::BadLine { line: lineno })?;
        if start == 0 || start > end {
            return Err(GffError::BadLine { line: lineno });
        }
        let attrs: Vec<(String, String)> = cols[8]
            .split(';')
            .filter_map(|kv| kv.trim().split_once('='))
            .map(|(k, v)| (k.trim().to_string(), v.trim().to_string()))
            .collect();
        let feat = Feat {
            line: lineno,
            kind: cols[2].to_string(),
            contig: cols[0].to_string(),
            strand: cols[6].trim().to_string(),
            attrs,
        };
        if is_child_kind(&feat.kind) {
            for p in feat.parents() {
                let c = children.entry(p).or_default();
                let list = match feat.kind.as_str() {
                    "exon" => &mut c.exons,
                    "CDS" => &mut c.cds,
                    "five_prime_UTR" => &mut c.utr5,
                    _ => &mut c.utr3,
                };
                list.push((start, end));
            }
        } else {
            feats.push(feat);
        }
    }

    // Pass 2: index parents by ID (first occurrence wins), then build genes.
    let mut by_id: HashMap<String, usize> = HashMap::new();
    for (i, f) in feats.iter().enumerate() {
        if let Some(id) = f.get("ID") {
            by_id.entry(id).or_insert(i);
        }
    }

    let mut genes: Vec<Gene> = Vec::new();
    let mut gene_of_feat: HashMap<usize, usize> = HashMap::new(); // feature index -> genes index
    let mut label_owner: HashMap<String, usize> = HashMap::new();

    for (i, f) in feats.iter().enumerate() {
        let Some(id) = f.get("ID") else { continue };
        if by_id.get(&id) != Some(&i) {
            continue; // repeated ID on a non-child feature: first occurrence only
        }
        let Some(ch) = children.get(&id) else { continue };
        let tx_strand = strand_of(f)?;
        let gene_fi = f
            .parents()
            .into_iter()
            .filter_map(|p| by_id.get(&p).copied())
            .find(|&p| p != i)
            .unwrap_or(i);
        let gi = match gene_of_feat.get(&gene_fi) {
            Some(&g) => g,
            None => {
                let gf = &feats[gene_fi];
                let strand = if gene_fi == i { tx_strand } else { strand_of(gf)? };
                let label = gf.label().unwrap_or_else(|| id.clone());
                let gi = genes.len();
                if label_owner.insert(label.clone(), gi).is_some() {
                    return Err(GffError::DuplicateGeneLabel(label));
                }
                genes.push(Gene { label, contig: gf.contig.clone(), strand, transcripts: Vec::new() });
                gene_of_feat.insert(gene_fi, gi);
                gi
            }
        };
        genes[gi].transcripts.push(Transcript {
            id,
            kind: f.kind.clone(),
            exons: finish(ch.exons.clone()),
            cds: finish(ch.cds.clone()),
            utr5: finish(ch.utr5.clone()),
            utr3: finish(ch.utr3.clone()),
        });
    }
    Ok(GeneModel { genes })
}
