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
    /// An exon/CDS/UTR line whose `Parent` is missing (`parent` empty) or matches no feature. `line` is 1-based.
    UnresolvedParent { line: usize, parent: String },
    /// A transcript's strand differs from the strand of its gene `parent`. `line` is the
    /// transcript's line (1-based).
    StrandMismatch { line: usize, parent: String },
    /// A transcript's (or child's) seqid differs from the seqid of its parent feature `parent`.
    /// `line` is the line of the transcript or child (1-based).
    ContigMismatch { line: usize, parent: String },
    /// An ID is repeated on features that differ in type, seqid or strand. (Lines of one
    /// discontinuous feature share an ID and agree in all three; that stays legal.) `line` is the
    /// later line (1-based).
    DuplicateId { id: String, line: usize },
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
            GffError::UnresolvedParent { line, parent } => {
                if parent.is_empty() {
                    write!(f, "GFF3 line {line}: child feature has no Parent")
                } else {
                    write!(f, "GFF3 line {line}: Parent '{parent}' matches no feature")
                }
            }
            GffError::StrandMismatch { line, parent } => {
                write!(f, "GFF3 line {line}: strand differs from the strand of its gene '{parent}'")
            }
            GffError::ContigMismatch { line, parent } => {
                write!(f, "GFF3 line {line}: seqid differs from the seqid of its parent '{parent}'")
            }
            GffError::DuplicateId { id, line } => write!(
                f,
                "GFF3 line {line}: ID '{id}' is used before on a feature of another type, seqid or strand"
            ),
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

fn is_exon_kind(k: &str) -> bool {
    k.ends_with("exon")
}

fn is_child_kind(k: &str) -> bool {
    is_exon_kind(k) || matches!(k, "CDS" | "five_prime_UTR" | "three_prime_UTR")
}

fn is_gene_kind(k: &str) -> bool {
    k.ends_with("gene")
}

struct ChildRec {
    line: usize,
    kind: String,
    contig: String,
    parents: Vec<String>,
    start: u64,
    end: u64,
}

struct Builder<'a> {
    feats: &'a [Feat],
    genes: Vec<Gene>,
    gene_of_feat: HashMap<usize, usize>,
    label_owner: HashMap<String, usize>,
}

impl Builder<'_> {
    /// Find or create the gene built from feature `fi`.
    fn gene_for(&mut self, fi: usize) -> Result<usize, GffError> {
        if let Some(&g) = self.gene_of_feat.get(&fi) {
            return Ok(g);
        }
        let f = &self.feats[fi];
        let strand = strand_of(f)?;
        let label = f.label().unwrap_or_else(|| format!("{}@line{}", f.kind, f.line));
        let gi = self.genes.len();
        if self.label_owner.insert(label.clone(), gi).is_some() {
            return Err(GffError::DuplicateGeneLabel(label));
        }
        self.genes.push(Gene { label, contig: f.contig.clone(), strand, transcripts: Vec::new() });
        self.gene_of_feat.insert(fi, gi);
        Ok(gi)
    }
}

pub fn parse_gff3<R: BufRead>(reader: R) -> Result<GeneModel, GffError> {
    // Pass 1: read lines into features and child records.
    let mut feats: Vec<Feat> = Vec::new();
    let mut child_recs: Vec<ChildRec> = Vec::new();
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
            let parents = feat.parents();
            child_recs.push(ChildRec { line: lineno, kind: feat.kind, contig: feat.contig, parents, start, end });
        } else {
            feats.push(feat);
        }
    }

    // Pass 2: index features by ID (first occurrence wins) and attach children.
    // A repeated ID is legal for the lines of one discontinuous feature (same type, seqid and
    // strand); otherwise it is `DuplicateId` (R29-6).
    let mut by_id: HashMap<String, usize> = HashMap::new();
    for (i, f) in feats.iter().enumerate() {
        if let Some(id) = f.get("ID") {
            match by_id.get(&id) {
                Some(&first) => {
                    let g = &feats[first];
                    if g.kind != f.kind || g.contig != f.contig || g.strand != f.strand {
                        return Err(GffError::DuplicateId { id, line: f.line });
                    }
                }
                None => {
                    by_id.insert(id, i);
                }
            }
        }
    }
    let mut children: HashMap<usize, Children> = HashMap::new();
    for c in &child_recs {
        if c.parents.is_empty() {
            return Err(GffError::UnresolvedParent { line: c.line, parent: String::new() });
        }
        for p in &c.parents {
            let Some(&fi) = by_id.get(p) else {
                return Err(GffError::UnresolvedParent { line: c.line, parent: p.clone() });
            };
            if feats[fi].contig != c.contig {
                return Err(GffError::ContigMismatch { line: c.line, parent: p.clone() });
            }
            let ch = children.entry(fi).or_default();
            let list = match c.kind.as_str() {
                "CDS" => &mut ch.cds,
                "five_prime_UTR" => &mut ch.utr5,
                "three_prime_UTR" => &mut ch.utr3,
                _ => &mut ch.exons,
            };
            list.push((c.start, c.end));
        }
    }

    // Pass 3: build genes and transcripts in file order.
    let mut b = Builder { feats: &feats, genes: Vec::new(), gene_of_feat: HashMap::new(), label_owner: HashMap::new() };
    for (i, f) in feats.iter().enumerate() {
        let id = f.get("ID");
        if let Some(id) = &id {
            if by_id.get(id) != Some(&i) {
                continue; // repeated ID: first occurrence only
            }
        }
        let has_children = children.contains_key(&i);
        let gene_parent = f
            .parents()
            .into_iter()
            .filter_map(|p| by_id.get(&p).copied())
            .find(|&p| p != i && is_gene_kind(&feats[p].kind));
        let gene_kind = is_gene_kind(&f.kind);
        if gene_kind {
            let gi = b.gene_for(i)?;
            if !has_children {
                continue;
            }
            // A gene-type feature that owns exons itself also yields a transcript.
            push_transcript(&mut b.genes[gi], f, id, children.get(&i));
            continue;
        }
        if !has_children && gene_parent.is_none() {
            continue; // region, repeat_region, or other non-transcript feature
        }
        if has_children {
            strand_of(f)?; // transcript strand first, then its gene's
        }
        let gi = match gene_parent {
            Some(gp) => b.gene_for(gp)?,
            None => b.gene_for(i)?, // own gene; also checks this feature's strand
        };
        if let Some(gp) = gene_parent {
            // R29-6: a transcript lies on its gene's seqid and strand. The gene's strand was
            // checked by `gene_for`; a childless transcript with no valid strand is not compared.
            let g = &feats[gp];
            let parent = || g.get("ID").unwrap_or_default();
            if g.contig != f.contig {
                return Err(GffError::ContigMismatch { line: f.line, parent: parent() });
            }
            if strand_of(f).is_ok_and(|s| s != b.genes[gi].strand) {
                return Err(GffError::StrandMismatch { line: f.line, parent: parent() });
            }
        }
        push_transcript(&mut b.genes[gi], f, id, children.get(&i));
    }
    Ok(GeneModel { genes: b.genes })
}

fn push_transcript(gene: &mut Gene, f: &Feat, id: Option<String>, ch: Option<&Children>) {
    let empty = Children::default();
    let ch = ch.unwrap_or(&empty);
    gene.transcripts.push(Transcript {
        id: id.unwrap_or_else(|| format!("{}@line{}", f.kind, f.line)),
        kind: f.kind.clone(),
        exons: finish(ch.exons.clone()),
        cds: finish(ch.cds.clone()),
        utr5: finish(ch.utr5.clone()),
        utr3: finish(ch.utr3.clone()),
    });
}
