//! `qrnamap de`: run the generated DESeq2 script with `Rscript`.
//!
//! Safety contract (rulings R18, R25): the script holds no paths or user text. The counts file,
//! the coldata file and the output prefix go to `Rscript` as separate argv elements through
//! `std::process::Command`. R cannot read zstd, so the counts go to a plain file in the scratch
//! directory first. Fastq paths are never written to the coldata file and never printed.

use std::fmt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use qrnamap_core::de::{check_design, coldata_tsv, deseq2_script, DeError, Shrink};
use qrnamap_core::{ColType, Contrast, CountMatrix, Design, SampleSheet};

use crate::exec::tail;
use crate::pipeline::{read_table, write_output, PipelineError, Scratch};

/// A problem in the count matrix file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MatrixReadError {
    /// The file has no lines.
    Empty,
    /// The header is not `gene<TAB>sample...` with at least one sample.
    BadHeader,
    /// A sample name occurs twice in the header.
    DuplicateSample(String),
    /// A gene id occurs in two rows.
    DuplicateGene(String),
    /// A row has another number of fields than the header. `line` is 1-based.
    BadRow { line: usize },
    /// A count is not a non-negative integer.
    BadCount { line: usize, column: String },
}

impl fmt::Display for MatrixReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatrixReadError::Empty => write!(f, "the count matrix is empty"),
            MatrixReadError::BadHeader => {
                write!(f, "the count matrix header must be 'gene', then one or more sample names, separated by tabs")
            }
            MatrixReadError::DuplicateSample(s) => write!(f, "sample '{s}' occurs twice in the count matrix header"),
            MatrixReadError::DuplicateGene(g) => write!(f, "gene '{g}' occurs in more than one row of the count matrix"),
            MatrixReadError::BadRow { line } => {
                write!(f, "count matrix line {line} has a different number of fields than the header")
            }
            MatrixReadError::BadCount { line, column } => {
                write!(f, "count matrix line {line}, column '{column}': the count is not a non-negative integer")
            }
        }
    }
}

impl std::error::Error for MatrixReadError {}

/// Parse the `counts.tsv` that `qrnamap matrix` writes. Blank lines and CRLF are accepted.
pub fn read_count_matrix(text: &str) -> Result<CountMatrix, MatrixReadError> {
    let mut lines = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.strip_suffix('\r').unwrap_or(l)))
        .filter(|(_, l)| !l.trim().is_empty());
    let (_, header) = lines.next().ok_or(MatrixReadError::Empty)?;
    let h: Vec<&str> = header.split('\t').collect();
    if h.len() < 2 || h[0] != "gene" || h[1..].iter().any(|s| s.is_empty()) {
        return Err(MatrixReadError::BadHeader);
    }
    let samples: Vec<String> = h[1..].iter().map(|s| s.to_string()).collect();
    let mut seen = std::collections::HashSet::new();
    for s in &samples {
        if !seen.insert(s.as_str()) {
            return Err(MatrixReadError::DuplicateSample(s.clone()));
        }
    }
    let mut genes = Vec::new();
    let mut counts = Vec::new();
    let mut seen_genes = std::collections::HashSet::new();
    for (line, row) in lines {
        let f: Vec<&str> = row.split('\t').collect();
        if f.len() != h.len() || f[0].is_empty() {
            return Err(MatrixReadError::BadRow { line });
        }
        let mut r = Vec::with_capacity(samples.len());
        for (j, cell) in f[1..].iter().enumerate() {
            if cell.is_empty() || !cell.bytes().all(|b| b.is_ascii_digit()) {
                return Err(MatrixReadError::BadCount { line, column: samples[j].clone() });
            }
            r.push(cell.parse::<u64>().map_err(|_| MatrixReadError::BadCount { line, column: samples[j].clone() })?);
        }
        if !seen_genes.insert(f[0].to_string()) {
            return Err(MatrixReadError::DuplicateGene(f[0].to_string()));
        }
        genes.push(f[0].to_string());
        counts.push(r);
    }
    Ok(CountMatrix { genes, samples, counts })
}

/// Arguments of one DE run.
#[derive(Debug, Clone)]
pub struct DeArgs {
    pub counts: PathBuf,
    pub samples: PathBuf,
    pub design: String,
    /// `factor,level,reference`.
    pub contrast: String,
    pub shrink: Option<Shrink>,
    /// Output prefix. The files are `<prefix>.results.csv` and so on.
    pub out: PathBuf,
    pub tmpdir: Option<PathBuf>,
}

/// What a successful run reports.
#[derive(Debug, Clone)]
pub struct DeOutput {
    pub deseq2_version: String,
    pub files: Vec<PathBuf>,
}

#[derive(Debug)]
pub enum DeRunError {
    /// The sample sheet is invalid. Holds one line per problem; it never holds fastq paths.
    Sheet(String),
    /// The design or contrast text is invalid.
    Design(String),
    Matrix(MatrixReadError),
    Pipeline(PipelineError),
    /// An input check failed before R started, or R reported a rank-deficient design.
    De(DeError),
    Io { path: PathBuf, source: std::io::Error },
    /// `Rscript` could not start.
    Spawn(std::io::Error),
    /// An R package needed for the run is not installed.
    MissingPackage(String),
    /// R exited with an error. `stderr_tail` is bounded.
    RFailed { status: String, stderr_tail: String },
    /// R exited with status 0 but did not write this file.
    MissingOutput(PathBuf),
}

impl fmt::Display for DeRunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DeRunError::Sheet(m) => write!(f, "{m}"),
            DeRunError::Design(m) => write!(f, "{m}"),
            DeRunError::Matrix(e) => write!(f, "{e}"),
            DeRunError::Pipeline(e) => write!(f, "{e}"),
            DeRunError::De(e) => write!(f, "{e}"),
            DeRunError::Io { path, source } => write!(f, "{}: {source}", path.display()),
            DeRunError::Spawn(e) => write!(f, "cannot start 'Rscript': {e} (is R on PATH?)"),
            DeRunError::MissingPackage(p) => write!(
                f,
                "the R package '{p}' is not installed; install it yourself and run again. Unshrunk results, if any, were kept"
            ),
            DeRunError::RFailed { status, stderr_tail } => {
                write!(f, "'Rscript' failed ({status}). stderr (last {} bytes at most):\n{stderr_tail}", crate::exec::TAIL_BYTES)
            }
            DeRunError::MissingOutput(p) => write!(f, "R exited with status 0 but did not write {}", p.display()),
        }
    }
}

impl std::error::Error for DeRunError {}

impl From<DeError> for DeRunError {
    fn from(e: DeError) -> Self {
        DeRunError::De(e)
    }
}

impl From<PipelineError> for DeRunError {
    fn from(e: PipelineError) -> Self {
        DeRunError::Pipeline(e)
    }
}

/// Removes the scratch subdirectory when dropped.
struct Cleanup(PathBuf);

impl Drop for Cleanup {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Classify R's stderr after a failure.
fn classify(status: String, stderr: &[u8]) -> DeRunError {
    let text = String::from_utf8_lossy(stderr);
    if let Some(i) = text.find("MISSING_PACKAGE:") {
        let rest = &text[i + "MISSING_PACKAGE:".len()..];
        let pkg: String = rest.chars().take_while(|c| c.is_ascii_alphanumeric() || *c == '.').collect();
        return DeRunError::MissingPackage(pkg);
    }
    if text.contains("there is no package called") && text.contains("DESeq2") {
        return DeRunError::MissingPackage("DESeq2".into());
    }
    if text.contains("not full rank") {
        return DeRunError::De(DeError::NotFullRank(tail(stderr)));
    }
    DeRunError::RFailed { status, stderr_tail: tail(stderr) }
}

fn abs(p: &Path) -> Result<PathBuf, DeRunError> {
    std::path::absolute(p).map_err(|source| DeRunError::Io { path: p.to_path_buf(), source })
}

/// Forward at most this many bytes of R's stderr after a successful run.
const FORWARD_BYTES: usize = 8192;

/// The last `FORWARD_BYTES` of `bytes` as text, cut on a UTF-8 boundary, each line prefixed `R: `.
fn forward_tail(bytes: &[u8]) -> String {
    let mut start = bytes.len().saturating_sub(FORWARD_BYTES);
    while start < bytes.len() && (bytes[start] & 0b1100_0000) == 0b1000_0000 {
        start += 1;
    }
    let cut = start > 0;
    let text = String::from_utf8_lossy(&bytes[start..]);
    // A cut usually lands inside a line; drop the partial first line.
    let text: &str = if cut { text.split_once('\n').map(|x| x.1).unwrap_or("") } else { &text };
    text.lines().filter(|l| !l.trim().is_empty()).map(|l| format!("R: {l}\n")).collect()
}

/// Remove the output files of an earlier run (and their `.partial` files) so that a run that
/// stops early leaves no stale results. Only NotFound is ignored.
pub fn remove_old_outputs(prefix: &Path) -> Result<(), DeRunError> {
    for suffix in [".results.csv", ".normalized.csv", ".vst.csv", ".versions.txt", ".results.shrunk.csv"] {
        for extra in ["", ".partial"] {
            let mut p = prefix.as_os_str().to_owned();
            p.push(suffix);
            p.push(extra);
            let p = PathBuf::from(p);
            match std::fs::remove_file(&p) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(source) => return Err(DeRunError::Io { path: p, source }),
            }
        }
    }
    Ok(())
}

pub fn run_de(args: &DeArgs) -> Result<DeOutput, DeRunError> {
    let prefix = abs(&args.out)?;
    remove_old_outputs(&prefix)?;
    let sheet_text = read_table(&args.samples)?;
    let sheet = SampleSheet::parse(&sheet_text).map_err(|errs| {
        let lines: Vec<String> = errs.iter().map(|e| format!("  {e}")).collect();
        DeRunError::Sheet(format!("sample sheet {}:\n{}", args.samples.display(), lines.join("\n")))
    })?;
    let design = Design::parse(&args.design, &sheet.columns).map_err(|e| DeRunError::Design(format!("design: {e}")))?;
    let contrast =
        Contrast::parse(&args.contrast, &sheet.columns).map_err(|e| DeRunError::Design(format!("contrast: {e}")))?;
    if !design.terms().contains(&contrast.factor) {
        return Err(DeError::ContrastNotInDesign(contrast.factor).into());
    }
    let counts_text = read_table(&args.counts)?;
    let matrix = read_count_matrix(&counts_text).map_err(DeRunError::Matrix)?;
    check_design(&sheet, &contrast, &matrix)?;

    let columns: Vec<(String, ColType)> = {
        let mut v: Vec<(String, ColType)> = Vec::new();
        for t in design.terms() {
            if !v.iter().any(|(n, _)| n == t) {
                let ty = sheet.columns.iter().find(|(n, _)| n == t).map(|(_, ty)| *ty).unwrap_or(ColType::Factor);
                v.push((t.clone(), ty));
            }
        }
        v
    };

    let scratch = Scratch::from_flag_or_env(args.tmpdir.clone());
    let dir = scratch.dir("de")?;
    let _guard = Cleanup(dir.clone());
    let counts_path = dir.join("counts.tsv");
    let coldata_path = dir.join("coldata.tsv");
    let script_path = dir.join("de.R");
    write_output(&counts_path, matrix.to_tsv().as_bytes())?;
    write_output(&coldata_path, coldata_tsv(&sheet, &design).as_bytes())?;
    write_output(&script_path, deseq2_script(&design, &contrast, &columns, args.shrink).as_bytes())?;

    if let Some(parent) = prefix.parent() {
        std::fs::create_dir_all(parent).map_err(|source| DeRunError::Io { path: parent.to_path_buf(), source })?;
    }
    let out = Command::new("Rscript")
        .arg("--vanilla")
        .arg(&script_path)
        .arg(&counts_path)
        .arg(&coldata_path)
        .arg(&prefix)
        .current_dir(&dir)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .map_err(DeRunError::Spawn)?;
    if !out.status.success() {
        return Err(classify(out.status.to_string(), &out.stderr));
    }
    eprint!("{}", forward_tail(&out.stderr));
    let with = |suffix: &str| {
        let mut s = prefix.as_os_str().to_owned();
        s.push(suffix);
        PathBuf::from(s)
    };
    let mut files = vec![with(".results.csv"), with(".normalized.csv"), with(".vst.csv"), with(".versions.txt")];
    if args.shrink.is_some() {
        files.push(with(".results.shrunk.csv"));
    }
    for f in &files {
        if !f.is_file() {
            return Err(DeRunError::MissingOutput(f.clone()));
        }
    }
    let versions = std::fs::read_to_string(with(".versions.txt"))
        .map_err(|source| DeRunError::Io { path: with(".versions.txt"), source })?;
    let deseq2_version = versions
        .lines()
        .find_map(|l| l.strip_prefix("DESeq2 "))
        .unwrap_or("unknown")
        .trim()
        .to_string();
    drop(_guard);
    scratch.cleanup();
    Ok(DeOutput { deseq2_version, files })
}
