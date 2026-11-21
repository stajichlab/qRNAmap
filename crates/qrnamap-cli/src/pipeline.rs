//! Native pipelines: transcripts, kallisto quant, HISAT2 + featureCounts, gene matrix and the
//! strand check.
//!
//! `quant_sample`, `align_sample` and `strand_check_sample` are independent per-sample
//! functions. Each takes everything it needs in a context struct, writes only
//! `<out>/<sample>/` and `<out>/<sample>.genecounts.tsv.zst`, and uses its own scratch
//! subdirectory, so a scheduler can run samples as separate jobs.
//!
//! Safety contract (ruling R18): FASTQ paths and sample names reach tools only as argv
//! elements (through `exec::run_call`). They are never pasted into shell or R text.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;
use std::fs::{self, File};
use std::io::{self, BufRead, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use qrnamap_core::{
    assemble, assigned_fraction, build_transcripts, feature_counts, gene_counts_from_abundance, gene_table,
    gene_counts_from_salmon, gene_table_tsv,
    hisat2_align, hisat2_build, kallisto_index, kallisto_quant, parse_featurecounts,
    parse_featurecounts_summary, saf, salmon_index, salmon_quant, splice_sites, tx2gene_tsv,
    unlisted_genes, write_fasta,
    BuildError, ConfigError, CountError, CountMatrix, Mode, Quantifier, RunConfig, Sample,
    StrandReport, Strandedness, ToolCall,
};
use qseq_io::{parse_gff3, read_fasta, GeneModel};

use crate::exec::{run_call, run_call_stderr, ExecError};
use crate::io::{file_format, open_input, Format, InputError};

/// The assigned (or pseudoaligned) fraction below which the CLI prints a warning after a run.
/// 0.5 is a chosen value. UNMEASURED: not calibrated on real libraries. A wrongly set
/// strandedness gives a fraction near 0 on the simulated data in the tests.
pub const LOW_ASSIGNED_FRACTION: f64 = 0.5;

/// At most this many gene names are listed in one warning.
const MAX_LISTED: usize = 20;

#[derive(Debug)]
pub enum PipelineError {
    /// A tool-call builder refused the sample (ruling R20). Raised before any tool runs.
    Config { sample: String, source: ConfigError },
    /// The sample name equals a name qrnamap writes in the output directory
    /// (`RESERVED_OUTPUT_NAMES`, or one of them with `.partial`). Raised before any tool runs.
    ReservedSampleName { sample: String, reserved: String },
    /// A tool failed. `sample` is `None` for index builds.
    Exec { sample: Option<String>, source: ExecError },
    /// Opening or reading an input file failed (includes truncated or corrupt compression).
    /// `sample` names the sample when the file is one of its FASTQ files.
    Input { sample: Option<String>, source: InputError },
    /// Reading or writing a file failed.
    Io { path: PathBuf, source: io::Error },
    /// The GFF3 file cannot be parsed.
    Gff { path: PathBuf, detail: String },
    /// Transcript building failed (for example, a gene on a contig that is not in the FASTA).
    Build(BuildError),
    /// A count table cannot be parsed.
    Count { path: PathBuf, source: CountError },
    /// A text table has a bad line. `line` is 1-based; 0 means the whole file.
    BadTable { path: PathBuf, line: usize, detail: String },
    /// A field is missing from a tool's report (`run_info.json`, featureCounts summary).
    MissingField { path: PathBuf, field: String },
    /// A temporary directory is needed but neither `--tmpdir` nor `$SCRATCH` is set.
    NoScratch,
    /// A gene label (or transcript id) contains a tab, newline or carriage return.
    BadGeneLabel(String),
    /// A path is not valid UTF-8, so it cannot become a tool argument.
    NonUtf8Path(PathBuf),
    /// The quantifier (`tool`: kallisto or salmon) mapped no read. kallisto 0.52.0 then exits
    /// with an error; salmon 2.8.0 exited 0 in a manual check. The usual causes are a wrong
    /// strandedness setting or reads that do not belong to the index.
    NoReadsPseudoaligned {
        sample: String,
        tool: &'static str,
        strandedness: Strandedness,
        stderr_tail: String,
    },
    /// A gene count file's header names another sample (for example a file copied or left
    /// from another run).
    CountFileSample { path: PathBuf, expected: String, found: String },
    /// A sample has no gene count file (its run failed or did not happen).
    MissingCountFile { sample: String, path: PathBuf },
    /// Mate 1 and mate 2 of one lane hold different numbers of FASTQ records. `lane` is 1-based.
    MateCountMismatch { sample: String, lane: usize, r1: u64, r2: u64 },
    /// A FASTQ file of the sample holds no record.
    NoReads { sample: String, file: PathBuf },
    /// The quantifier's processed count (kallisto `n_processed`, salmon `num_processed`)
    /// differs from the pairs (reads, for single-end) counted in the input files.
    PairCountMismatch { sample: String, input: u64, processed: u64, tool: &'static str },
    /// A sample has no `<sample>.summary.tsv.zst` (its run failed or did not happen).
    MissingSummaryFile { sample: String, path: PathBuf },
    /// A per-sample summary file names another sample, or is not a one-row summary table.
    SummaryFileSample { path: PathBuf, expected: String, found: String },
    /// Genes of the tx2gene table are not genes of the GFF3 given for the gene table, so the
    /// gene table would not match the count rows. `missing` holds the first few names (sorted),
    /// `n_missing` the total.
    GenesTableMismatch { missing: Vec<String>, n_missing: usize },
}

impl fmt::Display for PipelineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PipelineError::Config { sample, source } => {
                write!(f, "sample '{sample}': {source}")
            }
            PipelineError::ReservedSampleName { sample, reserved } => write!(
                f,
                "sample '{sample}': the name is reserved: qrnamap writes '{reserved}' in the \
                 output directory. Rename the sample in the sheet"
            ),
            PipelineError::Exec { sample: Some(s), source } => write!(f, "sample '{s}': {source}"),
            PipelineError::Exec { sample: None, source } => write!(f, "{source}"),
            PipelineError::Input { sample: Some(s), source } => write!(f, "sample '{s}': {source}"),
            PipelineError::Input { sample: None, source } => write!(f, "{source}"),
            PipelineError::Io { path, source } => write!(f, "{}: {source}", path.display()),
            PipelineError::Gff { path, detail } => {
                write!(f, "GFF3 {}: {detail}", path.display())
            }
            PipelineError::Build(e) => write!(f, "{e}"),
            PipelineError::Count { path, source } => write!(f, "{}: {source}", path.display()),
            PipelineError::BadTable { path, line: 0, detail } => {
                write!(f, "{}: {detail}", path.display())
            }
            PipelineError::BadTable { path, line, detail } => {
                write!(f, "{} line {line}: {detail}", path.display())
            }
            PipelineError::MissingField { path, field } => {
                write!(f, "{}: field '{field}' not found", path.display())
            }
            PipelineError::NoScratch => write!(
                f,
                "a temporary directory is needed: pass --tmpdir DIR, or set $SCRATCH or $TMPDIR"
            ),
            PipelineError::CountFileSample { path, expected, found } => write!(
                f,
                "{}: the header names sample '{found}', expected '{expected}'",
                path.display()
            ),
            PipelineError::MissingCountFile { sample, path } => write!(
                f,
                "sample '{sample}': gene count file {} does not exist (run quant or align for this sample)",
                path.display()
            ),
            PipelineError::MissingSummaryFile { sample, path } => write!(
                f,
                "sample '{sample}': summary file {} does not exist; rerun quant or align for sample '{sample}'",
                path.display()
            ),
            PipelineError::SummaryFileSample { path, expected, found } => write!(
                f,
                "{}: the file names sample '{found}', expected '{expected}'; rerun quant or align for sample '{expected}'",
                path.display()
            ),
            PipelineError::BadGeneLabel(l) => write!(
                f,
                "gene label or transcript id {l:?} contains a tab, newline or carriage return"
            ),
            PipelineError::NonUtf8Path(p) => {
                write!(f, "path {} is not valid UTF-8", p.display())
            }
            PipelineError::MateCountMismatch { sample, lane, r1, r2 } => write!(
                f,
                "sample '{sample}': lane {lane}: mate 1 holds {r1} FASTQ records and mate 2 holds {r2}; the mate files must hold the same reads"
            ),
            PipelineError::NoReads { sample, file } => {
                write!(f, "sample '{sample}': {} holds no FASTQ record", file.display())
            }
            PipelineError::PairCountMismatch { sample, input, processed, tool } => write!(
                f,
                "sample '{sample}': the input files hold {input} read pairs \
                 (reads, for single-end), but {tool} processed {processed}"
            ),
            PipelineError::NoReadsPseudoaligned { sample, tool, strandedness, stderr_tail } => {
                // kallisto text is unchanged; salmon calls it mapping.
                let verb = if *tool == "salmon" { "mapped" } else { "pseudoaligned" };
                write!(
                    f,
                    "sample '{sample}': {tool} {verb} 0 reads ({verb} fraction 0). \
                     Check the strandedness setting (now {strandedness:?}) with \
                     `qrnamap strand-check`, and check that the reads match the transcripts. \
                     {tool} stderr:\n{stderr_tail}"
                )
            }
            PipelineError::GenesTableMismatch { missing, n_missing } => {
                let more = if *n_missing > missing.len() { ", ..." } else { "" };
                write!(
                    f,
                    "{n_missing} gene(s) of the tx2gene table are not genes of the --gff GFF3: {}{more}. \
                     The gene table would not match the count rows; give the GFF3 the transcripts were built from",
                    missing.join(", ")
                )
            }
        }
    }
}

impl std::error::Error for PipelineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            PipelineError::Config { source, .. } => Some(source),
            PipelineError::Exec { source, .. } => Some(source),
            PipelineError::Input { source, .. } => Some(source),
            PipelineError::Io { source, .. } => Some(source),
            PipelineError::Build(e) => Some(e),
            PipelineError::Count { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// An input error that is not tied to one sample.
fn input_err(source: InputError) -> PipelineError {
    PipelineError::Input { sample: None, source }
}

fn io_err(path: &Path) -> impl FnOnce(io::Error) -> PipelineError + '_ {
    move |source| PipelineError::Io { path: path.to_path_buf(), source }
}

// ---------------------------------------------------------------------------------------------
// Scratch space

/// Where temporary files go. `dir(sub)` creates `<base>/qrnamap-<pid>/<sub>`.
#[derive(Debug, Clone)]
pub struct Scratch {
    base: Option<PathBuf>,
}

impl Scratch {
    pub fn new(base: impl Into<PathBuf>) -> Self {
        Scratch { base: Some(base.into()) }
    }

    /// No temporary directory. `dir` then fails with `NoScratch`.
    pub fn none() -> Self {
        Scratch { base: None }
    }

    /// Ruling R23: `--tmpdir` when given, otherwise `$SCRATCH`, otherwise `$TMPDIR` (each
    /// when set and non-empty), otherwise none.
    pub fn from_flag_or_env(flag: Option<PathBuf>) -> Self {
        let env = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        Scratch { base: flag.or_else(|| env("SCRATCH")).or_else(|| env("TMPDIR")) }
    }

    fn run_dir(&self) -> Result<PathBuf, PipelineError> {
        let base = self.base.as_ref().ok_or(PipelineError::NoScratch)?;
        Ok(base.join(format!("qrnamap-{}", std::process::id())))
    }

    /// Create and return `<base>/qrnamap-<pid>/<sub>`.
    pub fn dir(&self, sub: &str) -> Result<PathBuf, PipelineError> {
        let d = absolute(&self.run_dir()?.join(sub))?;
        fs::create_dir_all(&d).map_err(io_err(&d))?;
        Ok(d)
    }

    /// Remove `<base>/qrnamap-<pid>` when it is empty. Errors are ignored.
    pub fn cleanup(&self) {
        if let Ok(d) = self.run_dir() {
            let _ = fs::remove_dir(d);
        }
    }
}

/// Removes a scratch subdirectory when dropped, also on error paths.
struct ScratchGuard(Option<PathBuf>);

impl Drop for ScratchGuard {
    fn drop(&mut self) {
        if let Some(d) = self.0.take() {
            let _ = fs::remove_dir_all(d);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Small file helpers

fn absolute(p: &Path) -> Result<PathBuf, PipelineError> {
    std::path::absolute(p).map_err(io_err(p))
}

fn utf8(p: &Path) -> Result<&str, PipelineError> {
    p.to_str().ok_or_else(|| PipelineError::NonUtf8Path(p.to_path_buf()))
}

/// Read a whole plain, gzip or zstd text file.
pub fn read_table(path: &Path) -> Result<String, PipelineError> {
    let mut r = open_input(path).map_err(input_err)?;
    let mut s = String::new();
    r.read_to_string(&mut s).map_err(io_err(path))?;
    Ok(s)
}

enum Sink {
    Plain(BufWriter<File>),
    Gz(flate2::write::GzEncoder<BufWriter<File>>),
    Zst(zstd::Encoder<'static, BufWriter<File>>),
}

/// An output file, compressed by its extension: `.zst` zstd, `.gz` gzip, anything else plain.
/// Data goes to `<path>.partial`; `finish` renames it to `path`, so a failed run leaves no
/// complete-looking file.
pub struct Output {
    sink: Sink,
    path: PathBuf,
    partial: PathBuf,
}

impl Output {
    pub fn create(path: &Path) -> Result<Output, PipelineError> {
        let mut partial = path.as_os_str().to_owned();
        partial.push(".partial");
        let partial = PathBuf::from(partial);
        let file = BufWriter::new(File::create(&partial).map_err(io_err(&partial))?);
        let name = path.to_string_lossy();
        let sink = if name.ends_with(".zst") {
            Sink::Zst(zstd::Encoder::new(file, 3).map_err(io_err(path))?)
        } else if name.ends_with(".gz") {
            Sink::Gz(flate2::write::GzEncoder::new(file, flate2::Compression::default()))
        } else {
            Sink::Plain(file)
        };
        Ok(Output { sink, path: path.to_path_buf(), partial })
    }

    pub fn finish(self) -> Result<(), PipelineError> {
        let e = io_err(&self.path);
        let mut file = match self.sink {
            Sink::Plain(w) => w,
            Sink::Gz(w) => w.finish().map_err(e)?,
            Sink::Zst(w) => w.finish().map_err(e)?,
        };
        file.flush().map_err(io_err(&self.path))?;
        fs::rename(&self.partial, &self.path).map_err(io_err(&self.path))
    }
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.sink {
            Sink::Plain(w) => w.write(buf),
            Sink::Gz(w) => w.write(buf),
            Sink::Zst(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match &mut self.sink {
            Sink::Plain(w) => w.flush(),
            Sink::Gz(w) => w.flush(),
            Sink::Zst(w) => w.flush(),
        }
    }
}

/// Write `data` to `path`, compressed by extension.
pub fn write_output(path: &Path, data: &[u8]) -> Result<(), PipelineError> {
    let mut out = Output::create(path)?;
    out.write_all(data).map_err(io_err(path))?;
    out.finish()
}

/// Copy `src` (any compression) to `dst` (compressed by extension of `dst`).
fn copy_converted(src: &Path, dst: &Path) -> Result<(), PipelineError> {
    let mut r = open_input(src).map_err(input_err)?;
    let mut out = Output::create(dst)?;
    io::copy(&mut r, &mut out).map_err(io_err(src))?;
    out.finish()
}

fn check_label(l: &str) -> Result<(), PipelineError> {
    if l.contains(['\t', '\n', '\r']) {
        return Err(PipelineError::BadGeneLabel(l.to_string()));
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// Reference loading and `transcripts`

/// Contigs as (name, sequence), as `read_fasta` returns them.
type Genome = Vec<(String, Vec<u8>)>;

fn load_model(gff: &Path) -> Result<GeneModel, PipelineError> {
    let model = parse_gff3(open_input(gff).map_err(input_err)?).map_err(|e| {
        PipelineError::Gff { path: gff.to_path_buf(), detail: e.to_string() }
    })?;
    for g in &model.genes {
        check_label(&g.label)?;
        for t in &g.transcripts {
            check_label(&t.id)?;
        }
    }
    Ok(model)
}

fn load_reference(reference: &Path, gff: &Path) -> Result<(Genome, GeneModel), PipelineError> {
    let genome = read_fasta(open_input(reference).map_err(input_err)?)
        .map_err(io_err(reference))?;
    Ok((genome, load_model(gff)?))
}

/// Name of the gene table in a run directory.
pub const GENES_TABLE: &str = "genes.tsv.zst";

/// Write the gene table (class, biotype, mitochondrial flag per gene) to `path`. The suffix
/// picks `.zst`, `.gz` or plain. A gene label with a tab, newline or carriage return is
/// refused before anything is written. `mito_contigs` are extra contig names to flag; they
/// are only compared with the GFF seqids.
pub fn write_genes_table(path: &Path, model: &GeneModel, mito_contigs: &[String]) -> Result<(), PipelineError> {
    for g in &model.genes {
        check_label(&g.label)?;
    }
    let rows = gene_table(model, mito_contigs);
    write_output(path, gene_table_tsv(&rows).as_bytes())?;
    // R35-7: a name that matches no contig is a note, not an error (the match is exact).
    for c in unmatched_mito_contigs(model, mito_contigs) {
        eprintln!("note: --mito-contig '{c}' matches no contig of the GFF3 (the match is exact and case-sensitive); it flags no gene");
    }
    Ok(())
}

/// The `mito_contigs` names that are neither the contig of a gene nor the seqid of a
/// mitochondrial `region` feature of `model`, in the given order.
pub fn unmatched_mito_contigs(model: &GeneModel, mito_contigs: &[String]) -> Vec<String> {
    let known: HashSet<&str> =
        model.genes.iter().map(|g| g.contig.as_str()).chain(model.mito_contigs.iter().map(String::as_str)).collect();
    mito_contigs.iter().filter(|c| !known.contains(c.as_str())).cloned().collect()
}

/// At most this many gene names are listed in `GenesTableMismatch`.
const MAX_MISSING_LISTED: usize = 10;

/// Every gene of `tx2gene` must be a gene label of `model` (R35-3). More genes in the model is
/// fine: genes without intervals are absent from tx2gene (ruling R15).
pub fn check_tx2gene_in_model(tx2gene: &HashMap<String, String>, model: &GeneModel) -> Result<(), PipelineError> {
    let labels: HashSet<&str> = model.genes.iter().map(|g| g.label.as_str()).collect();
    let missing: BTreeSet<&str> = tx2gene.values().map(String::as_str).filter(|g| !labels.contains(g)).collect();
    if missing.is_empty() {
        return Ok(());
    }
    Err(PipelineError::GenesTableMismatch {
        n_missing: missing.len(),
        missing: missing.iter().take(MAX_MISSING_LISTED).map(|g| g.to_string()).collect(),
    })
}

/// Remove an old `<out>/genes.tsv.zst` (and its `.partial`) before a run (ruling R24).
pub fn remove_old_genes_table(out_dir: &Path) -> Result<(), PipelineError> {
    for name in [GENES_TABLE.to_string(), format!("{GENES_TABLE}.partial")] {
        let p = out_dir.join(name);
        match fs::remove_file(&p) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(PipelineError::Io { path: p, source: e }),
        }
    }
    Ok(())
}

/// Where `prepare_quant` takes the gene table from: the GFF3 and extra mitochondrial contigs.
#[derive(Debug, Clone, Copy)]
pub struct GenesSource<'a> {
    pub gff: &'a Path,
    pub mito_contigs: &'a [String],
}

/// `qrnamap transcripts`: transcript FASTA and tx2gene (and optionally the SAF and HISAT2
/// splice sites) from a genome and a GFF3. Inputs may be plain, gzip or zstd. Outputs are
/// compressed by extension. Returns the number of transcripts written.
pub fn run_transcripts(
    reference: &Path,
    gff: &Path,
    out: &Path,
    tx2gene: &Path,
    saf_out: Option<&Path>,
    splice_out: Option<&Path>,
) -> Result<usize, PipelineError> {
    run_transcripts_with_genes(reference, gff, out, tx2gene, saf_out, splice_out, None)
}

/// `run_transcripts`, and the gene table too when `genes` is `Some((path, mito_contigs))`.
pub fn run_transcripts_with_genes(
    reference: &Path,
    gff: &Path,
    out: &Path,
    tx2gene: &Path,
    saf_out: Option<&Path>,
    splice_out: Option<&Path>,
    genes: Option<(&Path, &[String])>,
) -> Result<usize, PipelineError> {
    let (genome, model) = load_reference(reference, gff)?;
    let records = build_transcripts(&genome, &model).map_err(PipelineError::Build)?;
    let mut fa = Output::create(out)?;
    write_fasta(&records, &mut fa).map_err(io_err(out))?;
    fa.finish()?;
    write_output(tx2gene, tx2gene_tsv(&records).as_bytes())?;
    if let Some(p) = saf_out {
        write_output(p, saf(&model).as_bytes())?;
    }
    if let Some(p) = splice_out {
        write_output(p, splice_sites(&model).as_bytes())?;
    }
    if let Some((p, mito)) = genes {
        write_genes_table(p, &model, mito)?;
    }
    Ok(records.len())
}

// ---------------------------------------------------------------------------------------------
// Sample checks and read staging

/// Run every tool-call builder for every sample with placeholder paths, so a bad sample is a
/// named error before any tool runs (ruling R20). Paths are not checked here. A sample name
/// in `RESERVED_OUTPUT_NAMES` is refused first.
pub fn check_samples(samples: &[Sample], cfg: &RunConfig, mode: Mode) -> Result<(), PipelineError> {
    for s in samples {
        if let Some(r) = reserved_name(&s.name) {
            return Err(PipelineError::ReservedSampleName { sample: s.name.clone(), reserved: r });
        }
        let r = match mode {
            Mode::Quant => quant_call(s, "idx", "out", cfg).map(|_| ()),
            Mode::Align => hisat2_align(s, "idx", "ss", "x.sam", cfg)
                .and_then(|_| feature_counts("x.sam", "x.saf", "fc", !s.fastq2.is_empty(), cfg))
                .map(|_| ()),
        };
        r.map_err(|source| PipelineError::Config { sample: s.name.clone(), source })?;
    }
    Ok(())
}

/// Names qrnamap writes directly in the output directory that are also valid sample names.
/// A sample with one of these names (or one of them with `.partial`) would collide with the
/// file or directory, because each sample gets the directory `<out>/<sample>/`. Checked in
/// every mode, so one output directory can hold both modes' files.
pub const RESERVED_OUTPUT_NAMES: &[&str] = &[
    "index",
    "salmon_index",
    "slurm",
    "transcripts.kidx",
    "transcripts.fa.zst",
    "tx2gene.tsv.zst",
    "genes.saf",
    "splice_sites.txt",
    GENES_TABLE,
    MAPPING_SUMMARY,
    "counts.tsv.zst",
    "strand_check.tsv",
    "deseq2.results.csv",
    "deseq2.normalized.csv",
    "deseq2.vst.csv",
    "deseq2.versions.txt",
    "deseq2.results.shrunk.csv",
];

/// The reserved output name that `sample` equals, if any.
fn reserved_name(sample: &str) -> Option<String> {
    let base = sample.strip_suffix(".partial").unwrap_or(sample);
    RESERVED_OUTPUT_NAMES.contains(&base).then(|| sample.to_string())
}

/// Writes every byte read through it to `out` (a staged copy), so one pass both decompresses a
/// file and counts its records.
struct Tee<'a> {
    inner: Box<dyn BufRead>,
    out: &'a mut Output,
}

impl Read for Tee<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.out.write_all(&buf[..n]).map_err(|e| {
            io::Error::new(e.kind(), format!("writing {}: {e}", self.out.path.display()))
        })?;
        Ok(n)
    }
}

/// Copy `src` (any compression) to the plain file `dst` and count its FASTQ records in the
/// same pass. A truncated or corrupt stream, or a bad record, is an `InputError` for `sample`.
fn copy_counting(src: &Path, dst: &Path, sample: &str) -> Result<u64, PipelineError> {
    let named = |source| PipelineError::Input { sample: Some(sample.to_string()), source };
    let r = open_input(src).map_err(named)?;
    let mut out = Output::create(dst)?;
    let n = crate::run::count_records(io::BufReader::new(Tee { inner: r, out: &mut out }), src, None)
        .map_err(named)?;
    out.finish()?;
    Ok(n)
}

/// FASTQ paths of a sample after staging, and (when validated) the record count of every
/// lane of mate 1 (`counts[0]`) and mate 2 (`counts[1]`).
struct Staged {
    sample: Sample,
    counts: [Vec<u64>; 2],
}

/// R29-11: every FASTQ file of every sample exists, is a file and can be opened for reading.
/// Called before any index build, so a typo in the sheet fails in seconds, not after the index.
/// The error names the sample and the file.
pub fn check_fastq_files(samples: &[Sample]) -> Result<(), PipelineError> {
    for s in samples {
        for f in s.fastq1.iter().chain(&s.fastq2) {
            let path = PathBuf::from(f);
            let r = File::open(&path).and_then(|h| h.metadata()).and_then(|m| {
                if m.is_file() {
                    Ok(())
                } else {
                    Err(io::Error::new(io::ErrorKind::InvalidInput, "not a regular file"))
                }
            });
            if let Err(source) = r {
                return Err(PipelineError::Input { sample: Some(s.name.clone()), source: InputError::Io { path, source } });
            }
        }
    }
    Ok(())
}

/// Make every FASTQ path absolute and decompress the ones the tool cannot read into
/// `scratch_dir`. `gzip_ok`: the tool reads gzip itself (kallisto). Plain files always pass.
/// zstd files (and gzip when `!gzip_ok`) are decompressed, because neither tool reads zstd and
/// native HISAT2 reading gzip is UNVERIFIED (spec 11f).
///
/// `validate` (final review C1): read every file to its end and count its records; a
/// decompressed file is counted while it is copied, so it is read once. A truncated or corrupt
/// stream or a bad record is a named `Input` error for the sample. The tool then never reads a
/// file that was not checked to the end (spec 11f: no reliance on tool-side decompression).
fn stage_reads(
    sample: &Sample,
    gzip_ok: bool,
    validate: bool,
    scratch: &Scratch,
    guard: &mut ScratchGuard,
) -> Result<Staged, PipelineError> {
    let named = |source| PipelineError::Input { sample: Some(sample.name.clone()), source };
    let mut staged = sample.clone();
    let mut counts: [Vec<u64>; 2] = [Vec::new(), Vec::new()];
    for (mate, list) in [(1, &mut staged.fastq1), (2, &mut staged.fastq2)] {
        for (lane, p) in list.iter_mut().enumerate() {
            let path = absolute(Path::new(p.as_str()))?;
            let fmt = file_format(&path).map_err(named)?;
            let keep = fmt == Format::Plain || (fmt == Format::Gzip && gzip_ok);
            let new = if keep {
                if validate {
                    counts[mate - 1].push(crate::run::count_pairs(&path, None).map_err(named)?);
                }
                path
            } else {
                let dir = match &guard.0 {
                    Some(d) => d.clone(),
                    None => {
                        let d = scratch.dir(&sample.name)?;
                        guard.0 = Some(d.clone());
                        d
                    }
                };
                let dst = dir.join(format!("lane{}_R{mate}.fq", lane + 1));
                // Counting is free while copying; it also names the sample in a truncation error.
                let n = copy_counting(&path, &dst, &sample.name)?;
                if validate {
                    counts[mate - 1].push(n);
                }
                dst
            };
            *p = utf8(&new)?.to_string();
        }
    }
    Ok(Staged { sample: staged, counts })
}

/// Check the record counts from `stage_reads` (final review C1): no file may be empty, and in
/// each lane mate 1 and mate 2 must hold the same number of records. Returns the total pairs
/// (reads, for single-end).
fn check_read_counts(sample: &Sample, counts: &[Vec<u64>; 2]) -> Result<u64, PipelineError> {
    for (files, c) in [(&sample.fastq1, &counts[0]), (&sample.fastq2, &counts[1])] {
        for (f, &n) in files.iter().zip(c) {
            if n == 0 {
                return Err(PipelineError::NoReads { sample: sample.name.clone(), file: PathBuf::from(f) });
            }
        }
    }
    if !sample.fastq2.is_empty() {
        for (lane, (&r1, &r2)) in counts[0].iter().zip(&counts[1]).enumerate() {
            if r1 != r2 {
                return Err(PipelineError::MateCountMismatch { sample: sample.name.clone(), lane: lane + 1, r1, r2 });
            }
        }
    }
    Ok(counts[0].iter().sum())
}

/// The quantifier's processed count must equal the pairs counted in the inputs (final review
/// C1). `tool` names the quantifier in the error.
fn check_processed(
    sample: &str,
    input: u64,
    processed: u64,
    tool: &'static str,
) -> Result<(), PipelineError> {
    if input == processed {
        Ok(())
    } else {
        Err(PipelineError::PairCountMismatch { sample: sample.to_string(), input, processed, tool })
    }
}

/// A plain copy of `path` in scratch when it is compressed; the path itself when plain.
fn plain_copy(path: &Path, scratch: &Scratch, sub: &str, name: &str) -> Result<(PathBuf, ScratchGuard), PipelineError> {
    let path = absolute(path)?;
    if file_format(&path).map_err(input_err)? == Format::Plain {
        return Ok((path, ScratchGuard(None)));
    }
    let dir = scratch.dir(sub)?;
    let dst = dir.join(name);
    copy_converted(&path, &dst)?;
    Ok((dst, ScratchGuard(Some(dir))))
}

// ---------------------------------------------------------------------------------------------
// Per-sample results

/// Per-sample mapping summary.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleSummary {
    /// kallisto `n_processed`, or HISAT2's `<N> reads; of these:` (pairs for paired data,
    /// reads for single-end).
    pub reads_processed: u64,
    /// kallisto `p_pseudoaligned / 100`, or HISAT2 `overall alignment rate / 100`.
    pub fraction_mapped: f64,
    /// kallisto `p_pseudoaligned / 100`, or featureCounts `Assigned` / input pairs (reads).
    pub fraction_assigned: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SampleResult {
    pub sample: String,
    pub counts: BTreeMap<String, u64>,
    pub summary: SampleSummary,
    /// Warnings for the user (printed by the CLI on stderr).
    pub warnings: Vec<String>,
}

fn genecounts_path(out_dir: &Path, sample: &str) -> PathBuf {
    out_dir.join(format!("{sample}.genecounts.tsv.zst"))
}

/// Ruling R24: remove every output a previous run of this sample may have left, before any
/// tool runs, so a failed rerun cannot leave old counts for `matrix` to use.
pub fn remove_old_outputs(out_dir: &Path, sample: &str) -> Result<(), PipelineError> {
    let sd = out_dir.join(sample);
    let gc = genecounts_path(out_dir, sample);
    let mut gc_partial = gc.as_os_str().to_owned();
    gc_partial.push(".partial");
    let summary = summary_path(out_dir, sample);
    let mut summary_partial = summary.as_os_str().to_owned();
    summary_partial.push(".partial");
    let files = [
        gc,
        PathBuf::from(gc_partial),
        summary,
        PathBuf::from(summary_partial),
        sd.join("abundance.tsv"),
        sd.join("run_info.json"),
        sd.join("quant.sf"),
        sd.join("quant.sf.partial"),
        sd.join("cmd_info.json"),
        sd.join("lib_format_counts.json"),
        sd.join("featurecounts.summary"),
        sd.join("featurecounts.summary.partial"),
        sd.join(format!("{sample}.sam.zst")),
        sd.join(format!("{sample}.sam.zst.partial")),
    ];
    for p in files {
        match fs::remove_file(&p) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(PipelineError::Io { path: p, source: e }),
        }
    }
    // Salmon's per-sample subdirectories (aux_info/meta_info.json, logs,
    // libParams/flenDist.txt and the other files in them).
    for sub in ["aux_info", "logs", "libParams"] {
        let p = sd.join(sub);
        match fs::remove_dir_all(&p) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(PipelineError::Io { path: p, source: e }),
        }
    }
    Ok(())
}

/// File name of the run-wide mapping summary in the output directory.
pub const MAPPING_SUMMARY: &str = "mapping_summary.tsv.zst";

/// Header of `<out>/<sample>.summary.tsv.zst`.
pub const SUMMARY_HEADER: &str = "sample\tmode\tquantifier\treads_processed\tfraction_mapped\tfraction_assigned\tstrandedness\tqrnamap_version";

/// The 7-column header written before the `quantifier` column was added.
const OLD_SUMMARY_HEADER: &str = "sample\tmode\treads_processed\tfraction_mapped\t\
                                  fraction_assigned\tstrandedness\tqrnamap_version";

fn summary_path(out_dir: &Path, sample: &str) -> PathBuf {
    out_dir.join(format!("{sample}.summary.tsv.zst"))
}

/// Remove an old `<out>/mapping_summary.tsv.zst` (and its `.partial`) before a step that
/// removes sample outputs (ruling R24).
pub fn remove_old_mapping_summary(out_dir: &Path) -> Result<(), PipelineError> {
    for name in [MAPPING_SUMMARY.to_string(), format!("{MAPPING_SUMMARY}.partial")] {
        remove_if_present(&out_dir.join(name))?;
    }
    Ok(())
}

fn remove_if_present(p: &Path) -> Result<(), PipelineError> {
    match fs::remove_file(p) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(PipelineError::Io { path: p.to_path_buf(), source: e }),
    }
}

/// Write `<out>/<sample>.summary.tsv.zst`: the header `SUMMARY_HEADER` and one row. Floats use
/// Rust's shortest round-trip formatting. `mode` is `quant` or `align`. `quantifier` is the
/// tool name for quant and `None` (written as `none`) for align. `strandedness` is
/// `unstranded`, `forward` or `reverse`.
pub fn write_sample_summary(
    out_dir: &Path,
    sample: &str,
    mode: Mode,
    quantifier: Option<Quantifier>,
    summary: &SampleSummary,
    strandedness: Strandedness,
) -> Result<(), PipelineError> {
    let mode = match mode {
        Mode::Quant => "quant",
        Mode::Align => "align",
    };
    let quantifier = quantifier.map_or("none", Quantifier::name);
    let strand = match strandedness {
        Strandedness::Unstranded => "unstranded",
        Strandedness::Forward => "forward",
        Strandedness::Reverse => "reverse",
    };
    let text = format!(
        "{SUMMARY_HEADER}\n{sample}\t{mode}\t{quantifier}\t{}\t{}\t{}\t{strand}\t{}\n",
        summary.reads_processed,
        summary.fraction_mapped,
        summary.fraction_assigned,
        env!("CARGO_PKG_VERSION")
    );
    write_output(&summary_path(out_dir, sample), text.as_bytes())
}

/// Read `<out>/<sample>.summary.tsv.zst` and return its data row (8 cells). A missing file is
/// `MissingSummaryFile`; a row for another sample is `SummaryFileSample`.
fn read_sample_summary(out_dir: &Path, sample: &str) -> Result<Vec<String>, PipelineError> {
    let path = summary_path(out_dir, sample);
    if !path.exists() {
        return Err(PipelineError::MissingSummaryFile { sample: sample.into(), path });
    }
    let text = read_table(&path)?;
    let bad = |line: usize, detail: &str| PipelineError::BadTable { path: path.clone(), line, detail: detail.into() };
    let mut lines = text.lines();
    match lines.next() {
        Some(SUMMARY_HEADER) => {}
        Some(OLD_SUMMARY_HEADER) => {
            return Err(bad(
                1,
                "the summary header has no quantifier column: the file was written by an older \
                 qrnamap; run that sample again",
            ))
        }
        _ => return Err(bad(1, "the summary header is not the expected one")),
    }
    let row: Vec<String> = match lines.next() {
        Some(l) => l.split('\t').map(String::from).collect(),
        None => return Err(bad(2, "the summary row is missing")),
    };
    if row.len() != 8 || lines.next().is_some() {
        return Err(bad(2, "expected one row of 8 tab-separated fields"));
    }
    if row[0] != sample {
        return Err(PipelineError::SummaryFileSample { path: path.clone(), expected: sample.into(), found: row[0].clone() });
    }
    Ok(row)
}

/// The `suggested` cell of each row of a `strand-check` table, by sample name.
fn read_strand_suggestions(path: &Path) -> Result<HashMap<String, String>, PipelineError> {
    let text = read_table(path)?;
    let bad = |line: usize, detail: &str| PipelineError::BadTable { path: path.to_path_buf(), line, detail: detail.into() };
    let mut lines = text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty());
    let header: Vec<&str> = lines.next().ok_or_else(|| bad(0, "the table is empty"))?.1.split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).ok_or_else(|| bad(1, &format!("column '{n}' not found")));
    let (si, gi) = (col("sample")?, col("suggested")?);
    let mut out = HashMap::new();
    for (i, l) in lines {
        let f: Vec<&str> = l.split('\t').collect();
        match (f.get(si), f.get(gi)) {
            (Some(s), Some(g)) => out.insert(s.to_string(), g.to_string()),
            _ => return Err(bad(i + 1, "row has too few fields")),
        };
    }
    Ok(out)
}

/// The mapping summary table: one row per sample in `samples` order, read from
/// `<out>/<sample>.summary.tsv.zst`, plus the column `strand_check_suggestion` (the
/// `suggested` cell of the `strand-check` table row for the sample, or empty).
pub fn build_mapping_summary(out_dir: &Path, samples: &[String], strand_check: Option<&Path>) -> Result<String, PipelineError> {
    let suggestions = match strand_check {
        Some(p) => read_strand_suggestions(p)?,
        None => HashMap::new(),
    };
    let mut text = format!("{SUMMARY_HEADER}\tstrand_check_suggestion\n");
    for s in samples {
        let row = read_sample_summary(out_dir, s)?;
        text.push_str(&row.join("\t"));
        text.push('\t');
        text.push_str(suggestions.get(s).map(String::as_str).unwrap_or(""));
        text.push('\n');
    }
    Ok(text)
}

/// Write a per-sample gene count table (`<out>/<sample>.genecounts.tsv.zst`).
///
/// Format: tab-separated text. The first line is the header `gene\t<sample>`; its second
/// field is the sample name, which `read_gene_counts` checks. Each later line is
/// `<gene>\t<count>` (count a non-negative integer), genes in byte order. A gene label with a
/// tab, newline or CR is `BadGeneLabel`.
pub fn write_gene_counts(path: &Path, sample: &str, counts: &BTreeMap<String, u64>) -> Result<(), PipelineError> {
    let mut s = format!("gene\t{sample}\n");
    for (g, c) in counts {
        check_label(g)?;
        s.push_str(&format!("{g}\t{c}\n"));
    }
    write_output(path, s.as_bytes())
}

/// Read a table written by `write_gene_counts` for `sample`. A missing file is
/// `MissingCountFile`; a header that names another sample is `CountFileSample`.
pub fn read_gene_counts(path: &Path, sample: &str) -> Result<BTreeMap<String, u64>, PipelineError> {
    if !path.exists() {
        return Err(PipelineError::MissingCountFile { sample: sample.into(), path: path.to_path_buf() });
    }
    let text = read_table(path)?;
    let bad = |line: usize, detail: &str| PipelineError::BadTable {
        path: path.to_path_buf(),
        line,
        detail: detail.to_string(),
    };
    let mut lines = text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty());
    match lines.next().map(|(_, h)| h.split('\t').collect::<Vec<_>>()) {
        Some(h) if h.len() == 2 && h[0] == "gene" => {
            if h[1] != sample {
                return Err(PipelineError::CountFileSample {
                    path: path.to_path_buf(),
                    expected: sample.into(),
                    found: h[1].into(),
                });
            }
        }
        _ => return Err(bad(1, "header 'gene\\t<sample>' expected")),
    }
    let mut out = BTreeMap::new();
    for (i, l) in lines {
        let f: Vec<&str> = l.split('\t').collect();
        let n = match f.as_slice() {
            [g, n] if !g.is_empty() => n.parse::<u64>().map_err(|_| bad(i + 1, "count is not a non-negative integer"))?,
            _ => return Err(bad(i + 1, "expected two tab-separated fields")),
        };
        if out.insert(f[0].to_string(), n).is_some() {
            return Err(bad(i + 1, "gene occurs twice"));
        }
    }
    Ok(out)
}

/// Find `"field"` in a JSON object text and parse the number after its colon. Enough for the
/// flat `run_info.json` that kallisto writes; not a general JSON parser.
fn json_number(text: &str, field: &str) -> Option<f64> {
    let key = format!("\"{field}\"");
    let rest = &text[text.find(&key)? + key.len()..];
    let rest = rest.trim_start().strip_prefix(':')?.trim_start();
    let end = rest
        .find(|c: char| !(c.is_ascii_digit() || matches!(c, '-' | '+' | '.' | 'e' | 'E')))
        .unwrap_or(rest.len());
    rest[..end].parse().ok().filter(|v: &f64| v.is_finite())
}

/// The raw value text after `"field":` in a JSON text, up to the next `,`, `}`, `]` or
/// white space. `Ok(None)` when the key is absent. A key that occurs more than once (for
/// example again in a nested object) is a `BadTable` error, so a nested value is never read
/// by mistake. Not a general JSON parser.
fn json_value_once<'a>(
    text: &'a str,
    field: &str,
    path: &Path,
) -> Result<Option<&'a str>, PipelineError> {
    let key = format!("\"{field}\"");
    let mut found = None;
    for (i, _) in text.match_indices(&key) {
        let Some(rest) = text[i + key.len()..].trim_start().strip_prefix(':') else {
            continue; // the text occurs as a value, not as a key
        };
        if found.is_some() {
            let detail = format!("field '{field}' occurs twice");
            return Err(PipelineError::BadTable { path: path.to_path_buf(), line: 0, detail });
        }
        let rest = rest.trim_start();
        let end = rest.find(|c: char| matches!(c, ',' | '}' | ']') || c.is_whitespace());
        found = Some(&rest[..end.unwrap_or(rest.len())]);
    }
    Ok(found)
}

/// `(num_processed, num_mapped, percent_mapped)` from salmon's `aux_info/meta_info.json`.
/// Each key must occur once. The counts must be non-negative integers; the percentage must
/// be a finite number from 0 to 100.
fn parse_salmon_meta(text: &str, path: &Path) -> Result<(u64, u64, f64), PipelineError> {
    let raw = |f: &str| {
        let missing = || PipelineError::MissingField { path: path.to_path_buf(), field: f.into() };
        json_value_once(text, f, path)?.ok_or_else(missing)
    };
    let bad = |f: &str, what: &str| PipelineError::BadTable {
        path: path.to_path_buf(),
        line: 0,
        detail: format!("field '{f}' is not {what}"),
    };
    let int = |f: &str| raw(f)?.parse::<u64>().map_err(|_| bad(f, "a non-negative integer"));
    let processed = int("num_processed")?;
    let mapped = int("num_mapped")?;
    let pct = raw("percent_mapped")?
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && (0.0..=100.0).contains(v))
        .ok_or_else(|| bad("percent_mapped", "a number from 0 to 100"))?;
    Ok((processed, mapped, pct))
}

// ---------------------------------------------------------------------------------------------
// kallisto and salmon

/// Everything `quant_sample` needs.
#[derive(Debug, Clone)]
pub struct QuantContext {
    /// The index of `config.quantifier`: the kallisto index file (`transcripts.kidx`) or the
    /// salmon index directory (`salmon_index`).
    pub index: PathBuf,
    pub tx2gene: HashMap<String, String>,
    pub out_dir: PathBuf,
    pub scratch: Scratch,
    pub config: RunConfig,
}

/// Read a `transcript\tgene` table (no header).
pub fn read_tx2gene(path: &Path) -> Result<HashMap<String, String>, PipelineError> {
    let text = read_table(path)?;
    let mut m = HashMap::new();
    for (i, l) in text.lines().enumerate() {
        if l.trim().is_empty() {
            continue;
        }
        match l.trim_end_matches('\r').split('\t').collect::<Vec<_>>().as_slice() {
            [t, g] if !t.is_empty() && !g.is_empty() => {
                if m.insert(t.to_string(), g.to_string()).is_some() {
                    return Err(PipelineError::BadTable { path: path.to_path_buf(), line: i + 1, detail: format!("transcript '{t}' occurs twice") });
                }
            }
            _ => {
                return Err(PipelineError::BadTable { path: path.to_path_buf(), line: i + 1, detail: "expected 'transcript<TAB>gene'".into() })
            }
        }
    }
    Ok(m)
}

/// Build the index of `config.quantifier` once: kallisto in `<out_dir>/transcripts.kidx`,
/// salmon in the directory `<out_dir>/salmon_index` (an old one is removed first). Only that
/// one index is built. A compressed transcript FASTA is decompressed into scratch first and
/// removed after the index is built. With `genes`, `<out_dir>/genes.tsv.zst` is written from
/// that GFF3 after the index.
pub fn prepare_quant(
    transcripts: &Path,
    tx2gene: &Path,
    out_dir: &Path,
    scratch: Scratch,
    config: RunConfig,
    genes: Option<GenesSource<'_>>,
) -> Result<QuantContext, PipelineError> {
    let out_dir = absolute(out_dir)?;
    fs::create_dir_all(&out_dir).map_err(io_err(&out_dir))?;
    let tx2gene = read_tx2gene(tx2gene)?;
    // Read and check the GFF3 before the index build; write the table after it, so a failed
    // prepare leaves no table.
    let model = genes.map(|g| load_model(g.gff)).transpose()?;
    if let Some(m) = &model {
        check_tx2gene_in_model(&tx2gene, m)?;
    }
    let exec_err = |source| PipelineError::Exec { sample: None, source };
    let index = match config.quantifier {
        Quantifier::Kallisto => {
            let (fa, _guard) =
                plain_copy(transcripts, &scratch, "kallisto-index", "transcripts.fa")?;
            let index = out_dir.join("transcripts.kidx");
            let call = kallisto_index(utf8(&fa)?, utf8(&index)?);
            run_call(&call, &out_dir).map_err(exec_err)?;
            index
        }
        Quantifier::Salmon => {
            let index = out_dir.join("salmon_index");
            let (fa, _guard) = plain_copy(transcripts, &scratch, "salmon-index", "transcripts.fa")?;
            let call = salmon_index(utf8(&fa)?, utf8(&index)?, &config)
                .map_err(|source| PipelineError::Config { sample: "(index)".into(), source })?;
            // R24: the inputs are checked above; remove a stale index, then build.
            match fs::remove_dir_all(&index) {
                Ok(()) => {}
                Err(e) if e.kind() == io::ErrorKind::NotFound => {}
                Err(e) => return Err(PipelineError::Io { path: index, source: e }),
            }
            run_call(&call, &out_dir).map_err(exec_err)?;
            index
        }
    };
    if let (Some(m), Some(g)) = (&model, genes) {
        write_genes_table(&out_dir.join(GENES_TABLE), m, g.mito_contigs)?;
    }
    Ok(QuantContext { index, tx2gene, out_dir, scratch, config })
}

/// The quant tool call of `cfg.quantifier` for `sample`.
fn quant_call(
    sample: &Sample,
    index: &str,
    out: &str,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    match cfg.quantifier {
        Quantifier::Kallisto => kallisto_quant(sample, index, out, cfg),
        Quantifier::Salmon => salmon_quant(sample, index, out, cfg),
    }
}

/// Quantify one sample with `ctx.config.quantifier`. kallisto writes
/// `<out>/<sample>/abundance.tsv` and `run_info.json`; salmon writes `<out>/<sample>/quant.sf`
/// and `aux_info/meta_info.json` (the tools' names). Both write
/// `<out>/<sample>.genecounts.tsv.zst` and `<out>/<sample>.summary.tsv.zst`.
pub fn quant_sample(sample: &Sample, ctx: &QuantContext) -> Result<SampleResult, PipelineError> {
    let name = &sample.name;
    remove_old_outputs(&ctx.out_dir, name)?;
    let cfg_err = |source| PipelineError::Config { sample: name.clone(), source };
    quant_call(sample, "idx", "out", &ctx.config).map_err(cfg_err)?;
    let sample_dir = ctx.out_dir.join(name);
    fs::create_dir_all(&sample_dir).map_err(io_err(&sample_dir))?;
    let mut guard = ScratchGuard(None);
    // Both tools read gzip themselves (salmon 2.8.0: checked with a real run in the tests).
    let staged = stage_reads(sample, true, true, &ctx.scratch, &mut guard)?;
    let input_pairs = check_read_counts(sample, &staged.counts)?;
    let staged = staged.sample;
    let call = quant_call(&staged, utf8(&ctx.index)?, utf8(&sample_dir)?, &ctx.config)
        .map_err(cfg_err)?;
    let quantifier = ctx.config.quantifier;
    let (counts, n, p) = match quantifier {
        Quantifier::Kallisto => kallisto_run(name, &call, &sample_dir, ctx)?,
        Quantifier::Salmon => salmon_run(name, &call, &sample_dir, ctx)?,
    };
    drop(guard);
    check_processed(name, input_pairs, n, quantifier.name())?;
    write_gene_counts(&genecounts_path(&ctx.out_dir, name), name, &counts)?;
    let warnings = zero_warning(name, &counts).into_iter().collect();
    let summary = SampleSummary { reads_processed: n, fraction_mapped: p, fraction_assigned: p };
    let strand = ctx.config.strandedness;
    write_sample_summary(&ctx.out_dir, name, Mode::Quant, Some(quantifier), &summary, strand)?;
    Ok(SampleResult { sample: name.clone(), counts, summary, warnings })
}

/// Gene counts, processed pairs and mapped fraction (0 to 1).
type QuantOutput = (BTreeMap<String, u64>, u64, f64);

/// Run kallisto quant and read `abundance.tsv` and `run_info.json`.
fn kallisto_run(
    name: &str,
    call: &ToolCall,
    sample_dir: &Path,
    ctx: &QuantContext,
) -> Result<QuantOutput, PipelineError> {
    run_call(call, sample_dir).map_err(|source| match source {
        // kallisto 0.52.0 exits with status 1 when no read pseudoaligns (observed in the tests).
        ExecError::Failed { ref stderr_tail, .. } if stderr_tail.contains("no reads pseudoaligned") => {
            PipelineError::NoReadsPseudoaligned {
                sample: name.to_string(),
                tool: "kallisto",
                strandedness: ctx.config.strandedness,
                stderr_tail: stderr_tail.clone(),
            }
        }
        source => PipelineError::Exec { sample: Some(name.to_string()), source },
    })?;
    let ab_path = sample_dir.join("abundance.tsv");
    let counts = gene_counts_from_abundance(&read_table(&ab_path)?, &ctx.tx2gene)
        .map_err(|source| PipelineError::Count { path: ab_path.clone(), source })?;
    let ri_path = sample_dir.join("run_info.json");
    let ri = read_table(&ri_path)?;
    let field = |f: &str| {
        json_number(&ri, f).ok_or_else(|| PipelineError::MissingField { path: ri_path.clone(), field: f.into() })
    };
    let n = field("n_processed")?;
    let p = field("p_pseudoaligned")? / 100.0;
    Ok((counts, n as u64, p))
}

/// True when `meta_info.json` at `path` exists, parses and has `num_mapped` 0.
fn salmon_meta_zero_mapped(path: &Path) -> bool {
    fs::read_to_string(path)
        .ok()
        .and_then(|t| parse_salmon_meta(&t, path).ok())
        .is_some_and(|(_, mapped, _)| mapped == 0)
}

/// Run salmon quant and read `quant.sf` and `aux_info/meta_info.json`. Zero mapped fragments
/// is `NoReadsPseudoaligned`. salmon MIGRATION.md says salmon 2.x writes its output and then
/// exits non-zero in that case, so a failed run is checked through `meta_info.json` first.
/// salmon 2.8.0 exited 0 in a manual zero-mapped check, so a successful run is checked too.
/// `remove_old_outputs` removed any old `meta_info.json`, so the file is from this run.
fn salmon_run(
    name: &str,
    call: &ToolCall,
    sample_dir: &Path,
    ctx: &QuantContext,
) -> Result<QuantOutput, PipelineError> {
    let meta_path = sample_dir.join("aux_info").join("meta_info.json");
    let no_reads = |stderr_tail: &str| PipelineError::NoReadsPseudoaligned {
        sample: name.to_string(),
        tool: "salmon",
        strandedness: ctx.config.strandedness,
        stderr_tail: stderr_tail.to_string(),
    };
    let stderr_tail = match run_call_stderr(call, sample_dir) {
        Ok(t) => t,
        Err(ExecError::Failed { program, status, stderr_tail, stdout_tail }) => {
            if salmon_meta_zero_mapped(&meta_path) {
                return Err(no_reads(&stderr_tail));
            }
            let source = ExecError::Failed { program, status, stderr_tail, stdout_tail };
            return Err(PipelineError::Exec { sample: Some(name.to_string()), source });
        }
        Err(source) => return Err(PipelineError::Exec { sample: Some(name.to_string()), source }),
    };
    let (n, mapped, pct) = parse_salmon_meta(&read_table(&meta_path)?, &meta_path)?;
    if mapped == 0 {
        return Err(no_reads(&stderr_tail));
    }
    let sf_path = sample_dir.join("quant.sf");
    let counts = gene_counts_from_salmon(&read_table(&sf_path)?, &ctx.tx2gene)
        .map_err(|source| PipelineError::Count { path: sf_path.clone(), source })?;
    Ok((counts, n, pct / 100.0))
}

fn zero_warning(sample: &str, counts: &BTreeMap<String, u64>) -> Option<String> {
    if counts.values().all(|&c| c == 0) {
        Some(format!("sample '{sample}': all-zero sample (no counts in any gene)"))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------------------------
// HISAT2 + featureCounts

/// Everything `align_sample` and `strand_check_sample` need.
#[derive(Debug, Clone)]
pub struct AlignContext {
    pub index_prefix: PathBuf,
    pub splice_sites: PathBuf,
    pub saf: PathBuf,
    pub out_dir: PathBuf,
    pub scratch: Scratch,
    pub config: RunConfig,
    /// Keep a zstd copy of the SAM as `<out>/<sample>/<sample>.sam.zst`.
    pub keep_sam: bool,
}

/// Build splice sites (`<out>/splice_sites.txt`), SAF (`<out>/genes.saf`) and the HISAT2
/// index (`<out>/index/genome.N.ht2`) once. These two tables stay plain because the tools read
/// them; they hold one line per intron or exon. A compressed reference is decompressed into
/// scratch for `hisat2-build` (reading gzip there is UNVERIFIED) and removed afterwards.
/// Also writes `<out>/genes.tsv.zst`, after the index is built.
/// The transcripts are built too, so a gene on a contig missing from the FASTA is an error.
pub fn prepare_align(
    reference: &Path,
    gff: &Path,
    out_dir: &Path,
    scratch: Scratch,
    config: RunConfig,
    keep_sam: bool,
    mito_contigs: &[String],
) -> Result<AlignContext, PipelineError> {
    prepare_align_files(reference, gff, out_dir, scratch, config, keep_sam, Some(mito_contigs))
}

/// `prepare_align`; the gene table is written only when `genes_table` is `Some(mito_contigs)`.
fn prepare_align_files(
    reference: &Path,
    gff: &Path,
    out_dir: &Path,
    scratch: Scratch,
    config: RunConfig,
    keep_sam: bool,
    genes_table: Option<&[String]>,
) -> Result<AlignContext, PipelineError> {
    let out_dir = absolute(out_dir)?;
    fs::create_dir_all(&out_dir).map_err(io_err(&out_dir))?;
    let (genome, model) = load_reference(reference, gff)?;
    build_transcripts(&genome, &model).map_err(PipelineError::Build)?;
    drop(genome);
    let splice_path = out_dir.join("splice_sites.txt");
    write_output(&splice_path, splice_sites(&model).as_bytes())?;
    let saf_path = out_dir.join("genes.saf");
    write_output(&saf_path, saf(&model).as_bytes())?;
    let index_dir = out_dir.join("index");
    fs::create_dir_all(&index_dir).map_err(io_err(&index_dir))?;
    let prefix = index_dir.join("genome");
    let (fa, _guard) = plain_copy(reference, &scratch, "hisat2-build", "genome.fa")?;
    let call = hisat2_build(utf8(&fa)?, utf8(&prefix)?, &config)
        .map_err(|source| PipelineError::Config { sample: "(index)".into(), source })?;
    run_call(&call, &index_dir).map_err(|source| PipelineError::Exec { sample: None, source })?;
    if let Some(mito) = genes_table {
        write_genes_table(&out_dir.join(GENES_TABLE), &model, mito)?;
    }
    Ok(AlignContext {
        index_prefix: prefix,
        splice_sites: splice_path,
        saf: saf_path,
        out_dir,
        scratch,
        config,
        keep_sam,
    })
}

/// `prepare_align` for `strand-check`: the same reference files and index, but no
/// `<out>/genes.tsv.zst` (R35-1: a finished run directory keeps its gene table).
pub fn prepare_strand_check(
    reference: &Path,
    gff: &Path,
    out_dir: &Path,
    scratch: Scratch,
    config: RunConfig,
) -> Result<AlignContext, PipelineError> {
    prepare_align_files(reference, gff, out_dir, scratch, config, false, None)
}

fn summary_field(m: &BTreeMap<String, u64>, path: &Path, f: &str) -> Result<u64, PipelineError> {
    m.get(f).copied().ok_or_else(|| PipelineError::MissingField { path: path.to_path_buf(), field: f.into() })
}

fn read_summary(path: &Path) -> Result<BTreeMap<String, u64>, PipelineError> {
    parse_featurecounts_summary(&read_table(path)?).map_err(|source| PipelineError::Count { path: path.to_path_buf(), source })
}

/// Input pairs (reads, for single-end) and the overall alignment rate (0 to 1) from HISAT2's
/// summary on stderr: the lines `<N> reads; of these:` and `<R>% overall alignment rate`.
/// The sum of the featureCounts summary is not used for the total: it counts each alignment
/// of a multi-mapping fragment (`Unassigned_MultiMapping`), so it can exceed the input count.
fn hisat2_summary(stderr: &str) -> Option<(u64, f64)> {
    let total = stderr
        .lines()
        .find_map(|l| l.trim().strip_suffix(" reads; of these:"))?
        .trim()
        .parse()
        .ok()?;
    let rate: f64 = stderr
        .lines()
        .find_map(|l| l.trim().strip_suffix("% overall alignment rate"))?
        .trim()
        .parse()
        .ok()?;
    Some((total, rate / 100.0))
}

/// HISAT2 then featureCounts for one sample. The SAM goes to scratch and is deleted after
/// counting (a zstd copy is kept with `keep_sam`). Writes `<out>/<sample>/featurecounts.summary`
/// and `<out>/<sample>.genecounts.tsv.zst`.
pub fn align_sample(sample: &Sample, ctx: &AlignContext) -> Result<SampleResult, PipelineError> {
    let name = &sample.name;
    remove_old_outputs(&ctx.out_dir, name)?;
    check_samples(std::slice::from_ref(sample), &ctx.config, Mode::Align)?;
    let exec_err = |source| PipelineError::Exec { sample: Some(name.clone()), source };
    let cfg_err = |source| PipelineError::Config { sample: name.clone(), source };
    let paired = !sample.fastq2.is_empty();
    let sample_dir = ctx.out_dir.join(name);
    fs::create_dir_all(&sample_dir).map_err(io_err(&sample_dir))?;
    let sdir = ctx.scratch.dir(name)?;
    let mut guard = ScratchGuard(Some(sdir.clone()));
    let staged = stage_reads(sample, false, false, &ctx.scratch, &mut guard)?.sample;

    let sam = sdir.join(format!("{name}.sam"));
    let call = hisat2_align(&staged, utf8(&ctx.index_prefix)?, utf8(&ctx.splice_sites)?, utf8(&sam)?, &ctx.config).map_err(cfg_err)?;
    let hisat2_stderr = run_call_stderr(&call, &sdir).map_err(exec_err)?;
    let (total, mapped_rate) = hisat2_summary(&hisat2_stderr).ok_or_else(|| PipelineError::MissingField {
        path: PathBuf::from("hisat2-align-s stderr"),
        field: "'N reads; of these:' and 'overall alignment rate'".into(),
    })?;
    let fc_out = sdir.join("featurecounts.tsv");
    let call = feature_counts(utf8(&sam)?, utf8(&ctx.saf)?, utf8(&fc_out)?, paired, &ctx.config).map_err(cfg_err)?;
    run_call(&call, &sdir).map_err(exec_err)?;

    let counts = parse_featurecounts(&read_table(&fc_out)?)
        .map_err(|source| PipelineError::Count { path: fc_out.clone(), source })?;
    let sum_path = sdir.join("featurecounts.tsv.summary");
    let summary = read_summary(&sum_path)?;
    copy_converted(&sum_path, &sample_dir.join("featurecounts.summary"))?;
    let assigned = summary_field(&summary, &sum_path, "Assigned")?;
    let frac = |n: u64| if total == 0 { 0.0 } else { n as f64 / total as f64 };

    let mut warnings = Vec::new();
    if counts.is_empty() {
        warnings.push(format!("sample '{name}': the featureCounts table has no gene rows (all-zero sample)"));
    } else if let Some(w) = zero_warning(name, &counts) {
        warnings.push(w);
    }
    if ctx.keep_sam {
        copy_converted(&sam, &sample_dir.join(format!("{name}.sam.zst")))?;
    }
    write_gene_counts(&genecounts_path(&ctx.out_dir, name), name, &counts)?;
    drop(guard);
    let summary = SampleSummary {
        reads_processed: total,
        fraction_mapped: mapped_rate,
        fraction_assigned: frac(assigned),
    };
    write_sample_summary(&ctx.out_dir, name, Mode::Align, None, &summary, ctx.config.strandedness)?;
    Ok(SampleResult { sample: name.clone(), counts, summary, warnings })
}

// ---------------------------------------------------------------------------------------------
// Strand check

/// Copy the first `n` FASTQ records (4 lines each, quality kept) from `inputs`, in order, to
/// a plain file. Returns the number copied.
fn copy_first_records(inputs: &[String], out: &Path, n: u64) -> Result<u64, PipelineError> {
    let mut w = BufWriter::new(File::create(out).map_err(io_err(out))?);
    let mut copied = 0u64;
    let mut line = String::new();
    'files: for p in inputs {
        let path = PathBuf::from(p);
        let mut r = open_input(&path).map_err(input_err)?;
        let mut lineno = 0usize;
        while copied < n {
            for k in 0..4 {
                line.clear();
                let got = r.read_line(&mut line).map_err(io_err(&path))?;
                lineno += 1;
                if got == 0 {
                    if k == 0 {
                        continue 'files;
                    }
                    return Err(PipelineError::BadTable { path, line: lineno, detail: "truncated FASTQ record".into() });
                }
                if k == 0 && !line.starts_with('@') {
                    return Err(PipelineError::BadTable { path, line: lineno, detail: "FASTQ header does not start with @".into() });
                }
                if !line.ends_with('\n') {
                    line.push('\n');
                }
                w.write_all(line.as_bytes()).map_err(io_err(out))?;
            }
            copied += 1;
        }
        break;
    }
    w.flush().map_err(io_err(out))?;
    Ok(copied)
}

/// Strand check for one sample: map the first `pairs` pairs (reads, for single-end) with
/// HISAT2 without strandedness, run featureCounts with `-s 1`, `-s 2` and `-s 0`, and pass
/// (`-s 1` Assigned, `-s 2` Assigned, `-s 0` Unassigned_Ambiguity) to `assigned_fraction`.
/// See `qrnamap_core::strand` for the definitions. Uses scratch only.
pub fn strand_check_sample(sample: &Sample, ctx: &AlignContext, pairs: u64) -> Result<StrandReport, PipelineError> {
    let name = &sample.name;
    check_samples(std::slice::from_ref(sample), &ctx.config, Mode::Align)?;
    let exec_err = |source| PipelineError::Exec { sample: Some(name.clone()), source };
    let cfg_err = |source| PipelineError::Config { sample: name.clone(), source };
    let paired = !sample.fastq2.is_empty();
    let sdir = ctx.scratch.dir(&format!("{name}.strand-check"))?;
    let _guard = ScratchGuard(Some(sdir.clone()));

    let abs_list = |v: &[String]| -> Result<Vec<String>, PipelineError> {
        v.iter().map(|p| absolute(Path::new(p)).and_then(|a| utf8(&a).map(str::to_string))).collect()
    };
    let mut sub = sample.clone();
    let r1 = sdir.join("sub_R1.fq");
    let n1 = copy_first_records(&abs_list(&sample.fastq1)?, &r1, pairs)?;
    sub.fastq1 = vec![utf8(&r1)?.to_string()];
    if paired {
        let r2 = sdir.join("sub_R2.fq");
        let n2 = copy_first_records(&abs_list(&sample.fastq2)?, &r2, pairs)?;
        if n1 != n2 {
            return Err(PipelineError::BadTable {
                path: PathBuf::from(&sample.fastq2[0]),
                line: 0,
                detail: format!("sample '{name}': mate files hold {n1} and {n2} records"),
            });
        }
        sub.fastq2 = vec![utf8(&r2)?.to_string()];
    }
    let mut cfg = ctx.config.clone();
    cfg.strandedness = Strandedness::Unstranded;
    let sam = sdir.join("sub.sam");
    let call = hisat2_align(&sub, utf8(&ctx.index_prefix)?, utf8(&ctx.splice_sites)?, utf8(&sam)?, &cfg).map_err(cfg_err)?;
    run_call(&call, &sdir).map_err(exec_err)?;
    let mut run = |s: Strandedness, tag: &str, field: &str| -> Result<u64, PipelineError> {
        cfg.strandedness = s;
        let out = sdir.join(format!("fc_{tag}.tsv"));
        let call = feature_counts(utf8(&sam)?, utf8(&ctx.saf)?, utf8(&out)?, paired, &cfg).map_err(cfg_err)?;
        run_call(&call, &sdir).map_err(exec_err)?;
        let sp = sdir.join(format!("fc_{tag}.tsv.summary"));
        summary_field(&read_summary(&sp)?, &sp, field)
    };
    let sense = run(Strandedness::Forward, "s1", "Assigned")?;
    let antisense = run(Strandedness::Reverse, "s2", "Assigned")?;
    let ambiguous = run(Strandedness::Unstranded, "s0", "Unassigned_Ambiguity")?;
    Ok(assigned_fraction((sense, antisense, ambiguous)))
}

// ---------------------------------------------------------------------------------------------
// Matrix

/// Source of the gene list for `matrix --genes-from`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenesFrom {
    /// `transcript\tgene` table: second column.
    Tx2Gene,
    /// SAF with header `GeneID...`: first column.
    Saf,
}

/// Unique gene ids in order of first appearance.
pub fn read_gene_list(kind: GenesFrom, path: &Path) -> Result<Vec<String>, PipelineError> {
    gene_list_from_text(kind, &read_table(path)?, path)
}

/// The gene list (SAF order) built from a GFF3 file, as `prepare_align` writes it to
/// `genes.saf`. Labels with a tab, newline or CR are `BadGeneLabel`.
pub fn gene_list_from_gff(gff: &Path) -> Result<Vec<String>, PipelineError> {
    let model = parse_gff3(open_input(gff).map_err(input_err)?)
        .map_err(|e| PipelineError::Gff { path: gff.to_path_buf(), detail: e.to_string() })?;
    for g in &model.genes {
        check_label(&g.label)?;
    }
    gene_list_from_text(GenesFrom::Saf, &saf(&model), gff)
}

fn gene_list_from_text(kind: GenesFrom, text: &str, path: &Path) -> Result<Vec<String>, PipelineError> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for (i, l) in text.lines().enumerate() {
        let l = l.trim_end_matches('\r');
        if l.trim().is_empty() || (kind == GenesFrom::Saf && l.starts_with("GeneID\t")) {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        let g = match (kind, f.as_slice()) {
            (GenesFrom::Tx2Gene, [_, g]) => *g,
            (GenesFrom::Saf, [g, _, _, _, _]) => *g,
            _ => return Err(PipelineError::BadTable { path: path.to_path_buf(), line: i + 1, detail: "unexpected number of columns".into() }),
        };
        if g.is_empty() {
            return Err(PipelineError::BadTable { path: path.to_path_buf(), line: i + 1, detail: "empty gene id".into() });
        }
        if seen.insert(g.to_string()) {
            out.push(g.to_string());
        }
    }
    Ok(out)
}

/// Build the gene-by-sample matrix and the warnings for the user.
///
/// Rows: `genes` when given, else the sorted union of the sample genes. Any gene label with a
/// tab, newline or CR is `BadGeneLabel`. Warnings: genes in a sample but not in `genes`
/// (dropped from the matrix), and each sample whose counts are all zero or empty.
pub fn build_matrix(
    per_sample: &[(String, BTreeMap<String, u64>)],
    genes: Option<&[String]>,
) -> Result<(CountMatrix, Vec<String>), PipelineError> {
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
                let more = if extra.len() > MAX_LISTED { ", ..." } else { "" };
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

/// `qrnamap matrix`: read `<out>/<sample>.genecounts.tsv.zst` for each sample in sheet order
/// and write `<out>/counts.tsv.zst` and `<out>/mapping_summary.tsv.zst` (from the per-sample
/// summaries and `<out>/strand_check.tsv` when it exists). Returns the matrix path and the
/// warnings.
pub fn run_matrix(
    out_dir: &Path,
    samples: &[Sample],
    genes_from: Option<(GenesFrom, &Path)>,
) -> Result<(PathBuf, Vec<String>), PipelineError> {
    let genes = match genes_from {
        Some((kind, p)) => Some(read_gene_list(kind, p)?),
        None => None,
    };
    run_matrix_with_genes(out_dir, samples, genes.as_deref())
}

/// `run_matrix` with the row list given directly (`None`: the sorted union of sample genes).
pub fn run_matrix_with_genes(
    out_dir: &Path,
    samples: &[Sample],
    genes: Option<&[String]>,
) -> Result<(PathBuf, Vec<String>), PipelineError> {
    // A failing matrix must not leave an old mapping summary or an old matrix (R24, R35-2).
    // The DESeq2 outputs are removed by `de` itself (R26).
    remove_old_mapping_summary(out_dir)?;
    remove_if_present(&out_dir.join("counts.tsv.zst"))?;
    remove_if_present(&out_dir.join("counts.tsv.zst.partial"))?;
    let mut per_sample = Vec::with_capacity(samples.len());
    for s in samples {
        per_sample.push((s.name.clone(), read_gene_counts(&genecounts_path(out_dir, &s.name), &s.name)?));
    }
    let (matrix, warnings) = build_matrix(&per_sample, genes)?;
    // Check the summaries before any matrix file is written; write the summary table last.
    let strand_check = out_dir.join("strand_check.tsv");
    let strand_check = strand_check.is_file().then_some(strand_check);
    let names: Vec<String> = samples.iter().map(|s| s.name.clone()).collect();
    let summary = build_mapping_summary(out_dir, &names, strand_check.as_deref())?;
    let path = out_dir.join("counts.tsv.zst");
    write_output(&path, matrix.to_tsv().as_bytes())?;
    write_output(&out_dir.join(MAPPING_SUMMARY), summary.as_bytes())?;
    Ok((path, warnings))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_number_reads_kallisto_fields() {
        let t = "{\n\t\"n_targets\": 6,\n\t\"n_processed\": 10000,\n\t\"p_pseudoaligned\": 99.8,\n\t\"x\": -1e2\n}";
        assert_eq!(json_number(t, "n_processed"), Some(10000.0));
        assert_eq!(json_number(t, "p_pseudoaligned"), Some(99.8));
        assert_eq!(json_number(t, "x"), Some(-100.0));
        assert_eq!(json_number(t, "missing"), None);
        assert_eq!(json_number("{\"a\": \"text\"}", "a"), None);
    }

    #[test]
    fn processed_count_must_match_input_pairs() {
        check_processed("S1", 10, 10, "kallisto").unwrap();
        let e = check_processed("S1", 10, 4, "kallisto").unwrap_err();
        assert!(
            matches!(
                e,
                PipelineError::PairCountMismatch {
                    ref sample, input: 10, processed: 4, tool: "kallisto"
                } if sample == "S1"
            ),
            "{e:?}"
        );
        assert!(e.to_string().contains("'S1'") && e.to_string().contains("kallisto"), "{e}");
    }

    #[test]
    fn salmon_pair_count_mismatch_names_tool() {
        let e = check_processed("s", 10, 9, "salmon").unwrap_err();
        let t = e.to_string();
        assert!(t.contains("salmon processed 9") && !t.contains("kallisto"), "{t}");
    }

    #[test]
    fn parse_salmon_meta_ok_and_errors() {
        let p = Path::new("meta_info.json");
        let good = "{\n  \"salmon_version\": \"2.8.0\",\n  \"num_processed\": 10000,\n  \
                    \"num_mapped\": 9950,\n  \"percent_mapped\": 99.5,\n  \
                    \"frag_dist\": [1, 2],\n  \"x\": {\"y\": 1}\n}";
        assert_eq!(parse_salmon_meta(good, p).unwrap(), (10000, 9950, 99.5));
        // Missing key.
        let e = parse_salmon_meta(&good.replace("\"num_mapped\"", "\"other\""), p).unwrap_err();
        let missing =
            matches!(&e, PipelineError::MissingField { field, .. } if field == "num_mapped");
        assert!(missing, "{e:?}");
        // The key occurs twice: a nested object repeats it.
        let dup = good.replace("\"x\": {", "\"x\": {\"num_processed\": 3, ");
        let e = parse_salmon_meta(&dup, p).unwrap_err();
        assert!(e.to_string().contains("num_processed") && e.to_string().contains("twice"), "{e}");
        // Non-numeric value, a negative count and a fractional count.
        for bad in ["\"many\"", "-1", "1.5"] {
            let text =
                good.replace("\"num_processed\": 10000", &format!("\"num_processed\": {bad}"));
            let e = parse_salmon_meta(&text, p).unwrap_err();
            assert!(matches!(&e, PipelineError::BadTable { .. }), "{bad}: {e:?}");
            assert!(e.to_string().contains("num_processed"), "{e}");
        }
        let e = parse_salmon_meta(&good.replace("99.5", "\"high\""), p).unwrap_err();
        assert!(e.to_string().contains("percent_mapped"), "{e}");
    }

    #[test]
    fn salmon_zero_mapped_check_reads_meta_file() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("meta_info.json");
        assert!(!salmon_meta_zero_mapped(&p), "missing file");
        let meta = |m: u64| {
            format!("{{\"num_processed\": 5, \"num_mapped\": {m}, \"percent_mapped\": 0.0}}")
        };
        fs::write(&p, meta(0)).unwrap();
        assert!(salmon_meta_zero_mapped(&p));
        fs::write(&p, meta(3)).unwrap();
        assert!(!salmon_meta_zero_mapped(&p));
        fs::write(&p, "{\"num_mapped\": 0}").unwrap();
        assert!(!salmon_meta_zero_mapped(&p), "a file that does not parse");
    }

    #[test]
    fn output_compression_follows_extension() {
        let d = tempfile::tempdir().unwrap();
        for (name, magic) in [("a.zst", &[0x28u8, 0xb5][..]), ("a.gz", &[0x1f, 0x8b][..]), ("a.tsv", &b"hi"[..])] {
            let p = d.path().join(name);
            write_output(&p, b"hi\n").unwrap();
            assert_eq!(&fs::read(&p).unwrap()[..2], magic, "{name}");
            assert_eq!(read_table(&p).unwrap(), "hi\n");
            assert!(!d.path().join(format!("{name}.partial")).exists());
        }
    }
}
