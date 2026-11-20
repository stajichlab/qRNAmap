//! `qrnamap run`: all steps in one command, or the SLURM scripts for them.
//!
//! Local order (`--step all`): parse the sheet and run every tool-call builder (a single-end
//! sample without a fragment length fails here), check the design and contrast, check reads
//! per sample (warnings only), remove the old outputs of ALL samples, build the reference files
//! and the index, write the completion marker `<out>/.prepared`, run each sample, build the
//! matrix, then DESeq2 when a contrast is given. `--step prepare` stops after the marker;
//! `--step final` does the matrix and DESeq2 only. The SLURM scripts call those two steps and
//! `qrnamap quant|align --task-id` in between.
//!
//! Files in the run directory: `transcripts.fa.zst`, `tx2gene.tsv.zst` and `transcripts.kidx`
//! (quant); `splice_sites.txt`, `genes.saf` and `index/` (align); `<sample>/` and
//! `<sample>.genecounts.tsv.zst`; `counts.tsv.zst`; `deseq2.*` (with a contrast); `slurm/`.

use std::fmt;
use std::fs;
use std::io::{self, BufRead};
use std::ops::Range;
use std::path::{Path, PathBuf};

use qrnamap_core::de::{check_design, DeError, Shrink};
use qrnamap_core::{Contrast, CountMatrix, Design, Mode, RunConfig, SampleSheet, Strandedness};
use qseq_io::FastqReader;

use crate::de::{run_de, DeArgs, DeRunError};
use crate::io::{open_input, InputError};
use crate::pipeline::{
    align_sample, check_samples, prepare_align, prepare_quant, quant_sample, read_table, read_tx2gene,
    gene_list_from_gff, read_gene_list, remove_old_outputs, run_matrix_with_genes, run_transcripts, write_output, AlignContext, GenesFrom, PipelineError,
    QuantContext, SampleResult, Scratch,
};
use crate::slurm::{write_slurm_scripts, Plan, SlurmError};

/// Completion marker of the prepare step, in the run directory.
pub const PREPARED_MARKER: &str = ".prepared";
/// DESeq2 output prefix in the run directory.
pub const DE_PREFIX: &str = "deseq2";
/// Default design when `--contrast` is given without `--design`.
pub const DEFAULT_DESIGN: &str = "~ condition";

/// Which part of the run to do.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    All,
    Prepare,
    Final,
}

#[derive(Debug, Clone)]
pub struct RunArgs {
    pub samples: PathBuf,
    pub reference: PathBuf,
    pub gff: PathBuf,
    pub mode: Mode,
    pub out: PathBuf,
    /// Strandedness, fragment length, threads, max intron and reads per sample.
    pub config: RunConfig,
    pub tmpdir: Option<PathBuf>,
    pub design: Option<String>,
    pub contrast: Option<String>,
    pub shrink: Option<Shrink>,
    pub keep_sam: bool,
    pub step: Step,
}

/// SLURM options of `run --slurm`.
#[derive(Debug, Clone)]
pub struct SlurmOptions {
    pub per_task: usize,
    pub cpus: u32,
    pub mem_gb: u32,
    pub time: String,
    pub partition: Option<String>,
    pub account: Option<String>,
    pub qrnamap_bin: PathBuf,
    /// `--tool-path` directories, put first on the jobs' PATH (ruling R28).
    pub tool_paths: Vec<PathBuf>,
}

#[derive(Debug)]
pub enum RunError {
    /// The sample sheet is invalid (one line per problem; never FASTQ paths).
    Sheet(String),
    /// The design or contrast text is invalid, or `--design` came without `--contrast`.
    Design(String),
    /// The design does not fit the sheet (checked before any tool runs).
    DesignCheck(DeError),
    Input(InputError),
    Pipeline(PipelineError),
    De(DeRunError),
    Slurm(SlurmError),
    /// `--task-id` points past the last sample.
    TaskOutOfRange { task_id: usize, per_task: usize, n_samples: usize },
    /// `--task-id` and `--step final` need a run directory prepared by
    /// `qrnamap run --step prepare` from the same inputs (ruling R27). `reason` names what
    /// differs: a missing marker, the mode, the version, or an input's path, size or mtime.
    NotPrepared { dir: PathBuf, reason: String },
    /// A path with a newline cannot be recorded in the marker.
    UnsupportedPath(PathBuf),
    /// `--step final`: the sample's gene count file is older than the `.prepared` marker, so it
    /// was made before the current prepare step (final review I3).
    StaleCounts { sample: String, path: PathBuf },
    /// `--slurm` with a `--step` other than `all`.
    SlurmWithStep,
    /// A path is not valid UTF-8, so it cannot go into the scripts.
    NonUtf8Path(PathBuf),
    Io { path: PathBuf, source: io::Error },
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RunError::Sheet(m) | RunError::Design(m) => write!(f, "{m}"),
            RunError::DesignCheck(e) => write!(f, "design check: {e}"),
            RunError::Input(e) => write!(f, "{e}"),
            RunError::Pipeline(e) => write!(f, "{e}"),
            RunError::De(e) => write!(f, "{e}"),
            RunError::Slurm(e) => write!(f, "{e}"),
            RunError::TaskOutOfRange { task_id, per_task, n_samples } => write!(
                f,
                "task {task_id} with {per_task} sample(s) per task starts past the end of the sheet ({n_samples} samples)"
            ),
            RunError::NotPrepared { dir, reason } => write!(
                f,
                "run directory {} is not prepared for these inputs: {reason}. Run `qrnamap run --step prepare` with the same inputs first (it writes {PREPARED_MARKER})",
                dir.display()
            ),
            RunError::UnsupportedPath(p) => write!(f, "path {:?} contains a newline and cannot be recorded", p),
            RunError::StaleCounts { sample, path } => write!(
                f,
                "sample '{sample}': gene count file {} is older than {PREPARED_MARKER}; it was made before the current prepare step. Run quant or align for this sample again",
                path.display()
            ),
            RunError::SlurmWithStep => write!(f, "--slurm writes the scripts for all steps; do not combine it with --step"),
            RunError::NonUtf8Path(p) => write!(f, "path {} is not valid UTF-8", p.display()),
            RunError::Io { path, source } => write!(f, "{}: {source}", path.display()),
        }
    }
}

impl std::error::Error for RunError {}

impl From<PipelineError> for RunError {
    fn from(e: PipelineError) -> Self {
        RunError::Pipeline(e)
    }
}

/// What a local run produced.
#[derive(Debug, Clone, Default)]
pub struct RunOutput {
    pub counts: Option<PathBuf>,
    pub deseq2_version: Option<String>,
    pub de_files: Vec<PathBuf>,
}

/// Parse and validate a sample sheet file.
pub fn load_sheet(path: &Path) -> Result<SampleSheet, RunError> {
    let text = read_table(path)?;
    SampleSheet::parse(&text).map_err(|errs| {
        let lines: Vec<String> = errs.iter().map(|e| format!("  {e}")).collect();
        RunError::Sheet(format!("sample sheet {}:\n{}", path.display(), lines.join("\n")))
    })
}

/// Check the design and contrast against the sheet before any tool runs: the grammar, that the
/// contrast factor is a design term and a factor column, that level and reference differ, occur
/// in the sheet, and have at least 2 samples each. The all-zero check needs counts and runs in
/// `de` later.
pub fn check_de_plan(sheet: &SampleSheet, design: Option<&str>, contrast: Option<&str>) -> Result<(), RunError> {
    let Some(contrast) = contrast else {
        if design.is_some() {
            return Err(RunError::Design("--design needs --contrast".into()));
        }
        return Ok(());
    };
    let design = Design::parse(design.unwrap_or(DEFAULT_DESIGN), &sheet.columns)
        .map_err(|e| RunError::Design(format!("design: {e}")))?;
    let contrast = Contrast::parse(contrast, &sheet.columns).map_err(|e| RunError::Design(format!("contrast: {e}")))?;
    if !design.terms().contains(&contrast.factor) {
        return Err(RunError::DesignCheck(DeError::ContrastNotInDesign(contrast.factor)));
    }
    // A one-gene stand-in matrix with a count of 1 for every sample, so only the sheet-level
    // checks can fail here.
    let names: Vec<String> = sheet.samples.iter().map(|s| s.name.clone()).collect();
    let stand_in = CountMatrix { genes: vec!["stand-in".into()], counts: vec![vec![1; names.len()]], samples: names };
    check_design(sheet, &contrast, &stand_in).map_err(RunError::DesignCheck)
}

/// Samples `[task_id * per_task, min((task_id + 1) * per_task, n))`, in sheet order.
pub fn task_slice(n_samples: usize, task_id: usize, per_task: usize) -> Result<Range<usize>, RunError> {
    let out_of_range = || RunError::TaskOutOfRange { task_id, per_task, n_samples };
    if per_task == 0 {
        return Err(out_of_range());
    }
    let start = task_id.checked_mul(per_task).ok_or_else(out_of_range)?;
    if start >= n_samples {
        return Err(out_of_range());
    }
    Ok(start..(start + per_task).min(n_samples))
}

fn mode_word(m: Mode) -> &'static str {
    match m {
        Mode::Quant => "quant",
        Mode::Align => "align",
    }
}

fn io_err(path: &Path) -> impl FnOnce(io::Error) -> RunError + '_ {
    move |source| RunError::Io { path: path.to_path_buf(), source }
}

fn remove_if_present(p: &Path) -> Result<(), RunError> {
    match fs::remove_file(p) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(RunError::Io { path: p.to_path_buf(), source: e }),
    }
}

/// First line of the marker file.
const MARKER_HEADER: &str = "qrnamap-prepared\t1";

/// One recorded input: role, canonical path, size in bytes, mtime (seconds, nanoseconds).
#[derive(Debug, Clone, PartialEq, Eq)]
struct InputStamp {
    role: String,
    path: PathBuf,
    size: u64,
    secs: u64,
    nanos: u32,
}

fn stamp(role: &str, path: &Path) -> Result<InputStamp, String> {
    let canon = fs::canonicalize(path).map_err(|e| format!("{role}: cannot read {}: {e}", path.display()))?;
    let meta = fs::metadata(&canon).map_err(|e| format!("{role}: cannot read {}: {e}", canon.display()))?;
    let t = meta
        .modified()
        .map_err(|e| format!("{role}: no mtime for {}: {e}", canon.display()))?
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|_| format!("{role}: mtime of {} is before 1970", canon.display()))?;
    Ok(InputStamp { role: role.to_string(), path: canon, size: meta.len(), secs: t.as_secs(), nanos: t.subsec_nanos() })
}

/// Write `<out>/.prepared` (ruling R27). Only the prepare step calls this, last, after the
/// index is built. The file is written through `.partial` and a rename.
///
/// Format (tab-separated lines):
/// ```text
/// qrnamap-prepared<TAB>1
/// version<TAB><qrnamap version>
/// mode<TAB>quant|align
/// input<TAB><role><TAB><size bytes><TAB><mtime seconds><TAB><mtime nanoseconds><TAB><canonical path>
/// ```
/// One `input` line per entry of `inputs`, in order. Roles: `reference`, `gff`, and in quant
/// mode `transcripts` and `tx2gene`. The path is the last field, so it may hold a tab; a path
/// with a newline is refused (`UnsupportedPath`).
pub fn write_prepared_marker(out: &Path, mode: Mode, inputs: &[(&str, &Path)]) -> Result<(), RunError> {
    let mut text = format!("{MARKER_HEADER}\nversion\t{}\nmode\t{}\n", env!("CARGO_PKG_VERSION"), mode_word(mode));
    for (role, path) in inputs {
        let st = stamp(role, path).map_err(|reason| RunError::NotPrepared { dir: out.to_path_buf(), reason })?;
        let ps = st.path.to_str().ok_or_else(|| RunError::NonUtf8Path(st.path.clone()))?;
        if ps.contains(['\n', '\r']) {
            return Err(RunError::UnsupportedPath(st.path.clone()));
        }
        text.push_str(&format!("input\t{}\t{}\t{}\t{}\t{ps}\n", st.role, st.size, st.secs, st.nanos));
    }
    write_output(&out.join(PREPARED_MARKER), text.as_bytes())?;
    Ok(())
}

/// Remove `<out>/.prepared`. Plain `quant` and `align` call this before they build an index in
/// `<out>`, so a marker never describes an index it did not see built.
pub fn remove_prepared_marker(out: &Path) -> Result<(), RunError> {
    remove_if_present(&out.join(PREPARED_MARKER))
}

fn parse_marker(text: &str) -> Result<(String, String, Vec<InputStamp>), String> {
    let mut lines = text.lines();
    if lines.next() != Some(MARKER_HEADER) {
        return Err(format!("{PREPARED_MARKER} has an unknown format"));
    }
    let (mut version, mut mode, mut inputs) = (None, None, Vec::new());
    for l in lines {
        let f: Vec<&str> = l.splitn(6, '\t').collect();
        match f.as_slice() {
            ["version", v] => version = Some(v.to_string()),
            ["mode", m] => mode = Some(m.to_string()),
            ["input", role, size, secs, nanos, path] => {
                let bad = || format!("{PREPARED_MARKER}: bad input line for {role}");
                inputs.push(InputStamp {
                    role: role.to_string(),
                    path: PathBuf::from(path),
                    size: size.parse().map_err(|_| bad())?,
                    secs: secs.parse().map_err(|_| bad())?,
                    nanos: nanos.parse().map_err(|_| bad())?,
                });
            }
            _ => return Err(format!("{PREPARED_MARKER} has an unknown line")),
        }
    }
    match (version, mode) {
        (Some(v), Some(m)) => Ok((v, m, inputs)),
        _ => Err(format!("{PREPARED_MARKER} lacks the version or mode line")),
    }
}

/// Compare `<out>/.prepared` with this run (ruling R27): the mode, the qrnamap version, and for
/// each input the role set, the canonical path, the current size and the current mtime. Any
/// difference is `NotPrepared` with a reason that names the input and what changed.
pub fn check_prepared(out: &Path, mode: Mode, inputs: &[(&str, &Path)]) -> Result<(), RunError> {
    let not = |reason: String| RunError::NotPrepared { dir: out.to_path_buf(), reason };
    let marker = out.join(PREPARED_MARKER);
    let text = match fs::read_to_string(&marker) {
        Ok(t) => t,
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Err(not(format!("marker {} is missing", marker.display()))),
        Err(e) => return Err(not(format!("marker {}: {e}", marker.display()))),
    };
    let (version, prep_mode, recorded) = parse_marker(&text).map_err(not)?;
    if prep_mode != mode_word(mode) {
        return Err(not(format!("mode differs: prepared for {prep_mode}, this run is {}", mode_word(mode))));
    }
    if version != env!("CARGO_PKG_VERSION") {
        return Err(not(format!("version differs: prepared by qrnamap {version}, this is qrnamap {}", env!("CARGO_PKG_VERSION"))));
    }
    let roles = |v: Vec<&str>| v.join(", ");
    let rec_roles: Vec<&str> = recorded.iter().map(|r| r.role.as_str()).collect();
    let given_roles: Vec<&str> = inputs.iter().map(|(r, _)| *r).collect();
    if rec_roles != given_roles {
        return Err(not(format!("inputs differ: prepared with [{}], given [{}]", roles(rec_roles), roles(given_roles))));
    }
    for (rec, (role, path)) in recorded.iter().zip(inputs) {
        let now = stamp(role, path).map_err(not)?;
        if now.path != rec.path {
            return Err(not(format!("{role}: path differs: prepared from {}, given {}", rec.path.display(), now.path.display())));
        }
        if now.size != rec.size {
            return Err(not(format!("{role}: size changed: {} bytes when prepared, now {} bytes", rec.size, now.size)));
        }
        if (now.secs, now.nanos) != (rec.secs, rec.nanos) {
            return Err(not(format!(
                "{role}: mtime changed: {}.{:09} when prepared, now {}.{:09}",
                rec.secs, rec.nanos, now.secs, now.nanos
            )));
        }
    }
    Ok(())
}

/// Final review I3: refuse a sample count file whose mtime is older than the marker. "Older
/// than" (not "not newer than") so a coarse mtime resolution cannot refuse a file made right
/// after the marker. A missing count file is left to the matrix step (`MissingCountFile`).
pub fn check_counts_not_stale(out: &Path, sheet: &SampleSheet) -> Result<(), RunError> {
    let marker = out.join(PREPARED_MARKER);
    let mtime = |p: &Path| fs::metadata(p).and_then(|m| m.modified());
    let prepared = mtime(&marker).map_err(io_err(&marker))?;
    for s in &sheet.samples {
        let path = out.join(format!("{}.genecounts.tsv.zst", s.name));
        match mtime(&path) {
            Ok(t) if t < prepared => return Err(RunError::StaleCounts { sample: s.name.clone(), path }),
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(e) => return Err(RunError::Io { path, source: e }),
        }
    }
    Ok(())
}

fn require_file(out: &Path, p: &Path) -> Result<(), RunError> {
    if p.is_file() {
        Ok(())
    } else {
        Err(RunError::NotPrepared { dir: out.to_path_buf(), reason: format!("{} is missing", p.display()) })
    }
}

/// The marker inputs of a run directory, in marker order.
fn marker_inputs<'a>(mode: Mode, reference: &'a Path, gff: &'a Path, tx: Option<(&'a Path, &'a Path)>) -> Vec<(&'static str, &'a Path)> {
    let mut v = vec![("reference", reference), ("gff", gff)];
    if let (Mode::Quant, Some((fa, t2g))) = (mode, tx) {
        v.extend([("transcripts", fa), ("tx2gene", t2g)]);
    }
    v
}

/// The kallisto context of a prepared run directory, without building anything. The marker
/// must match `reference`, `gff`, `transcripts` and `tx2gene` (ruling R27).
pub fn prepared_quant(
    out: &Path,
    reference: &Path,
    gff: &Path,
    transcripts: &Path,
    tx2gene: &Path,
    scratch: Scratch,
    config: RunConfig,
) -> Result<QuantContext, RunError> {
    let out = std::path::absolute(out).map_err(io_err(out))?;
    check_prepared(&out, Mode::Quant, &marker_inputs(Mode::Quant, reference, gff, Some((transcripts, tx2gene))))?;
    let index = out.join("transcripts.kidx");
    require_file(&out, &index)?;
    Ok(QuantContext { index, tx2gene: read_tx2gene(tx2gene)?, out_dir: out, scratch, config })
}

/// The HISAT2 + featureCounts context of a prepared run directory, without building anything.
/// The marker must match `reference` and `gff` (ruling R27).
pub fn prepared_align(
    out: &Path,
    reference: &Path,
    gff: &Path,
    scratch: Scratch,
    config: RunConfig,
    keep_sam: bool,
) -> Result<AlignContext, RunError> {
    let out = std::path::absolute(out).map_err(io_err(out))?;
    check_prepared(&out, Mode::Align, &marker_inputs(Mode::Align, reference, gff, None))?;
    let index_prefix = out.join("index").join("genome");
    let splice_sites = out.join("splice_sites.txt");
    let saf = out.join("genes.saf");
    for p in [&out.join("index").join("genome.1.ht2"), &splice_sites, &saf] {
        require_file(&out, p)?;
    }
    Ok(AlignContext { index_prefix, splice_sites, saf, out_dir: out, scratch, config, keep_sam })
}

/// Sheet, builder and design checks shared by every step. No file is written.
fn check_inputs(args: &RunArgs) -> Result<SampleSheet, RunError> {
    let sheet = load_sheet(&args.samples)?;
    let mut cfg = args.config.clone();
    cfg.mode = args.mode;
    check_samples(&sheet.samples, &cfg, args.mode)?;
    check_de_plan(&sheet, args.design.as_deref(), args.contrast.as_deref())?;
    Ok(sheet)
}

/// Remove what an earlier run left: the marker, the outputs of every sample, the matrix and the
/// DESeq2 files. Called before the first tool runs.
fn remove_old_run_outputs(out: &Path, sheet: &SampleSheet) -> Result<(), RunError> {
    remove_if_present(&out.join(PREPARED_MARKER))?;
    for s in &sheet.samples {
        remove_old_outputs(out, &s.name)?;
    }
    remove_if_present(&out.join("counts.tsv.zst"))?;
    remove_if_present(&out.join("counts.tsv.zst.partial"))?;
    crate::de::remove_old_outputs(&out.join(DE_PREFIX)).map_err(RunError::De)?;
    Ok(())
}

enum Prepared {
    Quant(QuantContext),
    Align(AlignContext),
}

fn prepare(args: &RunArgs, out: &Path, scratch: &Scratch, warn: &mut dyn FnMut(&str)) -> Result<(SampleSheet, Prepared), RunError> {
    let sheet = check_inputs(args)?;
    for w in check_reads_per_sample(&sheet, args.config.reads_per_sample).map_err(RunError::Input)? {
        warn(&w.to_string());
    }
    fs::create_dir_all(out).map_err(io_err(out))?;
    remove_old_run_outputs(out, &sheet)?;
    let mut cfg = args.config.clone();
    cfg.mode = args.mode;
    let ctx = match args.mode {
        Mode::Quant => {
            let (fa, t2g) = (out.join("transcripts.fa.zst"), out.join("tx2gene.tsv.zst"));
            run_transcripts(&args.reference, &args.gff, &fa, &t2g, None, None)?;
            Prepared::Quant(prepare_quant(&fa, &t2g, out, scratch.clone(), cfg)?)
        }
        Mode::Align => Prepared::Align(prepare_align(&args.reference, &args.gff, out, scratch.clone(), cfg, args.keep_sam)?),
    };
    let (fa, t2g) = (out.join("transcripts.fa.zst"), out.join("tx2gene.tsv.zst"));
    write_prepared_marker(out, args.mode, &marker_inputs(args.mode, &args.reference, &args.gff, Some((&fa, &t2g))))?;
    Ok((sheet, ctx))
}

fn finish(args: &RunArgs, out: &Path, sheet: &SampleSheet, warn: &mut dyn FnMut(&str)) -> Result<RunOutput, RunError> {
    // Rows: quant from the tx2gene table the marker checked; align from the --gff given.
    let genes = match args.mode {
        Mode::Quant => read_gene_list(GenesFrom::Tx2Gene, &out.join("tx2gene.tsv.zst"))?,
        Mode::Align => gene_list_from_gff(&args.gff)?,
    };
    let (counts, warnings) = run_matrix_with_genes(out, &sheet.samples, Some(&genes))?;
    for w in &warnings {
        warn(w);
    }
    let mut res = RunOutput { counts: Some(counts.clone()), ..RunOutput::default() };
    if let Some(contrast) = &args.contrast {
        let de = DeArgs {
            counts,
            samples: args.samples.clone(),
            design: args.design.clone().unwrap_or_else(|| DEFAULT_DESIGN.to_string()),
            contrast: contrast.clone(),
            shrink: args.shrink,
            out: out.join(DE_PREFIX),
            tmpdir: args.tmpdir.clone(),
        };
        let r = run_de(&de).map_err(RunError::De)?;
        res.deseq2_version = Some(r.deseq2_version);
        res.de_files = r.files;
    }
    Ok(res)
}

/// Run the steps of `args.step` on this machine. `on_sample` is called after each sample and
/// `warn` for each warning.
pub fn run_local(
    args: &RunArgs,
    on_sample: &mut dyn FnMut(&SampleResult),
    warn: &mut dyn FnMut(&str),
) -> Result<RunOutput, RunError> {
    let out = std::path::absolute(&args.out).map_err(io_err(&args.out))?;
    let scratch = Scratch::from_flag_or_env(args.tmpdir.clone());
    let result = match args.step {
        Step::Final => {
            let sheet = check_inputs(args)?;
            let (fa, t2g) = (out.join("transcripts.fa.zst"), out.join("tx2gene.tsv.zst"));
            check_prepared(&out, args.mode, &marker_inputs(args.mode, &args.reference, &args.gff, Some((&fa, &t2g))))?;
            check_counts_not_stale(&out, &sheet)?;
            finish(args, &out, &sheet, warn)
        }
        Step::Prepare => prepare(args, &out, &scratch, warn).map(|_| RunOutput::default()),
        Step::All => {
            let (sheet, ctx) = prepare(args, &out, &scratch, warn)?;
            for s in &sheet.samples {
                let r = match &ctx {
                    Prepared::Quant(c) => quant_sample(s, c)?,
                    Prepared::Align(c) => align_sample(s, c)?,
                };
                on_sample(&r);
            }
            finish(args, &out, &sheet, warn)
        }
    };
    scratch.cleanup();
    result
}

fn utf8(p: &Path) -> Result<String, RunError> {
    p.to_str().map(str::to_string).ok_or_else(|| RunError::NonUtf8Path(p.to_path_buf()))
}

fn strand_word(s: Strandedness) -> &'static str {
    match s {
        Strandedness::Unstranded => "unstranded",
        Strandedness::Forward => "forward",
        Strandedness::Reverse => "reverse",
    }
}

/// The programs the jobs start by name: kallisto (quant); hisat2-build (a Python script, so also
/// `python`, and the `hisat2-build-s` it runs), hisat2-align-s and featureCounts (align);
/// Rscript with a contrast. The qrnamap binary is called by its absolute path.
pub fn required_tools(mode: Mode, with_de: bool) -> Vec<&'static str> {
    let mut v = match mode {
        Mode::Quant => vec!["kallisto"],
        Mode::Align => vec!["hisat2-build", "hisat2-build-s", "hisat2-align-s", "python", "featureCounts"],
    };
    if with_de {
        v.push("Rscript");
    }
    v
}

/// The SLURM plan for `args`. Paths become absolute. Sample names and FASTQ paths are not in it.
pub fn slurm_plan(args: &RunArgs, opts: &SlurmOptions, n_samples: usize) -> Result<Plan, RunError> {
    let abs = |p: &Path| std::path::absolute(p).map_err(io_err(p));
    let (samples, reference, gff, out) = (abs(&args.samples)?, abs(&args.reference)?, abs(&args.gff)?, abs(&args.out)?);
    let work_dir = std::env::current_dir().map_err(io_err(Path::new(".")))?;
    let c = &args.config;
    let mut common = vec!["--strandedness".to_string(), strand_word(c.strandedness).into()];
    if let Some((m, sd)) = c.fragment_len {
        common.extend(["--fragment-length".into(), format!("{m},{sd}")]);
    }
    common.extend(["--threads".into(), c.threads.to_string()]);

    let mut run_args = vec![
        "--samples".to_string(),
        utf8(&samples)?,
        "--reference".into(),
        utf8(&reference)?,
        "--gff".into(),
        utf8(&gff)?,
        "--mode".into(),
        mode_word(args.mode).into(),
        "--out".into(),
        utf8(&out)?,
    ];
    run_args.extend(common.iter().cloned());
    run_args.extend([
        "--max-intron".into(),
        c.max_intron.to_string(),
        "--reads-per-sample".into(),
        c.reads_per_sample.to_string(),
    ]);
    if let Some(ct) = &args.contrast {
        run_args.extend(["--design".into(), args.design.clone().unwrap_or_else(|| DEFAULT_DESIGN.into()), "--contrast".into(), ct.clone()]);
        if let Some(s) = args.shrink {
            run_args.extend(["--shrink".into(), s.name().into()]);
        }
    }
    if args.keep_sam {
        run_args.push("--keep-sam".into());
    }

    let mut task_args = vec!["--samples".to_string(), utf8(&samples)?];
    match args.mode {
        Mode::Quant => task_args.extend([
            "--reference".into(),
            utf8(&reference)?,
            "--gff".into(),
            utf8(&gff)?,
            "--transcripts".into(),
            utf8(&out.join("transcripts.fa.zst"))?,
            "--tx2gene".into(),
            utf8(&out.join("tx2gene.tsv.zst"))?,
        ]),
        Mode::Align => task_args.extend(["--reference".into(), utf8(&reference)?, "--gff".into(), utf8(&gff)?]),
    }
    task_args.extend(["--out".into(), utf8(&out)?]);
    task_args.extend(common);
    if args.mode == Mode::Align {
        task_args.extend(["--max-intron".into(), c.max_intron.to_string()]);
        if args.keep_sam {
            task_args.push("--keep-sam".into());
        }
    }
    let tools = required_tools(args.mode, args.contrast.is_some());
    let tool_dirs = crate::slurm::resolve_tool_dirs(&tools, &opts.tool_paths, std::env::var_os("PATH").as_deref())
        .map_err(RunError::Slurm)?;
    Ok(Plan {
        n_samples,
        per_task: opts.per_task,
        cpus: opts.cpus,
        mem_gb: opts.mem_gb,
        time: opts.time.clone(),
        partition: opts.partition.clone(),
        account: opts.account.clone(),
        run_dir: out,
        work_dir,
        qrnamap_bin: abs(&opts.qrnamap_bin)?,
        mode: args.mode,
        run_args,
        task_args,
        tool_dirs,
        tools: tools.iter().map(|t| t.to_string()).collect(),
    })
}

/// `run --slurm`: check the inputs (sheet, builders, design), then write the scripts. No tool
/// runs and no read is counted (the prepare job does the reads check). Returns `submit.sh`.
pub fn run_slurm(args: &RunArgs, opts: &SlurmOptions) -> Result<PathBuf, RunError> {
    if args.step != Step::All {
        return Err(RunError::SlurmWithStep);
    }
    let sheet = check_inputs(args)?;
    let plan = slurm_plan(args, opts, sheet.samples.len())?;
    write_slurm_scripts(&plan).map_err(RunError::Slurm)
}

/// Turn a read error into an `InputError`, keeping a wrapped `Truncated` or `Corrupt` error.
pub(crate) fn input_error(path: &Path, e: io::Error) -> InputError {
    if e.get_ref().is_some_and(|inner| inner.is::<InputError>()) {
        return *e.into_inner().and_then(|b| b.downcast::<InputError>().ok()).expect("checked with is::<InputError>()");
    }
    InputError::Io { path: path.to_path_buf(), source: e }
}

/// Count the FASTQ records in `reader`. With `limit`, stop as soon as the count is above the
/// limit and return `limit + 1`, so a large file is not read to its end. A truncated or bad
/// record is an error. `path` is used in error messages.
pub fn count_records<R: BufRead>(reader: R, path: &Path, limit: Option<u64>) -> Result<u64, InputError> {
    let mut n = 0u64;
    for rec in FastqReader::new(reader) {
        rec.map_err(|e| input_error(path, e))?;
        n += 1;
        if limit.is_some_and(|l| n > l) {
            break;
        }
    }
    Ok(n)
}

/// Count the records (read pairs, when `path` is a mate 1 file) of a plain, gzip or zstd FASTQ.
/// See `count_records` for `limit`.
pub fn count_pairs(path: &Path, limit: Option<u64>) -> Result<u64, InputError> {
    count_records(open_input(path)?, path, limit)
}

/// A sample with more read pairs than `reads_per_sample`. In the CLI this is a warning, never
/// an error (spec 11a).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Warning {
    pub sample: String,
    pub limit: u64,
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "sample '{}' has more than {} read pairs (reads_per_sample = {}); the run continues. Raise --reads-per-sample to silence this, or check the SLURM resources",
            self.sample, self.limit, self.limit
        )
    }
}

/// A warning for each sample whose `fastq1` files (all lanes, summed) hold more than `limit`
/// records. Reading stops as soon as a sample is over the limit.
pub fn check_reads_per_sample(sheet: &SampleSheet, limit: u64) -> Result<Vec<Warning>, InputError> {
    let mut out = Vec::new();
    for s in &sheet.samples {
        let mut total = 0u64;
        for lane in &s.fastq1 {
            total += count_pairs(Path::new(lane), Some(limit - total))?;
            if total > limit {
                out.push(Warning { sample: s.name.clone(), limit });
                break;
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn task_slice_covers_the_sheet_in_order() {
        assert_eq!(task_slice(7, 0, 3).unwrap(), 0..3);
        assert_eq!(task_slice(7, 1, 3).unwrap(), 3..6);
        assert_eq!(task_slice(7, 2, 3).unwrap(), 6..7);
        assert!(matches!(task_slice(7, 3, 3), Err(RunError::TaskOutOfRange { .. })));
        assert!(matches!(task_slice(6, 2, 3), Err(RunError::TaskOutOfRange { .. })));
        assert!(matches!(task_slice(6, 0, 0), Err(RunError::TaskOutOfRange { .. })));
        assert!(matches!(task_slice(6, usize::MAX, 2), Err(RunError::TaskOutOfRange { .. })));
        assert_eq!(task_slice(1, 0, 1).unwrap(), 0..1);
    }

    #[test]
    fn design_without_contrast_is_an_error() {
        let sheet = SampleSheet::parse("sample\tfastq1\tcondition\na\t/x\tc\n").unwrap();
        assert!(check_de_plan(&sheet, None, None).is_ok());
        assert!(matches!(check_de_plan(&sheet, Some("~ condition"), None), Err(RunError::Design(_))));
    }
}
