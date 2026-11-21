//! `qrnamap` command line. Argument parsing and printing only; the work is in `pipeline`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand, ValueEnum};
use qrnamap_cli::pipeline::{
    self, align_sample, check_fastq_files, check_samples, prepare_align, prepare_quant,
    quant_sample, run_matrix, run_transcripts_with_genes, strand_check_sample, GenesFrom,
    GenesSource, SampleResult, Scratch, LOW_ASSIGNED_FRACTION,
};
use qrnamap_core::{Mode, Quantifier, RunConfig, Sample, SampleSheet, Strandedness};

#[derive(Parser)]
#[command(
    name = "qrnamap",
    version,
    about = "Bulk RNA-seq quantification (kallisto, salmon, or HISAT2 + featureCounts)"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Download a genome FASTA and its GFF3 (NCBI from the genome catalog, or FungiDB).
    Fetch {
        #[arg(long, value_enum)]
        source: SourceArg,
        /// NCBI: a catalog accession (for example GCF_000146045.2). FungiDB: an organism name.
        #[arg(long)]
        id: String,
        /// Output directory. NCBI files keep the catalog names; FungiDB files get a .gz name.
        #[arg(long)]
        out: PathBuf,
        /// FungiDB API key file. It wins over $EUPATHDB_KEY when both are given.
        #[arg(long)]
        eupathdb_key_file: Option<PathBuf>,
        /// FungiDB site URL (FungiDB only). Default: https://fungidb.org. Tests point it at a
        /// local server.
        #[arg(long)]
        base_url: Option<String>,
        /// Test only: read the genome catalog from this JSON file instead of the built-in one.
        #[arg(long, hide = true)]
        catalog: Option<PathBuf>,
    },
    /// Transcript FASTA and tx2gene (optionally SAF and HISAT2 splice sites) from a genome and GFF3.
    Transcripts {
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        gff: PathBuf,
        /// Transcript FASTA. Compressed by extension: .zst, .gz, else plain.
        #[arg(long, default_value = "transcripts.fa.zst")]
        out: PathBuf,
        #[arg(long, default_value = "tx2gene.tsv.zst")]
        tx2gene: PathBuf,
        #[arg(long)]
        saf: Option<PathBuf>,
        #[arg(long)]
        splice_sites: Option<PathBuf>,
        /// Gene table (class, biotype, mitochondrial flag per gene). Compressed by extension.
        #[arg(long)]
        genes_table: Option<PathBuf>,
        /// A contig to flag mitochondrial in genes.tsv.zst (repeatable). Added to the contigs
        /// the GFF3 marks with a `region` feature `genome=mitochondrion`.
        #[arg(long = "mito-contig", value_name = "NAME", allow_hyphen_values = true, value_parser = parse_contig_name)]
        mito_contig: Vec<String>,
    },
    /// kallisto or salmon index and quant per sample; writes DIR/<sample>.genecounts.tsv.zst.
    Quant {
        #[arg(long)]
        samples: PathBuf,
        /// Pseudo-alignment tool. Default: kallisto. With --task-id it must match the
        /// quantifier the run directory was prepared with.
        #[arg(long, value_enum)]
        quantifier: Option<QuantArg>,
        /// Salmon only: pass --softclipOverhangs to salmon quant. Needs --quantifier salmon.
        #[arg(long)]
        salmon_softclip_overhangs: bool,
        /// Salmon only: pass --minScoreFraction F to salmon quant. F is greater than 0 and at
        /// most 1. Needs --quantifier salmon.
        #[arg(long, value_name = "F", value_parser = parse_min_score_fraction, allow_hyphen_values = true)]
        salmon_min_score_fraction: Option<f64>,
        #[arg(long)]
        transcripts: PathBuf,
        #[arg(long)]
        tx2gene: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[command(flatten)]
        run: RunFlags,
        #[command(flatten)]
        task: TaskFlags,
        /// With --task-id (required then): the genome FASTA the run directory was prepared
        /// from. Checked against OUT/.prepared.
        #[arg(long, requires = "task_id")]
        reference: Option<PathBuf>,
        /// With --task-id (required then): the GFF3 the run directory was prepared from.
        /// Without --task-id: optional; when given, OUT/genes.tsv.zst is written from it.
        #[arg(long)]
        gff: Option<PathBuf>,
        /// A contig to flag mitochondrial in genes.tsv.zst (repeatable). Added to the contigs
        /// the GFF3 marks with a `region` feature `genome=mitochondrion`.
        #[arg(long = "mito-contig", value_name = "NAME", allow_hyphen_values = true, value_parser = parse_contig_name)]
        mito_contig: Vec<String>,
    },
    /// HISAT2 index, then HISAT2 and featureCounts per sample; writes DIR/<sample>.genecounts.tsv.zst.
    Align {
        #[arg(long)]
        samples: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        gff: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[command(flatten)]
        run: RunFlags,
        #[command(flatten)]
        task: TaskFlags,
        /// HISAT2 --max-intronlen. The default 5000 is a chosen value, UNVERIFIED for fungi.
        #[arg(long, default_value_t = 5000)]
        max_intron: u32,
        /// Keep a zstd copy of each SAM as DIR/<sample>/<sample>.sam.zst.
        #[arg(long)]
        keep_sam: bool,
        /// A contig to flag mitochondrial in genes.tsv.zst (repeatable). Added to the contigs
        /// the GFF3 marks with a `region` feature `genome=mitochondrion`.
        #[arg(long = "mito-contig", value_name = "NAME", allow_hyphen_values = true, value_parser = parse_contig_name)]
        mito_contig: Vec<String>,
    },
    /// All steps: reference files and index, per-sample counts, matrix, and DESeq2 when
    /// --contrast is given. With --slurm, writes SLURM scripts to OUT/slurm/ and runs no tool.
    Run {
        #[arg(long)]
        samples: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        gff: PathBuf,
        #[arg(long, value_enum)]
        mode: ModeArg,
        /// Pseudo-alignment tool (quant mode only; an error with --mode align). Default:
        /// kallisto.
        #[arg(long, value_enum)]
        quantifier: Option<QuantArg>,
        /// Salmon only: pass --softclipOverhangs to salmon quant. Needs --quantifier salmon.
        #[arg(long)]
        salmon_softclip_overhangs: bool,
        /// Salmon only: pass --minScoreFraction F to salmon quant. F is greater than 0 and at
        /// most 1. Needs --quantifier salmon.
        #[arg(long, value_name = "F", value_parser = parse_min_score_fraction, allow_hyphen_values = true)]
        salmon_min_score_fraction: Option<f64>,
        /// Run directory. On a cluster it must be on shared storage (for example /bigdata).
        #[arg(long)]
        out: PathBuf,
        #[command(flatten)]
        run: RunFlags,
        /// HISAT2 --max-intronlen (align mode). The default 5000 is a chosen value, UNVERIFIED for fungi.
        #[arg(long, default_value_t = 5000)]
        max_intron: u32,
        /// Design formula over sample sheet columns. Default with --contrast: '~ condition'.
        #[arg(long, requires = "contrast")]
        design: Option<String>,
        /// FACTOR,LEVEL,REFERENCE. Runs DESeq2 after the matrix; writes OUT/deseq2.*.
        #[arg(long)]
        contrast: Option<String>,
        /// Also write OUT/deseq2.results.shrunk.csv with lfcShrink. normal needs no extra
        /// package. apeglm and ashr need the R packages apeglm and ashr, which the pixi
        /// environment installs.
        #[arg(long, value_enum, requires = "contrast")]
        shrink: Option<ShrinkArg>,
        /// Read pairs per sample above which a warning is printed (spec 11a: a warning in the CLI).
        #[arg(long, default_value_t = 20_000_000)]
        reads_per_sample: u64,
        /// Keep a zstd copy of each SAM (align mode).
        #[arg(long)]
        keep_sam: bool,
        /// Run only a part: prepare (reference files, index, removal of old outputs) or final
        /// (matrix and DESeq2). The SLURM scripts use these.
        #[arg(long, value_enum, default_value = "all")]
        step: StepArg,
        /// Write SLURM scripts (prepare, sample array, final, submit.sh) instead of running.
        #[arg(long)]
        slurm: bool,
        /// Samples per array task. Default 1 until a measured per-sample time exists (aim for
        /// about 1 to 1.5 hours of work per task).
        #[arg(long, default_value_t = 1, requires = "slurm")]
        samples_per_task: usize,
        /// CPUs per SLURM task. Default: the --threads value. UNMEASURED placeholder.
        #[arg(long, requires = "slurm")]
        cpus: Option<u32>,
        /// Memory per SLURM task in GB. The default is an UNMEASURED placeholder and is small on
        /// purpose: on the AMB data HISAT2 used about 40 MiB (41,208 KB) and featureCounts about 35 MiB (36,144 KB) resident, and two
        /// 4 GB jobs gave identical counts. No job ran at 8 GB (tested: 32 GB and 4 GB); the 4 GB runs
        /// took 10:32 and 19:26 against 6:28 at 32 GB on different nodes, so the effect on run time is
        /// open (validation/amb/README.md). Yeast align tasks peaked at 6.7 to 8.0 GB MaxRSS. Not
        /// measured for genomes near 50 Mb or for 20 million read pairs per sample: raise it for those.
        #[arg(long, default_value_t = qrnamap_cli::slurm::DEFAULT_MEM_GB, requires = "slurm")]
        mem_gb: u32,
        /// SLURM time limit per task. UNMEASURED placeholder.
        #[arg(long, default_value = qrnamap_cli::slurm::DEFAULT_TIME, requires = "slurm")]
        time: String,
        #[arg(long, requires = "slurm")]
        partition: Option<String>,
        #[arg(long, requires = "slurm")]
        account: Option<String>,
        /// The qrnamap binary the jobs call. Default: this program's absolute path.
        #[arg(long, requires = "slurm")]
        qrnamap_bin: Option<PathBuf>,
        /// A directory with tools for the jobs (repeatable). These go first on the jobs' PATH,
        /// then the directories where the needed tools are found on the current PATH. The
        /// scripts do not depend on the PATH of the shell that runs submit.sh.
        #[arg(long, requires = "slurm")]
        tool_path: Vec<PathBuf>,
        /// A contig to flag mitochondrial in genes.tsv.zst (repeatable). Added to the contigs
        /// the GFF3 marks with a `region` feature `genome=mitochondrion`.
        #[arg(long = "mito-contig", value_name = "NAME", allow_hyphen_values = true, value_parser = parse_contig_name)]
        mito_contig: Vec<String>,
    },
    /// Gene-by-sample matrix DIR/counts.tsv.zst from DIR/<sample>.genecounts.tsv.zst.
    Matrix {
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        samples: PathBuf,
        /// Row set and order from a tx2gene table or a SAF file: `--genes-from tx2gene FILE`.
        #[arg(long, num_args = 2, value_names = ["KIND", "FILE"])]
        genes_from: Option<Vec<String>>,
    },
    /// DESeq2 differential expression (Wald test) with Rscript.
    De {
        /// counts.tsv from `qrnamap matrix` (plain, .gz or .zst).
        #[arg(long)]
        counts: PathBuf,
        /// Sample sheet (TSV with a header row). Its covariate columns are the variables for
        /// --design and --contrast. Factor values use only [A-Za-z0-9_.-]. A header `name:numeric`
        /// makes a numeric column.
        #[arg(long)]
        samples: PathBuf,
        /// Design formula over sample sheet columns.
        #[arg(long, default_value = "~ condition")]
        design: String,
        /// FACTOR,LEVEL,REFERENCE.
        #[arg(long)]
        contrast: String,
        /// Also write PREFIX.results.shrunk.csv with lfcShrink. Default: none. normal needs no
        /// extra package. apeglm and ashr need the R packages apeglm and ashr, which the pixi
        /// environment installs.
        #[arg(long, value_enum)]
        shrink: Option<ShrinkArg>,
        /// Output prefix: PREFIX.results.csv, .normalized.csv, .vst.csv, .versions.txt.
        #[arg(long)]
        out: PathBuf,
        /// Temporary directory. Default: $SCRATCH, then $TMPDIR (ruling R23).
        #[arg(long)]
        tmpdir: Option<PathBuf>,
    },
    /// Estimate library strandedness from the first N pairs of each sample.
    StrandCheck {
        #[arg(long)]
        samples: PathBuf,
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        gff: PathBuf,
        #[arg(long)]
        out: PathBuf,
        /// Read pairs (reads, for single-end) taken from the start of each sample.
        #[arg(long, default_value_t = 100_000)]
        pairs: u64,
        #[arg(long, default_value_t = 5000)]
        max_intron: u32,
        #[arg(long, default_value_t = 1)]
        threads: u32,
        /// Temporary directory. Default: $SCRATCH, then $TMPDIR (ruling R23).
        #[arg(long)]
        tmpdir: Option<PathBuf>,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum ModeArg {
    Quant,
    Align,
}

#[derive(Clone, Copy, ValueEnum)]
enum QuantArg {
    Kallisto,
    Salmon,
}

impl QuantArg {
    fn core(self) -> Quantifier {
        match self {
            QuantArg::Kallisto => Quantifier::Kallisto,
            QuantArg::Salmon => Quantifier::Salmon,
        }
    }
}

/// The quantifier of `--quantifier`, or the default (kallisto) when the flag is not given.
fn quantifier_or_default(q: Option<QuantArg>) -> Quantifier {
    q.map_or(Quantifier::Kallisto, QuantArg::core)
}

/// Parse `--salmon-min-score-fraction`: finite and in (0, 1]. This runs before any command
/// work, so an invalid value is refused before any file is written (R20).
fn parse_min_score_fraction(x: &str) -> Result<f64, String> {
    let v: f64 = x
        .trim()
        .parse()
        .map_err(|_| format!("--salmon-min-score-fraction: {x:?} is not a number"))?;
    if v.is_finite() && v > 0.0 && v <= 1.0 {
        Ok(v)
    } else {
        Err(format!(
            "--salmon-min-score-fraction must be finite and in (0, 1], got {x}"
        ))
    }
}

/// Set the Salmon options on `cfg`. A Salmon option needs `--quantifier salmon` (quant mode).
/// `align_mode` is true for `run --mode align`.
fn apply_salmon_opts(
    cfg: &mut RunConfig,
    quantifier: Option<QuantArg>,
    softclip: bool,
    min_score_fraction: Option<f64>,
    align_mode: bool,
) -> Result<(), String> {
    let mut given = Vec::new();
    if softclip {
        given.push("--salmon-softclip-overhangs");
    }
    if min_score_fraction.is_some() {
        given.push("--salmon-min-score-fraction");
    }
    if given.is_empty() {
        return Ok(());
    }
    let names = given.join(" and ");
    if align_mode {
        return Err(format!(
            "{names} cannot be used with --mode align; it applies to --mode quant \
             with --quantifier salmon"
        ));
    }
    if quantifier_or_default(quantifier) != Quantifier::Salmon {
        return Err(format!("{names} needs --quantifier salmon"));
    }
    cfg.salmon_softclip_overhangs = softclip;
    cfg.salmon_min_score_fraction = min_score_fraction;
    Ok(())
}

#[derive(Clone, Copy, ValueEnum)]
enum StepArg {
    All,
    Prepare,
    Final,
}

/// SLURM array task flags of `quant` and `align`.
#[derive(Args)]
struct TaskFlags {
    /// Process only the samples of this array task: samples [ID*K, (ID+1)*K) in sheet order.
    /// The index is not built: it is taken from a run directory prepared by
    /// `qrnamap run --step prepare` (it needs OUT/.prepared).
    #[arg(long, value_name = "ID")]
    task_id: Option<usize>,
    /// K, the samples per array task (with --task-id).
    #[arg(long, default_value_t = 1, requires = "task_id")]
    samples_per_task: usize,
}

#[derive(Clone, Copy, ValueEnum)]
enum SourceArg {
    Ncbi,
    Fungidb,
}

#[derive(Clone, Copy, ValueEnum)]
enum StrandArg {
    Unstranded,
    Forward,
    Reverse,
}

#[derive(Clone, Copy, ValueEnum)]
enum ShrinkArg {
    Normal,
    Apeglm,
    Ashr,
}

#[derive(Args)]
struct RunFlags {
    #[arg(long, value_enum, default_value = "unstranded")]
    strandedness: StrandArg,
    /// Fragment length mean and SD for single-end kallisto or salmon, as MEAN,SD.
    #[arg(long, value_parser = parse_fragment_length)]
    fragment_length: Option<(f64, f64)>,
    #[arg(long, default_value_t = 1)]
    threads: u32,
    /// Temporary directory. Default: $SCRATCH, then $TMPDIR (ruling R23). With `run --slurm`
    /// the directory is written into the scripts instead of the $SCRATCH check. The jobs
    /// create it (qrnamap runs create_dir_all, with a `qrnamap-<pid>` subdirectory) when
    /// the parent is writable. Use storage that every compute node can reach (for example
    /// /bigdata). A /scratch or $SCRATCH path gives a warning: it is job-specific.
    #[arg(long)]
    tmpdir: Option<PathBuf>,
}

impl ShrinkArg {
    fn core(self) -> qrnamap_core::de::Shrink {
        match self {
            ShrinkArg::Normal => qrnamap_core::de::Shrink::Normal,
            ShrinkArg::Apeglm => qrnamap_core::de::Shrink::Apeglm,
            ShrinkArg::Ashr => qrnamap_core::de::Shrink::Ashr,
        }
    }
}

impl TaskFlags {
    /// The slice of the sheet this invocation processes: all samples without --task-id.
    fn range(&self, n: usize) -> Result<std::ops::Range<usize>, String> {
        match self.task_id {
            None => Ok(0..n),
            Some(id) => qrnamap_cli::run::task_slice(n, id, self.samples_per_task)
                .map_err(|x| x.to_string()),
        }
    }
}

/// A `--mito-contig` value: any non-empty text. It is only compared with GFF seqids.
fn parse_contig_name(s: &str) -> Result<String, String> {
    if s.is_empty() {
        return Err("a contig name must not be empty".into());
    }
    Ok(s.to_string())
}

fn parse_fragment_length(s: &str) -> Result<(f64, f64), String> {
    let (a, b) = s.split_once(',').ok_or("expected MEAN,SD")?;
    let p = |x: &str| {
        x.trim()
            .parse::<f64>()
            .map_err(|_| format!("{x:?} is not a number"))
    };
    Ok((p(a)?, p(b)?))
}

impl RunFlags {
    /// The quantifier is not set here: `RunFlags` is also flattened into `align`, which has no
    /// `--quantifier`. The `quant` and `run` call sites set `quantifier` after this call.
    fn config(&self, mode: Mode) -> RunConfig {
        RunConfig {
            mode,
            strandedness: match self.strandedness {
                StrandArg::Unstranded => Strandedness::Unstranded,
                StrandArg::Forward => Strandedness::Forward,
                StrandArg::Reverse => Strandedness::Reverse,
            },
            fragment_len: self.fragment_length,
            threads: self.threads,
            ..RunConfig::default()
        }
    }
}

fn load_sheet(path: &Path) -> Result<Vec<Sample>, String> {
    let text = pipeline::read_table(path).map_err(|e| e.to_string())?;
    SampleSheet::parse(&text)
        .map(|s| s.samples)
        .map_err(|errs| {
            let lines: Vec<String> = errs.iter().map(|e| format!("  {e}")).collect();
            format!("sample sheet {}:\n{}", path.display(), lines.join("\n"))
        })
}

/// Remove the old outputs of every sheet sample (ruling R24, final review I2). Plain `quant`
/// and `align` call this before they build an index, as `run` does.
fn remove_all_old_outputs(out: &Path, samples: &[Sample]) -> Result<(), pipeline::PipelineError> {
    if !out.is_dir() {
        return Ok(());
    }
    pipeline::remove_old_genes_table(out)?;
    pipeline::remove_old_mapping_summary(out)?;
    samples
        .iter()
        .try_for_each(|s| pipeline::remove_old_outputs(out, &s.name))
}

/// Label of the mapped fraction in the per-sample stdout line: kallisto pseudoaligns, salmon maps.
fn quant_label(q: Quantifier) -> &'static str {
    match q {
        Quantifier::Kallisto => "pseudoaligned",
        Quantifier::Salmon => "mapped",
    }
}

fn report(r: &SampleResult, what: &str) {
    let s = &r.summary;
    println!(
        "{}\treads_processed={}\t{what}={:.2}%\tassigned_fraction={:.4}",
        r.sample,
        s.reads_processed,
        s.fraction_mapped * 100.0,
        s.fraction_assigned
    );
    for w in &r.warnings {
        eprintln!("warning: {w}");
    }
    if s.fraction_assigned < LOW_ASSIGNED_FRACTION {
        eprintln!(
            "warning: sample '{}': assigned fraction {:.4} is below {LOW_ASSIGNED_FRACTION}; check the strandedness setting (qrnamap strand-check)",
            r.sample, s.fraction_assigned
        );
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let e = |x: pipeline::PipelineError| x.to_string();
    match cli.cmd {
        Cmd::Fetch {
            source,
            id,
            out,
            eupathdb_key_file,
            base_url,
            catalog,
        } => {
            use qrnamap_cli::fetch_cmd::{run_fetch, FetchArgs, Source};
            let source = match source {
                SourceArg::Ncbi => Source::Ncbi,
                SourceArg::Fungidb => Source::Fungidb,
            };
            let args = FetchArgs {
                source,
                id,
                out,
                key_file: eupathdb_key_file,
                base_url,
                catalog,
            };
            let r = run_fetch(&args).map_err(|x| x.to_string())?;
            for w in &r.warnings {
                eprintln!("warning: {w}");
            }
            println!("fasta\t{}", r.fetched.fasta.display());
            if let Some(g) = &r.fetched.gff {
                println!("gff\t{}", g.display());
            }
        }
        Cmd::Transcripts {
            reference,
            gff,
            out,
            tx2gene,
            saf,
            splice_sites,
            genes_table,
            mito_contig,
        } => {
            if genes_table.is_none() && !mito_contig.is_empty() {
                eprintln!("note: --mito-contig has no effect without --genes-table (no gene table is written)");
            }
            let genes = genes_table.as_deref().map(|p| (p, mito_contig.as_slice()));
            let n = run_transcripts_with_genes(
                &reference,
                &gff,
                &out,
                &tx2gene,
                saf.as_deref(),
                splice_sites.as_deref(),
                genes,
            )
            .map_err(e)?;
            eprintln!("wrote {n} transcripts to {}", out.display());
        }
        Cmd::Quant {
            samples,
            quantifier,
            salmon_softclip_overhangs,
            salmon_min_score_fraction,
            transcripts,
            tx2gene,
            out,
            run,
            task,
            reference,
            gff,
            mito_contig,
        } => {
            let samples = load_sheet(&samples)?;
            let mut cfg = run.config(Mode::Quant);
            cfg.quantifier = quantifier_or_default(quantifier);
            apply_salmon_opts(
                &mut cfg,
                quantifier,
                salmon_softclip_overhangs,
                salmon_min_score_fraction,
                false,
            )?;
            check_samples(&samples, &cfg, Mode::Quant).map_err(e)?;
            let range = task.range(samples.len())?;
            let scratch = Scratch::from_flag_or_env(run.tmpdir.clone());
            if task.task_id.is_none() && gff.is_none() {
                if !mito_contig.is_empty() {
                    return Err(
                        "quant: --mito-contig needs --gff (the gene table comes from the GFF3)"
                            .into(),
                    );
                }
                eprintln!(
                    "note: no --gff given, so {}/genes.tsv.zst is not written",
                    out.display()
                );
            }
            if task.task_id.is_some() && !mito_contig.is_empty() {
                eprintln!("note: --mito-contig has no effect with --task-id (the prepare step writes genes.tsv.zst)");
            }
            let ctx = if task.task_id.is_some() {
                let (Some(r), Some(g)) = (&reference, &gff) else {
                    return Err("quant --task-id needs --reference and --gff (the inputs the run directory was prepared from)".into());
                };
                let ctx = qrnamap_cli::run::prepared_quant(
                    &out,
                    r,
                    g,
                    &transcripts,
                    &tx2gene,
                    scratch.clone(),
                    cfg,
                )
                .map_err(|x| x.to_string())?;
                check_fastq_files(&samples).map_err(e)?;
                ctx
            } else {
                check_fastq_files(&samples).map_err(e)?;
                // Ruling R27: an index built here was not made by `run --step prepare`.
                qrnamap_cli::run::remove_prepared_marker(&out).map_err(|x| x.to_string())?;
                // Final review I2: a failing sample must not leave later samples' old counts.
                remove_all_old_outputs(&out, &samples).map_err(e)?;
                let genes = gff.as_deref().map(|g| GenesSource {
                    gff: g,
                    mito_contigs: &mito_contig,
                });
                prepare_quant(&transcripts, &tx2gene, &out, scratch.clone(), cfg, genes)
                    .map_err(e)?
            };
            for s in &samples[range] {
                report(
                    &quant_sample(s, &ctx).map_err(e)?,
                    quant_label(ctx.config.quantifier),
                );
            }
            scratch.cleanup();
        }
        Cmd::Align {
            samples,
            reference,
            gff,
            out,
            run,
            task,
            max_intron,
            keep_sam,
            mito_contig,
        } => {
            let samples = load_sheet(&samples)?;
            let mut cfg = run.config(Mode::Align);
            cfg.max_intron = max_intron;
            check_samples(&samples, &cfg, Mode::Align).map_err(e)?;
            let range = task.range(samples.len())?;
            let scratch = Scratch::from_flag_or_env(run.tmpdir.clone());
            if task.task_id.is_some() && !mito_contig.is_empty() {
                eprintln!("note: --mito-contig has no effect with --task-id (the prepare step writes genes.tsv.zst)");
            }
            let ctx = if task.task_id.is_some() {
                let ctx = qrnamap_cli::run::prepared_align(
                    &out,
                    &reference,
                    &gff,
                    scratch.clone(),
                    cfg,
                    keep_sam,
                )
                .map_err(|x| x.to_string())?;
                check_fastq_files(&samples).map_err(e)?;
                ctx
            } else {
                check_fastq_files(&samples).map_err(e)?;
                // Ruling R27: an index built here was not made by `run --step prepare`.
                qrnamap_cli::run::remove_prepared_marker(&out).map_err(|x| x.to_string())?;
                // Final review I2: a failing sample must not leave later samples' old counts.
                remove_all_old_outputs(&out, &samples).map_err(e)?;
                prepare_align(
                    &reference,
                    &gff,
                    &out,
                    scratch.clone(),
                    cfg,
                    keep_sam,
                    &mito_contig,
                )
                .map_err(e)?
            };
            for s in &samples[range] {
                report(&align_sample(s, &ctx).map_err(e)?, "mapped");
            }
            scratch.cleanup();
        }
        Cmd::Run {
            samples,
            reference,
            gff,
            mode,
            quantifier,
            salmon_softclip_overhangs,
            salmon_min_score_fraction,
            out,
            run,
            max_intron,
            design,
            contrast,
            shrink,
            reads_per_sample,
            keep_sam,
            step,
            slurm,
            samples_per_task,
            cpus,
            mem_gb,
            time,
            partition,
            account,
            qrnamap_bin,
            tool_path,
            mito_contig,
        } => {
            use qrnamap_cli::run::{
                run_local, run_slurm, scratch_out_warning, scratch_tmpdir_warning, RunArgs,
                SlurmOptions, Step,
            };
            // Spec decision 7: --quantifier with align mode is an error, not ignored.
            if matches!(mode, ModeArg::Align) && quantifier.is_some() {
                let m = "--quantifier cannot be used with --mode align; \
                         it applies to --mode quant only";
                return Err(m.into());
            }
            let align_mode = matches!(mode, ModeArg::Align);
            let abs_out =
                std::path::absolute(&out).map_err(|x| format!("{}: {x}", out.display()))?;
            if let Some(w) = scratch_out_warning(&abs_out, std::env::var_os("SCRATCH").as_deref()) {
                eprintln!("warning: {w}");
            }
            if slurm {
                if let Some(t) = &run.tmpdir {
                    let abs_t =
                        std::path::absolute(t).map_err(|x| format!("{}: {x}", t.display()))?;
                    if let Some(w) =
                        scratch_tmpdir_warning(&abs_t, std::env::var_os("SCRATCH").as_deref())
                    {
                        eprintln!("warning: {w}");
                    }
                }
            }
            let mode = match mode {
                ModeArg::Quant => Mode::Quant,
                ModeArg::Align => Mode::Align,
            };
            let mut config = run.config(mode);
            config.quantifier = quantifier_or_default(quantifier);
            apply_salmon_opts(
                &mut config,
                quantifier,
                salmon_softclip_overhangs,
                salmon_min_score_fraction,
                align_mode,
            )?;
            config.max_intron = max_intron;
            config.reads_per_sample = reads_per_sample;
            let step = match step {
                StepArg::All => Step::All,
                StepArg::Prepare => Step::Prepare,
                StepArg::Final => Step::Final,
            };
            let args = RunArgs {
                samples,
                reference,
                gff,
                mode,
                out,
                config,
                tmpdir: run.tmpdir.clone(),
                design,
                contrast,
                shrink: shrink.map(ShrinkArg::core),
                keep_sam,
                step,
                mito_contigs: mito_contig,
            };
            if slurm {
                let qrnamap_bin = match qrnamap_bin {
                    Some(p) => p,
                    None => std::env::current_exe()
                        .map_err(|x| format!("cannot find the path of this program: {x}"))?,
                };
                let opts = SlurmOptions {
                    per_task: samples_per_task,
                    cpus: cpus.unwrap_or(run.threads),
                    mem_gb,
                    time,
                    partition,
                    account,
                    qrnamap_bin,
                    tool_paths: tool_path,
                };
                let submit = run_slurm(&args, &opts).map_err(|x| x.to_string())?;
                let q = qrnamap_cli::slurm::sh_quote(&submit.to_string_lossy())
                    .map_err(|x| x.to_string())?;
                println!(
                    "wrote SLURM scripts to {}",
                    submit.parent().unwrap_or(&submit).display()
                );
                println!("submit with:\n  bash {q}");
                println!("The --cpus, --mem-gb and --time values are UNMEASURED placeholders; check them before submitting.");
            } else {
                let what = if mode == Mode::Quant {
                    quant_label(args.config.quantifier)
                } else {
                    "mapped"
                };
                let r = run_local(&args, &mut |r| report(r, what), &mut |w| {
                    eprintln!("warning: {w}")
                })
                .map_err(|x| x.to_string())?;
                if let Some(p) = &r.counts {
                    eprintln!("wrote {}", p.display());
                }
                if let Some(v) = &r.deseq2_version {
                    println!("DESeq2 {v}");
                }
                for f in &r.de_files {
                    eprintln!("wrote {}", f.display());
                }
            }
        }
        Cmd::Matrix {
            out,
            samples,
            genes_from,
        } => {
            let samples = load_sheet(&samples)?;
            let gf = match genes_from.as_deref() {
                None => None,
                Some([kind, file]) => {
                    let k = match kind.as_str() {
                        "tx2gene" => GenesFrom::Tx2Gene,
                        "saf" => GenesFrom::Saf,
                        other => {
                            return Err(format!(
                                "--genes-from: kind must be tx2gene or saf, not {other:?}"
                            ))
                        }
                    };
                    Some((k, PathBuf::from(file)))
                }
                Some(_) => return Err("--genes-from takes KIND FILE".into()),
            };
            let (path, warnings) =
                run_matrix(&out, &samples, gf.as_ref().map(|(k, p)| (*k, p.as_path())))
                    .map_err(e)?;
            for w in &warnings {
                eprintln!("warning: {w}");
            }
            eprintln!("wrote {}", path.display());
        }
        Cmd::De {
            counts,
            samples,
            design,
            contrast,
            shrink,
            out,
            tmpdir,
        } => {
            let shrink = shrink.map(ShrinkArg::core);
            let args = qrnamap_cli::de::DeArgs {
                counts,
                samples,
                design,
                contrast,
                shrink,
                out,
                tmpdir,
            };
            let r = qrnamap_cli::de::run_de(&args).map_err(|x| x.to_string())?;
            println!("DESeq2 {}", r.deseq2_version);
            for f in &r.files {
                eprintln!("wrote {}", f.display());
            }
        }
        Cmd::StrandCheck {
            samples,
            reference,
            gff,
            out,
            pairs,
            max_intron,
            threads,
            tmpdir,
        } => {
            let samples = load_sheet(&samples)?;
            let cfg = RunConfig {
                mode: Mode::Align,
                max_intron,
                threads,
                ..RunConfig::default()
            };
            check_samples(&samples, &cfg, Mode::Align).map_err(e)?;
            check_fastq_files(&samples).map_err(e)?;
            // Ruling R33 (R24 order: check inputs, remove stale outputs, then build): an old
            // table must not survive a failed rerun and be joined into the mapping summary.
            for name in ["strand_check.tsv", "strand_check.tsv.partial"] {
                match std::fs::remove_file(out.join(name)) {
                    Ok(()) => {}
                    Err(x) if x.kind() == std::io::ErrorKind::NotFound => {}
                    Err(x) => {
                        return Err(e(pipeline::PipelineError::Io {
                            path: out.join(name),
                            source: x,
                        }))
                    }
                }
            }
            let scratch = Scratch::from_flag_or_env(tmpdir);
            // Final review I1: the index, genes.saf and splice_sites.txt built below replace the
            // prepared ones in OUT, so the marker no longer describes them.
            qrnamap_cli::run::remove_prepared_marker(&out).map_err(|x| x.to_string())?;
            let ctx = pipeline::prepare_strand_check(&reference, &gff, &out, scratch.clone(), cfg)
                .map_err(e)?;
            let mut table = String::from(
                "sample\tsense\tantisense\tambiguous\tsense_fraction\tantisense_fraction\tambiguous_fraction\tsuggested\tlow_fraction_warning\n",
            );
            for s in &samples {
                let r = strand_check_sample(s, &ctx, pairs).map_err(e)?;
                let suggested = format!("{:?}", r.suggested).to_lowercase();
                let row = format!(
                    "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{suggested}\t{}\n",
                    s.name,
                    r.sense,
                    r.antisense,
                    r.ambiguous,
                    r.sense_fraction,
                    r.antisense_fraction,
                    r.ambiguous_fraction,
                    r.low_fraction_warning
                );
                print!("{row}");
                table.push_str(&row);
                if r.low_fraction_warning {
                    eprintln!(
                        "warning: sample '{}': neither strand holds more than {} of the stranded-assigned fragments (sense {:.4}, antisense {:.4}); the library looks unstranded, or too few reads were assigned",
                        s.name, qrnamap_core::STRAND_THRESHOLD, r.sense_fraction, r.antisense_fraction
                    );
                }
            }
            pipeline::write_output(&out.join("strand_check.tsv"), table.as_bytes()).map_err(e)?;
            scratch.cleanup();
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::FAILURE
        }
    }
}
