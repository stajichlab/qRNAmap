//! `qrnamap` command line. Argument parsing and printing only; the work is in `pipeline`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand, ValueEnum};
use qrnamap_cli::pipeline::{
    self, align_sample, check_samples, prepare_align, prepare_quant, quant_sample, run_matrix,
    run_transcripts, strand_check_sample, GenesFrom, SampleResult, Scratch,
    LOW_ASSIGNED_FRACTION,
};
use qrnamap_core::{Mode, RunConfig, Sample, SampleSheet, Strandedness};

#[derive(Parser)]
#[command(name = "qrnamap", version, about = "Bulk RNA-seq quantification (kallisto or HISAT2 + featureCounts)")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
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
    },
    /// kallisto index and quant per sample; writes DIR/<sample>.genecounts.tsv.zst.
    Quant {
        #[arg(long)]
        samples: PathBuf,
        #[arg(long)]
        transcripts: PathBuf,
        #[arg(long)]
        tx2gene: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[command(flatten)]
        run: RunFlags,
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
        /// HISAT2 --max-intronlen. The default 5000 is a chosen value, UNVERIFIED for fungi.
        #[arg(long, default_value_t = 5000)]
        max_intron: u32,
        /// Keep a zstd copy of each SAM as DIR/<sample>/<sample>.sam.zst.
        #[arg(long)]
        keep_sam: bool,
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
        #[arg(long)]
        tmpdir: Option<PathBuf>,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum StrandArg {
    Unstranded,
    Forward,
    Reverse,
}

#[derive(Args)]
struct RunFlags {
    #[arg(long, value_enum, default_value = "unstranded")]
    strandedness: StrandArg,
    /// Fragment length mean and SD for single-end kallisto, as MEAN,SD.
    #[arg(long, value_parser = parse_fragment_length)]
    fragment_length: Option<(f64, f64)>,
    #[arg(long, default_value_t = 1)]
    threads: u32,
    /// Temporary directory. Default: $SCRATCH.
    #[arg(long)]
    tmpdir: Option<PathBuf>,
}

fn parse_fragment_length(s: &str) -> Result<(f64, f64), String> {
    let (a, b) = s.split_once(',').ok_or("expected MEAN,SD")?;
    let p = |x: &str| x.trim().parse::<f64>().map_err(|_| format!("{x:?} is not a number"));
    Ok((p(a)?, p(b)?))
}

impl RunFlags {
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
    SampleSheet::parse(&text).map(|s| s.samples).map_err(|errs| {
        let lines: Vec<String> = errs.iter().map(|e| format!("  {e}")).collect();
        format!("sample sheet {}:\n{}", path.display(), lines.join("\n"))
    })
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
        Cmd::Transcripts { reference, gff, out, tx2gene, saf, splice_sites } => {
            let n = run_transcripts(&reference, &gff, &out, &tx2gene, saf.as_deref(), splice_sites.as_deref()).map_err(e)?;
            eprintln!("wrote {n} transcripts to {}", out.display());
        }
        Cmd::Quant { samples, transcripts, tx2gene, out, run } => {
            let samples = load_sheet(&samples)?;
            let cfg = run.config(Mode::Quant);
            check_samples(&samples, &cfg, Mode::Quant).map_err(e)?;
            let scratch = Scratch::from_flag_or_env(run.tmpdir.clone());
            let ctx = prepare_quant(&transcripts, &tx2gene, &out, scratch.clone(), cfg).map_err(e)?;
            for s in &samples {
                report(&quant_sample(s, &ctx).map_err(e)?, "pseudoaligned");
            }
            scratch.cleanup();
        }
        Cmd::Align { samples, reference, gff, out, run, max_intron, keep_sam } => {
            let samples = load_sheet(&samples)?;
            let mut cfg = run.config(Mode::Align);
            cfg.max_intron = max_intron;
            check_samples(&samples, &cfg, Mode::Align).map_err(e)?;
            let scratch = Scratch::from_flag_or_env(run.tmpdir.clone());
            let ctx = prepare_align(&reference, &gff, &out, scratch.clone(), cfg, keep_sam).map_err(e)?;
            for s in &samples {
                report(&align_sample(s, &ctx).map_err(e)?, "mapped");
            }
            scratch.cleanup();
        }
        Cmd::Matrix { out, samples, genes_from } => {
            let samples = load_sheet(&samples)?;
            let gf = match genes_from.as_deref() {
                None => None,
                Some([kind, file]) => {
                    let k = match kind.as_str() {
                        "tx2gene" => GenesFrom::Tx2Gene,
                        "saf" => GenesFrom::Saf,
                        other => return Err(format!("--genes-from: kind must be tx2gene or saf, not {other:?}")),
                    };
                    Some((k, PathBuf::from(file)))
                }
                Some(_) => return Err("--genes-from takes KIND FILE".into()),
            };
            let (path, warnings) =
                run_matrix(&out, &samples, gf.as_ref().map(|(k, p)| (*k, p.as_path()))).map_err(e)?;
            for w in &warnings {
                eprintln!("warning: {w}");
            }
            eprintln!("wrote {}", path.display());
        }
        Cmd::StrandCheck { samples, reference, gff, out, pairs, max_intron, threads, tmpdir } => {
            let samples = load_sheet(&samples)?;
            let cfg = RunConfig { mode: Mode::Align, max_intron, threads, ..RunConfig::default() };
            check_samples(&samples, &cfg, Mode::Align).map_err(e)?;
            let scratch = Scratch::from_flag_or_env(tmpdir);
            let ctx = prepare_align(&reference, &gff, &out, scratch.clone(), cfg, false).map_err(e)?;
            let mut table = String::from(
                "sample\tsense\tantisense\tambiguous\tsense_fraction\tantisense_fraction\tambiguous_fraction\tsuggested\tlow_fraction_warning\n",
            );
            for s in &samples {
                let r = strand_check_sample(s, &ctx, pairs).map_err(e)?;
                let suggested = format!("{:?}", r.suggested).to_lowercase();
                let row = format!(
                    "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{suggested}\t{}\n",
                    s.name, r.sense, r.antisense, r.ambiguous, r.sense_fraction, r.antisense_fraction,
                    r.ambiguous_fraction, r.low_fraction_warning
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
