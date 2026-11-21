//! Tool calls as data. Builders make argv lists. Nothing here runs a program or does I/O.
//!
//! Safety contract (ruling R18): paths and sample names become single argv
//! elements only. No builder or lowering function makes a shell string, so a
//! space, quote or `$` in a name is never interpreted. Path parameters are
//! `&str` so that every argv element is valid UTF-8 by construction.

use crate::config::{RunConfig, Strandedness};
use crate::sheet::Sample;
use std::fmt;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    KallistoIndex,
    KallistoQuant,
    SalmonIndex,
    SalmonQuant,
    Hisat2Build,
    Hisat2Align,
    FeatureCounts,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ToolCall {
    pub tool: Tool,
    /// Arguments after the program name.
    pub args: Vec<String>,
    pub inputs: Vec<PathBuf>,
    pub outputs: Vec<PathBuf>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ConfigError {
    /// Single-end kallisto and salmon need a fragment length mean and SD.
    FragmentLengthRequired { sample: String },
    /// Fragment length mean or SD is not finite or not > 0 (single-end kallisto and salmon).
    InvalidFragmentLength { sample: String },
    /// `threads` is 0.
    InvalidThreads,
    /// Paired sample whose R1 and R2 lists differ in length.
    MateCountMismatch { sample: String },
    /// Sample has no R1 files.
    NoReads { sample: String },
    /// A FASTQ path contains a comma. HISAT2 lanes are comma-joined, so it would split.
    CommaInPath { sample: String, path: String },
    /// Salmon `--minScoreFraction` is not finite or not in (0, 1].
    InvalidMinScoreFraction { value: f64 },
}

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConfigError::FragmentLengthRequired { sample } => write!(
                f,
                "sample '{sample}' is single-end: the quantifier needs a fragment length mean and SD"
            ),
            ConfigError::InvalidFragmentLength { sample } => write!(
                f,
                "sample '{sample}': fragment length mean and SD must be finite and greater than 0"
            ),
            ConfigError::InvalidThreads => write!(f, "threads must be at least 1"),
            ConfigError::NoReads { sample } => write!(f, "sample '{sample}' has no R1 files"),
            ConfigError::CommaInPath { sample, path } => write!(
                f,
                "sample '{sample}': FASTQ path {path:?} contains a comma, which HISAT2 cannot take"
            ),
            ConfigError::InvalidMinScoreFraction { value } => write!(
                f,
                "--salmon-min-score-fraction must be finite and in (0, 1], got {value}"
            ),
            ConfigError::MateCountMismatch { sample } => {
                write!(f, "sample '{sample}' has different numbers of R1 and R2 files")
            }
        }
    }
}

impl std::error::Error for ConfigError {}

fn program(tool: Tool) -> &'static str {
    match tool {
        Tool::KallistoIndex | Tool::KallistoQuant => "kallisto",
        Tool::SalmonIndex | Tool::SalmonQuant => "salmon",
        Tool::Hisat2Build => "hisat2-build",
        // Not the `hisat2` wrapper: that Perl script joins its arguments into one string and
        // runs it with `system`, so a shell reads paths and names (a space splits a path, `$HOME`
        // is expanded; seen in the Task 11 tests). `hisat2-align-s` is the aligner binary the
        // wrapper calls. The `-s` (small index) binary matches the `.ht2` files that
        // `hisat2_build` lists; see the genome size note there.
        Tool::Hisat2Align => "hisat2-align-s",
        Tool::FeatureCounts => "featureCounts",
    }
}

fn check_threads(cfg: &RunConfig) -> Result<(), ConfigError> {
    if cfg.threads == 0 {
        return Err(ConfigError::InvalidThreads);
    }
    Ok(())
}

/// R1 must be non-empty. When R2 is non-empty, lane counts must match.
fn check_reads(sample: &Sample) -> Result<(), ConfigError> {
    if sample.fastq1.is_empty() {
        return Err(ConfigError::NoReads {
            sample: sample.name.clone(),
        });
    }
    if !sample.fastq2.is_empty() && sample.fastq1.len() != sample.fastq2.len() {
        return Err(ConfigError::MateCountMismatch {
            sample: sample.name.clone(),
        });
    }
    Ok(())
}

/// Fragment length mean and SD for a single-end sample. Both must be finite and > 0.
fn fragment_len(sample: &Sample, cfg: &RunConfig) -> Result<(f64, f64), ConfigError> {
    let (mean, sd) = cfg
        .fragment_len
        .ok_or_else(|| ConfigError::FragmentLengthRequired {
            sample: sample.name.clone(),
        })?;
    if !(mean.is_finite() && sd.is_finite() && mean > 0.0 && sd > 0.0) {
        return Err(ConfigError::InvalidFragmentLength {
            sample: sample.name.clone(),
        });
    }
    Ok((mean, sd))
}

fn st(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| a.to_string()).collect()
}

/// `kallisto index -i <index> <transcripts.fa>`.
pub fn kallisto_index(transcripts_fasta: &str, index: &str) -> ToolCall {
    ToolCall {
        tool: Tool::KallistoIndex,
        args: st(&["index", "-i", index, transcripts_fasta]),
        inputs: vec![transcripts_fasta.into()],
        outputs: vec![index.into()],
    }
}

/// `kallisto quant`. Paired lanes are passed as adjacent pairs in lane order (R19).
pub fn kallisto_quant(
    sample: &Sample,
    index: &str,
    out: &str,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    check_threads(cfg)?;
    check_reads(sample)?;
    let paired = !sample.fastq2.is_empty();
    let mut args = st(&["quant", "-i", index, "-o", out, "--plaintext", "-t"]);
    args.push(cfg.threads.to_string());
    match cfg.strandedness {
        Strandedness::Forward => args.push("--fr-stranded".into()),
        Strandedness::Reverse => args.push("--rf-stranded".into()),
        Strandedness::Unstranded => {}
    }
    let mut files: Vec<&String> = Vec::new();
    if paired {
        for (a, b) in sample.fastq1.iter().zip(&sample.fastq2) {
            files.push(a);
            files.push(b);
        }
    } else {
        let (mean, sd) = fragment_len(sample, cfg)?;
        // f64 Display gives the shortest round-trip form: 200.0 -> "200".
        args.extend([
            "--single".into(),
            "-l".into(),
            mean.to_string(),
            "-s".into(),
            sd.to_string(),
        ]);
        files.extend(sample.fastq1.iter());
    }
    args.extend(files.iter().map(|f| f.to_string()));
    let mut inputs: Vec<PathBuf> = files.iter().map(|f| PathBuf::from(f.as_str())).collect();
    inputs.push(index.into());
    Ok(ToolCall {
        tool: Tool::KallistoQuant,
        args,
        inputs,
        outputs: vec![
            PathBuf::from(out).join("abundance.tsv"),
            PathBuf::from(out).join("run_info.json"),
        ],
    })
}

/// `salmon index -t <transcripts.fa> -i <index_dir> -p <threads>
/// --keepDuplicates`. The index is a directory.
/// `--keepDuplicates` keeps transcripts with identical sequences, so each keeps its own row.
pub fn salmon_index(
    transcripts_fasta: &str,
    index_dir: &str,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    check_threads(cfg)?;
    Ok(ToolCall {
        tool: Tool::SalmonIndex,
        args: vec![
            "index".into(),
            "-t".into(),
            transcripts_fasta.into(),
            "-i".into(),
            index_dir.into(),
            "-p".into(),
            cfg.threads.to_string(),
            "--keepDuplicates".into(),
        ],
        inputs: vec![transcripts_fasta.into()],
        outputs: vec![index_dir.into()],
    })
}

/// `salmon quant`. Library type comes from strandedness and layout. Lanes are passed as an
/// R1 list after `-1` and an R2 list after `-2`, in lane order.
pub fn salmon_quant(
    sample: &Sample,
    index_dir: &str,
    out_dir: &str,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    check_threads(cfg)?;
    check_reads(sample)?;
    if let Some(v) = cfg.salmon_min_score_fraction {
        if !(v.is_finite() && v > 0.0 && v <= 1.0) {
            return Err(ConfigError::InvalidMinScoreFraction { value: v });
        }
    }
    let paired = !sample.fastq2.is_empty();
    let lib = match (paired, cfg.strandedness) {
        (true, Strandedness::Unstranded) => "IU",
        (true, Strandedness::Forward) => "ISF",
        (true, Strandedness::Reverse) => "ISR",
        (false, Strandedness::Unstranded) => "U",
        (false, Strandedness::Forward) => "SF",
        (false, Strandedness::Reverse) => "SR",
    };
    let mut args = st(&["quant", "-i", index_dir, "-l", lib]);
    let files: Vec<&String> = if paired {
        args.push("-1".into());
        args.extend(sample.fastq1.iter().cloned());
        args.push("-2".into());
        args.extend(sample.fastq2.iter().cloned());
        sample.fastq1.iter().chain(&sample.fastq2).collect()
    } else {
        let (mean, sd) = fragment_len(sample, cfg)?;
        args.push("-r".into());
        args.extend(sample.fastq1.iter().cloned());
        args.extend([
            "--fldMean".into(),
            mean.to_string(),
            "--fldSD".into(),
            sd.to_string(),
        ]);
        sample.fastq1.iter().collect()
    };
    if cfg.salmon_softclip_overhangs {
        args.push("--softclipOverhangs".into());
    }
    if let Some(v) = cfg.salmon_min_score_fraction {
        args.extend(["--minScoreFraction".into(), v.to_string()]);
    }
    args.extend([
        "-p".into(),
        cfg.threads.to_string(),
        "-o".into(),
        out_dir.into(),
    ]);
    let mut inputs: Vec<PathBuf> = files.iter().map(|f| PathBuf::from(f.as_str())).collect();
    inputs.push(index_dir.into());
    Ok(ToolCall {
        tool: Tool::SalmonQuant,
        args,
        inputs,
        outputs: vec![
            PathBuf::from(out_dir).join("quant.sf"),
            PathBuf::from(out_dir)
                .join("aux_info")
                .join("meta_info.json"),
        ],
    })
}

/// `hisat2-build -p <threads> <reference.fa> <index_prefix>`.
///
/// The eight-file output list `<prefix>.1.ht2` to `.8.ht2` is correct only for genomes
/// below about 4 Gb. HISAT2 writes `.ht2l` files above that. The 50 Mb cap of the
/// genome catalog keeps the list valid.
pub fn hisat2_build(
    reference_fasta: &str,
    index_prefix: &str,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    check_threads(cfg)?;
    Ok(ToolCall {
        tool: Tool::Hisat2Build,
        args: vec![
            "-p".into(),
            cfg.threads.to_string(),
            reference_fasta.into(),
            index_prefix.into(),
        ],
        inputs: vec![reference_fasta.into()],
        outputs: (1..=8)
            .map(|i| PathBuf::from(format!("{index_prefix}.{i}.ht2")))
            .collect(),
    })
}

/// `hisat2` alignment with known splice sites only. Lanes are comma-joined lists.
/// `index_prefix` is listed as an input as given. The executor expands `.N.ht2`.
pub fn hisat2_align(
    sample: &Sample,
    index_prefix: &str,
    splice_sites: &str,
    sam: &str,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    check_threads(cfg)?;
    check_reads(sample)?;
    for path in sample.fastq1.iter().chain(&sample.fastq2) {
        if path.contains(',') {
            return Err(ConfigError::CommaInPath {
                sample: sample.name.clone(),
                path: path.clone(),
            });
        }
    }
    let paired = !sample.fastq2.is_empty();
    let mut args = vec!["-p".to_string(), cfg.threads.to_string()];
    args.extend(st(&[
        "-x",
        index_prefix,
        "--known-splicesite-infile",
        splice_sites,
        "--no-temp-splicesite",
        "--max-intronlen",
    ]));
    args.push(cfg.max_intron.to_string());
    let strand = match (cfg.strandedness, paired) {
        (Strandedness::Unstranded, _) => None,
        (Strandedness::Forward, false) => Some("F"),
        (Strandedness::Forward, true) => Some("FR"),
        (Strandedness::Reverse, false) => Some("R"),
        (Strandedness::Reverse, true) => Some("RF"),
    };
    if let Some(s) = strand {
        args.push("--rna-strandness".into());
        args.push(s.into());
    }
    if paired {
        args.extend([
            "-1".into(),
            sample.fastq1.join(","),
            "-2".into(),
            sample.fastq2.join(","),
        ]);
    } else {
        args.extend(["-U".into(), sample.fastq1.join(",")]);
    }
    args.extend(["-S".into(), sam.into()]);
    let mut inputs: Vec<PathBuf> = sample
        .fastq1
        .iter()
        .chain(&sample.fastq2)
        .map(PathBuf::from)
        .collect();
    inputs.push(splice_sites.into());
    inputs.push(index_prefix.into());
    Ok(ToolCall {
        tool: Tool::Hisat2Align,
        args,
        inputs,
        outputs: vec![sam.into()],
    })
}

/// `featureCounts` on a SAF built from our gene model.
pub fn feature_counts(
    sam: &str,
    saf: &str,
    out: &str,
    paired: bool,
    cfg: &RunConfig,
) -> Result<ToolCall, ConfigError> {
    check_threads(cfg)?;
    let mut args = vec!["-T".to_string(), cfg.threads.to_string()];
    args.extend(st(&["-F", "SAF", "-a", saf]));
    if paired {
        args.extend(st(&["-p", "--countReadPairs"]));
    }
    let s = match cfg.strandedness {
        Strandedness::Unstranded => "0",
        Strandedness::Forward => "1",
        Strandedness::Reverse => "2",
    };
    args.extend(st(&["-s", s, "-o", out, sam]));
    Ok(ToolCall {
        tool: Tool::FeatureCounts,
        args,
        inputs: vec![sam.into(), saf.into()],
        outputs: vec![out.into(), format!("{out}.summary").into()],
    })
}

/// Full argv for native execution, program name first.
pub fn argv_native(call: &ToolCall) -> Vec<String> {
    let mut v = vec![program(call.tool).to_string()];
    v.extend(call.args.iter().cloned());
    v
}
