//! SLURM scripts for `qrnamap run --slurm`.
//!
//! Four files go to `<run_dir>/slurm/`:
//! * `prepare.sbatch`: `qrnamap run ... --step prepare`. Builds the transcripts (quant) or the
//!   SAF and splice sites (align) and the index once, in the run directory on shared storage,
//!   and removes the old outputs of every sample first.
//! * `samples.sbatch`: an array of `ceil(n / per_task)` tasks. Task `i` runs
//!   `qrnamap quant|align ... --task-id i --samples-per-task K` and reuses the prepared index.
//! * `final.sbatch`: `qrnamap run ... --step final` (matrix, then DESeq2 when a contrast is set).
//! * `submit.sh`: submits the three with `sbatch --parsable`; the array depends on prepare and
//!   final on the array, both with `afterok`.
//!
//! Safety contract (ruling R18): every path and argument is written through `sh_quote`. Sample
//! names and FASTQ paths are never written: a task picks its samples by `$SLURM_ARRAY_TASK_ID`
//! and the CLI reads the sample sheet itself. Text with a newline or NUL is refused, not quoted.
//! The scripts never use `BASH_SOURCE` or `dirname`: Rust writes absolute paths in. Temporary
//! files go to node-local `${SCRATCH:?}`, which fails the job when it is not set, unless the plan
//! has a `tmpdir` (`--tmpdir DIR`): then the scripts pass that quoted directory and never use
//! `SCRATCH`.

use std::ffi::OsStr;
use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};

use qrnamap_core::Mode;

/// Default `--time`. UNMEASURED placeholder: no per-sample timing exists yet.
pub const DEFAULT_TIME: &str = "02:00:00";
/// Default `--mem-gb`. Placeholder, small on purpose. Measured on the AMB data set (see
/// `validation/amb/README.md`): HISAT2 used 41,208 KB and featureCounts used 36,144 KB resident
/// memory. The 13 GB that `sacct` reported as MaxRSS was page cache of the staged scratch files.
/// A job with `--mem 4G` finished with identical counts. NOT measured for genomes near 50 Mb or
/// for 20 million read pairs per sample: raise `--mem-gb` for those.
pub const DEFAULT_MEM_GB: u32 = 8;

/// Message of the `${SCRATCH:?...}` check. In the generated script `$SCRATCH` is written as
/// `\$SCRATCH` so that it stays literal in the parameter expansion.
const SCRATCH_MSG: &str = "SCRATCH is not set: this job needs node-local scratch. Submit from a login shell where \\$SCRATCH is set, or regenerate the scripts with --tmpdir DIR.";

/// Everything needed to write the scripts. It holds the sample COUNT only, never the names.
#[derive(Debug, Clone)]
pub struct Plan {
    pub n_samples: usize,
    /// Samples per array task (`--samples-per-task`).
    pub per_task: usize,
    pub cpus: u32,
    pub mem_gb: u32,
    /// SLURM time limit, for example `02:00:00` or `1-00:00:00`.
    pub time: String,
    pub partition: Option<String>,
    pub account: Option<String>,
    /// Absolute run directory (`--out`) on shared storage.
    pub run_dir: PathBuf,
    /// Absolute directory the jobs work in (relative FASTQ paths in the sheet resolve here).
    pub work_dir: PathBuf,
    /// Absolute path of the `qrnamap` binary.
    pub qrnamap_bin: PathBuf,
    pub mode: Mode,
    /// Arguments after `qrnamap run` for the prepare and final jobs (without `--step`,
    /// `--tmpdir` and `--slurm`).
    pub run_args: Vec<String>,
    /// Arguments after `qrnamap quant|align` for the array tasks (without `--task-id`,
    /// `--samples-per-task` and `--tmpdir`).
    pub task_args: Vec<String>,
    /// Directories written into `export PATH=...:"$PATH"` at the top of every .sbatch, in
    /// order (ruling R28). Absolute.
    pub tool_dirs: Vec<PathBuf>,
    /// `--tmpdir DIR` of `run --slurm`: absolute, baked into the scripts as `--tmpdir <DIR>`.
    /// The qrnamap binary creates it (create_dir_all) when its parent is writable; it should be on
    /// storage every compute node can reach. When `None`, the jobs
    /// use node-local `${SCRATCH:?}`.
    pub tmpdir: Option<PathBuf>,
    /// Tool names `prepare.sbatch` checks with `command -v` before it starts qrnamap.
    pub tools: Vec<String>,
}

/// The generated script texts.
#[derive(Debug, Clone)]
pub struct SlurmScripts {
    pub prepare: String,
    pub samples: String,
    pub final_: String,
    pub submit: String,
    pub n_tasks: usize,
}

#[derive(Debug)]
pub enum SlurmError {
    /// The sample sheet has no samples, so there is no array.
    NoSamples,
    /// A value that goes into an `#SBATCH` line or a count is not allowed.
    BadValue { what: &'static str, value: String },
    /// Text that would go into a script holds a newline or NUL. It is refused, not quoted.
    NewlineOrNul(String),
    NonUtf8Path(PathBuf),
    /// A path that must be absolute is relative.
    RelativePath(PathBuf),
    /// The run directory holds `%` or `\`, which SLURM reads as a pattern in `--output`.
    LogPathPattern(PathBuf),
    /// Tools the jobs need that are neither on PATH nor in a `--tool-path` directory.
    MissingTools(Vec<String>),
    /// A `--tool-path` directory (or a tool directory) is not usable.
    ToolPath { path: PathBuf, detail: String },
    Io { path: PathBuf, source: std::io::Error },
}

impl fmt::Display for SlurmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SlurmError::NoSamples => write!(f, "the sample sheet has no samples, so there is no SLURM array to write"),
            SlurmError::BadValue { what, value } => write!(f, "{what}: value {value:?} is not allowed"),
            SlurmError::NewlineOrNul(s) => {
                write!(f, "{s:?} contains a newline or NUL character and cannot be written into a shell script")
            }
            SlurmError::NonUtf8Path(p) => write!(f, "path {} is not valid UTF-8", p.display()),
            SlurmError::RelativePath(p) => write!(f, "path {} must be absolute", p.display()),
            SlurmError::MissingTools(t) => write!(
                f,
                "required tools not found on PATH or in a --tool-path directory: {}. Put them on PATH or pass --tool-path DIR",
                t.join(", ")
            ),
            SlurmError::ToolPath { path, detail } => write!(f, "--tool-path {}: {detail}", path.display()),
            SlurmError::LogPathPattern(p) => write!(
                f,
                "run directory {} contains '%' or '\\', which SLURM reads as a file name pattern; choose another --out",
                p.display()
            ),
            SlurmError::Io { path, source } => write!(f, "{}: {source}", path.display()),
        }
    }
}

impl std::error::Error for SlurmError {}

/// Quote `s` as one shell word: wrap it in single quotes and write each `'` as `'\''`. Inside
/// single quotes the shell reads every other byte as itself (`$`, backticks, `\`, `"`, `;`, `*`).
/// A newline or NUL is refused with `NewlineOrNul`.
pub fn sh_quote(s: &str) -> Result<String, SlurmError> {
    if s.contains(['\n', '\0']) {
        return Err(SlurmError::NewlineOrNul(s.to_string()));
    }
    Ok(format!("'{}'", s.replace('\'', "'\\''")))
}

fn path_str(p: &Path) -> Result<&str, SlurmError> {
    if !p.is_absolute() {
        return Err(SlurmError::RelativePath(p.to_path_buf()));
    }
    p.to_str().ok_or_else(|| SlurmError::NonUtf8Path(p.to_path_buf()))
}

fn quote_path(p: &Path) -> Result<String, SlurmError> {
    sh_quote(path_str(p)?)
}

fn quote_all(words: &[String]) -> Result<String, SlurmError> {
    Ok(words.iter().map(|w| sh_quote(w)).collect::<Result<Vec<_>, _>>()?.join(" "))
}

fn check_sbatch_value(what: &'static str, v: &str, ok: impl Fn(char) -> bool) -> Result<(), SlurmError> {
    if v.is_empty() || !v.chars().all(ok) {
        return Err(SlurmError::BadValue { what, value: v.to_string() });
    }
    Ok(())
}

/// True when the path can go into an `#SBATCH` line without quoting.
fn plain_path(s: &str) -> bool {
    s.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '/' | '.' | '_' | '-' | '+'))
}

fn n_tasks(plan: &Plan) -> Result<usize, SlurmError> {
    if plan.n_samples == 0 {
        return Err(SlurmError::NoSamples);
    }
    if plan.per_task == 0 {
        return Err(SlurmError::BadValue { what: "--samples-per-task", value: "0".into() });
    }
    Ok(plan.n_samples.div_ceil(plan.per_task))
}

struct Checked<'a> {
    plan: &'a Plan,
    n_tasks: usize,
    logs: String,
    /// `export PATH=<tool dirs>:"$PATH"` line, or empty.
    export: String,
    /// `command -v` check of the required tools (prepare.sbatch), or empty.
    tool_check: String,
    /// Quoted `--tmpdir` directory, or `None` to use `${SCRATCH:?}`.
    tmpdir: Option<String>,
}

fn check(plan: &Plan) -> Result<Checked<'_>, SlurmError> {
    let n_tasks = n_tasks(plan)?;
    if plan.cpus == 0 {
        return Err(SlurmError::BadValue { what: "--cpus", value: "0".into() });
    }
    if plan.mem_gb == 0 {
        return Err(SlurmError::BadValue { what: "--mem-gb", value: "0".into() });
    }
    check_sbatch_value("--time", &plan.time, |c| c.is_ascii_digit() || matches!(c, ':' | '-'))?;
    let name_char = |c: char| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | ',' | '-');
    if let Some(p) = &plan.partition {
        check_sbatch_value("--partition", p, name_char)?;
    }
    if let Some(a) = &plan.account {
        check_sbatch_value("--account", a, name_char)?;
    }
    let run_dir = path_str(&plan.run_dir)?;
    sh_quote(run_dir)?;
    if run_dir.contains(['%', '\\']) {
        return Err(SlurmError::LogPathPattern(plan.run_dir.clone()));
    }
    path_str(&plan.work_dir)?;
    path_str(&plan.qrnamap_bin)?;
    let mut dirs = Vec::new();
    for d in &plan.tool_dirs {
        let s = path_str(d)?;
        if s.contains(':') {
            return Err(SlurmError::ToolPath { path: d.clone(), detail: "contains ':', which separates PATH entries".into() });
        }
        dirs.push(sh_quote(s)?);
    }
    let export = if dirs.is_empty() { String::new() } else { format!("export PATH={}:\"$PATH\"\n", dirs.join(":")) };
    let tool_check = if plan.tools.is_empty() {
        String::new()
    } else {
        let names = plan.tools.iter().map(|t| sh_quote(t)).collect::<Result<Vec<_>, _>>()?.join(" ");
        format!(
            "missing=\"\"\nfor tool in {names}; do\n  command -v \"$tool\" >/dev/null 2>&1 || missing=\"$missing $tool\"\ndone\nif [ -n \"$missing\" ]; then\n  echo \"qrnamap: required tools not found on PATH:$missing\" >&2\n  exit 1\nfi\n"
        )
    };
    let tmpdir = plan.tmpdir.as_deref().map(quote_path).transpose()?;
    let logs = plan.run_dir.join("slurm").join("logs");
    Ok(Checked { plan, n_tasks, logs: path_str(&logs)?.to_string(), export, tool_check, tmpdir })
}

impl Checked<'_> {
    fn header(&self, job: &str, array: bool, what: &str) -> String {
        let p = self.plan;
        let mut h = String::from("#!/bin/bash\n");
        h.push_str(&format!("#SBATCH --job-name=qrnamap-{job}\n"));
        if array {
            h.push_str(&format!("#SBATCH --array=0-{}\n", self.n_tasks - 1));
        }
        h.push_str(&format!("#SBATCH --cpus-per-task={}\n", p.cpus));
        h.push_str(&format!("#SBATCH --mem={}G\n", p.mem_gb));
        h.push_str(&format!("#SBATCH --time={}\n", p.time));
        if let Some(x) = &p.partition {
            h.push_str(&format!("#SBATCH --partition={x}\n"));
        }
        if let Some(x) = &p.account {
            h.push_str(&format!("#SBATCH --account={x}\n"));
        }
        // submit.sh also passes --output (quoted) on the command line, which wins over this
        // line. The line is written only when the path needs no quoting.
        if plain_path(&self.logs) {
            h.push_str(&format!("#SBATCH --output={}/{}\n", self.logs, log_pattern(array)));
        }
        h.push_str("#\n");
        h.push_str(&format!("# Generated by `qrnamap run --slurm`. {what}\n"));
        h.push_str("# The --cpus-per-task, --mem and --time values are UNMEASURED placeholders (no per-sample\n");
        h.push_str("# timing exists yet). Set them from a measured run with --cpus, --mem-gb and --time.\n");
        h.push_str("# The --mem default is small: AMB data used about 40 MiB (41,208 KB) resident (not measured for large genomes).\n");
        h.push_str("# Sample names and FASTQ paths are not in this script: qrnamap reads the sample sheet.\n");
        h.push_str("# The tool directories below were found when the scripts were written (ruling R28), so the\n");
        h.push_str("# jobs do not depend on the PATH of the shell that runs submit.sh.\n");
        h.push_str("set -euo pipefail\n");
        h.push_str(&self.export);
        if self.tmpdir.is_none() {
            h.push_str(&format!(": \"${{SCRATCH:?{SCRATCH_MSG}}}\"\n"));
        }
        h
    }

    fn command(&self, sub: &str, args: &[String], extra: &str) -> Result<String, SlurmError> {
        let p = self.plan;
        let tmp = self.tmpdir.as_deref().unwrap_or("\"${SCRATCH:?}\"");
        Ok(format!(
            "cd -- {}\nexec {} {} {}{extra} --tmpdir {tmp}\n",
            quote_path(&p.work_dir)?,
            quote_path(&p.qrnamap_bin)?,
            sh_quote(sub)?,
            quote_all(args)?,
        ))
    }
}

fn log_pattern(array: bool) -> &'static str {
    if array {
        "%x_%A_%a.out"
    } else {
        "%x_%j.out"
    }
}

fn mode_word(m: Mode) -> &'static str {
    match m {
        Mode::Quant => "quant",
        Mode::Align => "align",
    }
}

/// Write the four script texts. Nothing is written to disk.
pub fn slurm_scripts(plan: &Plan) -> Result<SlurmScripts, SlurmError> {
    let c = check(plan)?;
    let k = plan.per_task;
    let prepare = c.header(
        "prepare",
        false,
        "Builds the reference files and the index once in the run directory, and removes the old outputs of all samples.",
    ) + &c.tool_check
        + &c.command("run", &plan.run_args, " '--step' 'prepare'")?;
    let samples = c.header(
        "samples",
        true,
        &format!("Array task i processes samples [i*{k}, (i+1)*{k}) of the sample sheet, in sheet order."),
    ) + ": \"${SLURM_ARRAY_TASK_ID:?this script runs as a SLURM array task}\"\n"
        + &c.command(
            mode_word(plan.mode),
            &plan.task_args,
            &format!(" '--task-id' \"$SLURM_ARRAY_TASK_ID\" '--samples-per-task' {}", sh_quote(&k.to_string())?),
        )?;
    let final_ = c.header("final", false, "Gene matrix, then DESeq2 when a contrast was given.")
        + &c.command("run", &plan.run_args, " '--step' 'final'")?;

    let sdir = plan.run_dir.join("slurm");
    let q = |name: &str| quote_path(&sdir.join(name));
    let out = |array: bool| sh_quote(&format!("{}/{}", c.logs, log_pattern(array)));
    let mut submit = String::from("#!/bin/bash\n");
    submit.push_str("# Generated by `qrnamap run --slurm`. Submits three jobs with sbatch:\n");
    submit.push_str("#   prepare: reference files and index, once; removes old outputs of all samples\n");
    submit.push_str(&format!("#   samples: array of {} task(s), {k} sample(s) per task; starts after prepare succeeds\n", c.n_tasks));
    submit.push_str("#   final:   gene matrix (and DESeq2 when a contrast was given); starts after every task succeeds\n");
    submit.push_str("# All paths are written in by qrnamap. This script does not look up its own location.\n");
    submit.push_str("# The samples and final jobs use --kill-on-invalid-dep=yes: when a job they depend on fails,\n");
    submit.push_str("# SLURM cancels them instead of leaving them pending.\n");
    submit.push_str("#\n");
    submit.push_str("# Recovery when some array tasks fail (the final job is then cancelled, or still pending):\n");
    submit.push_str("#   1. If the final job is still pending, cancel it (its id is printed by this script):\n");
    submit.push_str("#        scancel FINAL_JOB_ID\n");
    submit.push_str("#   2. Fix the cause, then run only the failed task ids again (TASK_IDS such as 3,5):\n");
    submit.push_str(&format!("#        sbatch --array=TASK_IDS --output={} {}\n", out(true)?, q("samples.sbatch")?));
    submit.push_str("#   3. Submit the final job again, after the new array job (NEW_ARRAY_JOB_ID) succeeds:\n");
    submit.push_str(&format!(
        "#        sbatch --kill-on-invalid-dep=yes --dependency=afterok:NEW_ARRAY_JOB_ID --output={} {}\n",
        out(false)?,
        q("final.sbatch")?
    ));
    submit.push_str("set -euo pipefail\n");
    submit.push_str(&format!("prepare_id=$(sbatch --parsable --output={} {})\n", out(false)?, q("prepare.sbatch")?));
    submit.push_str("prepare_id=${prepare_id%%;*}\n");
    submit.push_str(&format!(
        "samples_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:\"$prepare_id\" --output={} {})\n",
        out(true)?,
        q("samples.sbatch")?
    ));
    submit.push_str("samples_id=${samples_id%%;*}\n");
    submit.push_str(&format!(
        "final_id=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:\"$samples_id\" --output={} {})\n",
        out(false)?,
        q("final.sbatch")?
    ));
    submit.push_str("final_id=${final_id%%;*}\n");
    submit.push_str("echo \"prepare job: $prepare_id\"\n");
    submit.push_str(&format!("echo \"samples array job: $samples_id ({} task(s))\"\n", c.n_tasks));
    submit.push_str("echo \"final job: $final_id\"\n");
    submit.push_str(&format!("echo \"logs: \"{}\n", sh_quote(&c.logs)?));
    Ok(SlurmScripts { prepare, samples, final_, submit, n_tasks: c.n_tasks })
}

fn is_executable_file(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    fs::metadata(p).map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0).unwrap_or(false)
}

/// The tool directories to bake into the scripts (ruling R28). `explicit` (`--tool-path`) come
/// first: each must exist and is canonicalized. Then, for each tool not found in them, the
/// first directory of `path_env` (the current PATH) that holds it, made absolute. Duplicates
/// are dropped; the order is stable. A tool found nowhere is listed in `MissingTools`.
pub fn resolve_tool_dirs(tools: &[&str], explicit: &[PathBuf], path_env: Option<&OsStr>) -> Result<Vec<PathBuf>, SlurmError> {
    let mut dirs: Vec<PathBuf> = Vec::new();
    for d in explicit {
        let canon = fs::canonicalize(d).map_err(|e| SlurmError::ToolPath { path: d.clone(), detail: e.to_string() })?;
        if !canon.is_dir() {
            return Err(SlurmError::ToolPath { path: d.clone(), detail: "is not a directory".into() });
        }
        if !dirs.contains(&canon) {
            dirs.push(canon);
        }
    }
    let path_dirs: Vec<PathBuf> = path_env.map(|p| std::env::split_paths(p).collect()).unwrap_or_default();
    let mut missing = Vec::new();
    for t in tools {
        if dirs.iter().any(|d| is_executable_file(&d.join(t))) {
            continue;
        }
        match path_dirs.iter().find(|d| !d.as_os_str().is_empty() && is_executable_file(&d.join(t))) {
            Some(d) => {
                let abs = std::path::absolute(d).map_err(|e| SlurmError::ToolPath { path: d.clone(), detail: e.to_string() })?;
                if !dirs.contains(&abs) {
                    dirs.push(abs);
                }
            }
            None => missing.push(t.to_string()),
        }
    }
    if !missing.is_empty() {
        return Err(SlurmError::MissingTools(missing));
    }
    Ok(dirs)
}

/// The array script (`samples.sbatch`) text.
pub fn slurm_array_script(plan: &Plan) -> Result<String, SlurmError> {
    Ok(slurm_scripts(plan)?.samples)
}

/// Write the scripts into `<run_dir>/slurm/` (and create `<run_dir>/slurm/logs/`). Returns
/// the path of `submit.sh`.
pub fn write_slurm_scripts(plan: &Plan) -> Result<PathBuf, SlurmError> {
    use std::os::unix::fs::PermissionsExt;
    let s = slurm_scripts(plan)?;
    let sdir = plan.run_dir.join("slurm");
    let logs = sdir.join("logs");
    fs::create_dir_all(&logs).map_err(|source| SlurmError::Io { path: logs.clone(), source })?;
    let submit = sdir.join("submit.sh");
    for (name, text) in [("prepare.sbatch", &s.prepare), ("samples.sbatch", &s.samples), ("final.sbatch", &s.final_), ("submit.sh", &s.submit)] {
        let p = sdir.join(name);
        fs::write(&p, text).map_err(|source| SlurmError::Io { path: p.clone(), source })?;
        fs::set_permissions(&p, fs::Permissions::from_mode(0o755)).map_err(|source| SlurmError::Io { path: p.clone(), source })?;
    }
    Ok(submit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sh_quote_unit_cases() {
        for (input, want) in [
            ("", "''"),
            ("a b", "'a b'"),
            ("$HOME", "'$HOME'"),
            ("`x`", "'`x`'"),
            ("'", "''\\'''"),
            ("\"", "'\"'"),
            ("\\", "'\\'"),
            (";", "';'"),
            ("&", "'&'"),
            ("*", "'*'"),
        ] {
            assert_eq!(sh_quote(input).unwrap(), want, "{input:?}");
        }
        assert!(matches!(sh_quote("a\nb"), Err(SlurmError::NewlineOrNul(_))));
        assert!(matches!(sh_quote("\0"), Err(SlurmError::NewlineOrNul(_))));
    }

    #[test]
    fn task_count_rounds_up() {
        let mut p = Plan {
            n_samples: 7,
            per_task: 3,
            cpus: 1,
            mem_gb: 1,
            time: "1:00:00".into(),
            partition: None,
            account: None,
            run_dir: "/r".into(),
            work_dir: "/w".into(),
            qrnamap_bin: "/b".into(),
            mode: Mode::Align,
            run_args: vec![],
            task_args: vec![],
            tool_dirs: vec![],
            tmpdir: None,
            tools: vec![],
        };
        assert_eq!(n_tasks(&p).unwrap(), 3);
        p.n_samples = 6;
        assert_eq!(n_tasks(&p).unwrap(), 2);
        p.n_samples = 0;
        assert!(matches!(n_tasks(&p), Err(SlurmError::NoSamples)));
    }
}
