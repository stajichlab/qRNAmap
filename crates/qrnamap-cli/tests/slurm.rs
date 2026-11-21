//! SLURM script generation. No real SLURM job is submitted: `submit.sh` runs against a stub
//! `sbatch` on a temporary PATH that records its arguments.

mod common;

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

use qrnamap_cli::slurm::{sh_quote, DEFAULT_MEM_GB, slurm_array_script, slurm_scripts, write_slurm_scripts, Plan, SlurmError};
use qrnamap_core::{Mode, Quantifier};

fn plan(n: usize, k: usize, run_dir: &Path) -> Plan {
    Plan {
        n_samples: n,
        per_task: k,
        cpus: 4,
        mem_gb: DEFAULT_MEM_GB,
        time: "02:00:00".into(),
        partition: Some("batch".into()),
        account: None,
        run_dir: run_dir.to_path_buf(),
        work_dir: PathBuf::from("/bigdata/lab/work"),
        qrnamap_bin: PathBuf::from("/opt/qrnamap/bin/qrnamap"),
        mode: Mode::Quant,
        quantifier: Some(Quantifier::Kallisto),
        run_args: vec!["--samples".into(), "/bigdata/lab/sheet.tsv".into(), "--out".into(), run_dir.to_str().unwrap().into()],
        task_args: vec!["--samples".into(), "/bigdata/lab/sheet.tsv".into(), "--out".into(), run_dir.to_str().unwrap().into()],
        tool_dirs: vec![PathBuf::from("/opt/tools/bin")],
        tools: vec!["kallisto".into()],
        tmpdir: None,
    }
}

fn all_scripts(p: &Plan) -> Vec<(String, String)> {
    let s = slurm_scripts(p).unwrap();
    vec![
        ("prepare.sbatch".into(), s.prepare),
        ("samples.sbatch".into(), s.samples),
        ("final.sbatch".into(), s.final_),
        ("submit.sh".into(), s.submit),
    ]
}

#[test]
fn slurm_script_uses_scratch_with_required_check() {
    let p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    for (name, text) in all_scripts(&p) {
        if name.ends_with(".sbatch") {
            assert!(text.contains("\"${SCRATCH:?}\""), "{name}:\n{text}");
            assert!(text.contains("--tmpdir \"${SCRATCH:?}\""), "{name}:\n{text}");
            assert!(text.contains("set -euo pipefail"), "{name}");
        }
    }
}

#[test]
fn slurm_script_has_no_bash_source() {
    let p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    for (name, text) in all_scripts(&p) {
        assert!(!text.contains("BASH_SOURCE"), "{name}:\n{text}");
        assert!(!text.contains("dirname"), "{name}:\n{text}");
    }
}

#[test]
fn slurm_script_has_no_hardcoded_scratch_path() {
    let p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    for (name, text) in all_scripts(&p) {
        assert!(!text.contains("/scratch/"), "{name}:\n{text}");
    }
}

#[test]
fn slurm_resources_are_marked_unmeasured() {
    let p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    for (name, text) in all_scripts(&p) {
        if name.ends_with(".sbatch") {
            assert!(text.contains("UNMEASURED"), "{name}");
            for want in ["#SBATCH --cpus-per-task=4", "#SBATCH --mem=8G", "#SBATCH --time=02:00:00", "#SBATCH --partition=batch", "#SBATCH --job-name="] {
                assert!(text.contains(want), "{name} lacks {want}:\n{text}");
            }
            assert!(!text.contains("--account"), "{name}");
        }
    }
}

#[test]
fn slurm_array_range_matches_task_count() {
    let d = Path::new("/bigdata/lab/run1");
    for (n, k, want) in [(7, 3, "--array=0-2"), (6, 3, "--array=0-1"), (1, 1, "--array=0-0"), (1, 3, "--array=0-0"), (48, 1, "--array=0-47")] {
        let s = slurm_array_script(&plan(n, k, d)).unwrap();
        assert!(s.contains(&format!("#SBATCH {want}\n")), "{n}/{k}:\n{s}");
        assert!(s.contains(&format!("'--samples-per-task' '{k}'")), "{s}");
        assert!(s.contains("'--task-id' \"$SLURM_ARRAY_TASK_ID\""), "{s}");
    }
    assert!(matches!(slurm_array_script(&plan(0, 3, d)), Err(SlurmError::NoSamples)));
    assert!(matches!(slurm_array_script(&plan(3, 0, d)), Err(SlurmError::BadValue { .. })));
}

#[test]
fn slurm_rejects_unsafe_sbatch_values() {
    let d = Path::new("/bigdata/lab/run1");
    let mut p = plan(3, 1, d);
    p.partition = Some("batch\n#SBATCH --x".into());
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::BadValue { .. })));
    let mut p = plan(3, 1, d);
    p.account = Some("lab acct".into());
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::BadValue { .. })));
    let mut p = plan(3, 1, d);
    p.time = "1:00:00; rm".into();
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::BadValue { .. })));
    let mut p = plan(3, 1, d);
    p.run_dir = PathBuf::from("relative/run");
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::RelativePath(_))));
    let mut p = plan(3, 1, d);
    p.run_dir = PathBuf::from("/bigdata/run%x");
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::LogPathPattern(_))));
    let mut p = plan(3, 1, d);
    p.run_dir = PathBuf::from("/bigdata/new\nline");
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::NewlineOrNul(_))));
}

#[test]
fn slurm_index_job_precedes_array() {
    let p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    let s = slurm_scripts(&p).unwrap();
    let sub = &s.submit;
    assert!(sub.starts_with("#!/bin/bash\n"), "{sub}");
    assert!(sub.contains("set -euo pipefail"));
    let line = |key: &str| sub.lines().find(|l| l.starts_with(key)).unwrap_or_else(|| panic!("{key} in\n{sub}")).to_string();
    let prep = line("prepare_id=$(sbatch --parsable");
    let arr = line("samples_id=$(sbatch --parsable");
    let fin = line("final_id=$(sbatch --parsable");
    assert!(!prep.contains("--dependency"), "{prep}");
    assert!(arr.contains("--dependency=afterok:\"$prepare_id\""), "{arr}");
    assert!(fin.contains("--dependency=afterok:\"$samples_id\""), "{fin}");
    let pos = |l: &str| sub.find(l).unwrap();
    assert!(pos(&prep) < pos(&arr) && pos(&arr) < pos(&fin));
    // The prepare job builds the index; the array reuses it.
    assert!(s.prepare.contains("'run'") && s.prepare.contains("'--step' 'prepare'"), "{}", s.prepare);
    assert!(s.final_.contains("'--step' 'final'"), "{}", s.final_);
    assert!(s.samples.contains("'/opt/qrnamap/bin/qrnamap' 'quant'"), "{}", s.samples);
}

// ---------------------------------------------------------------------------------------------
// sh_quote

const HOSTILE: [&str; 12] = [
    "with space",
    "$HOME",
    "`id`",
    "it's",
    "say \"hi\"",
    "back\\slash",
    "a;b",
    "a&b",
    "*",
    "",
    "$(touch /nonexistent/x)",
    "'''",
];

#[test]
fn sh_quote_round_trips_through_sh() {
    for s in HOSTILE {
        let q = sh_quote(s).unwrap();
        assert!(q.starts_with('\'') && q.ends_with('\''), "{q}");
        let out = Command::new("/bin/sh").arg("-c").arg(format!("printf %s {q}")).output().unwrap();
        assert!(out.status.success(), "{s:?}");
        assert_eq!(String::from_utf8(out.stdout).unwrap(), s, "quoted as {q}");
    }
}

#[test]
fn sh_quote_rejects_newline_and_nul() {
    assert!(matches!(sh_quote("a\nb"), Err(SlurmError::NewlineOrNul(_))));
    assert!(matches!(sh_quote("a\0b"), Err(SlurmError::NewlineOrNul(_))));
}

#[test]
fn sh_quote_exact_forms() {
    assert_eq!(sh_quote("").unwrap(), "''");
    assert_eq!(sh_quote("abc").unwrap(), "'abc'");
    assert_eq!(sh_quote("it's").unwrap(), "'it'\\''s'");
}

// ---------------------------------------------------------------------------------------------
// Hostile run directory, executed with stubs

fn write_exe(path: &Path, body: &str) {
    fs::write(path, body).unwrap();
    fs::set_permissions(path, fs::Permissions::from_mode(0o755)).unwrap();
}

/// Arguments of each stub call, split on the `CALL` marker line.
fn calls(log: &Path) -> Vec<Vec<String>> {
    let text = fs::read_to_string(log).unwrap();
    let mut out: Vec<Vec<String>> = Vec::new();
    for l in text.lines() {
        if l == "CALL" {
            out.push(Vec::new());
        } else {
            out.last_mut().unwrap().push(l.to_string());
        }
    }
    out
}

#[test]
fn slurm_script_with_hostile_run_dir_round_trips() {
    let tmp = tempfile::tempdir().unwrap();
    let run_dir = tmp.path().join("run dir $HOME it's `echo pwned`");
    let work_dir = tmp.path().join("work $USER");
    fs::create_dir_all(&work_dir).unwrap();
    let stub = tmp.path().join("stub");
    fs::create_dir_all(&stub).unwrap();
    let log = tmp.path().join("calls.log");
    // Stub sbatch: record the arguments, print a job id 1001, 1002, ... in --parsable form.
    write_exe(
        &stub.join("sbatch"),
        "#!/bin/sh\nn=$(grep -c '^CALL$' \"$STUB_LOG\" 2>/dev/null || true)\nn=$((${n:-0}+1))\n{ echo CALL; printf '%s\\n' \"$@\"; } >> \"$STUB_LOG\"\necho \"100$n;cluster\"\n",
    );
    // Stub qrnamap: record the arguments and the working directory.
    let fake_bin = stub.join("qrnamap stub");
    write_exe(&fake_bin, "#!/bin/sh\n{ echo CALL; pwd; printf '%s\\n' \"$@\"; } >> \"$STUB_LOG\"\n");

    let mut p = plan(7, 3, &run_dir);
    p.work_dir = work_dir.clone();
    p.qrnamap_bin = fake_bin.clone();
    let rd = run_dir.to_str().unwrap().to_string();
    p.run_args = vec!["--out".into(), rd.clone()];
    p.task_args = vec!["--out".into(), rd.clone()];
    let submit = write_slurm_scripts(&p).unwrap();
    let sdir = run_dir.join("slurm");
    assert_eq!(submit, sdir.join("submit.sh"));
    assert!(sdir.join("logs").is_dir());
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch", "submit.sh"] {
        let o = Command::new("bash").arg("-n").arg(sdir.join(f)).output().unwrap();
        assert!(o.status.success(), "bash -n {f}: {}", String::from_utf8_lossy(&o.stderr));
    }

    let path = format!("{}:{}", stub.display(), std::env::var("PATH").unwrap());
    let o = Command::new("bash").arg(&submit).env("PATH", &path).env("STUB_LOG", &log).env("HOME", "/home/EXPANDED").output().unwrap();
    let stdout = String::from_utf8_lossy(&o.stdout).into_owned();
    assert!(o.status.success(), "{stdout}\n{}", String::from_utf8_lossy(&o.stderr));
    assert!(stdout.contains("1001") && stdout.contains("1002") && stdout.contains("1003"), "{stdout}");
    let c = calls(&log);
    assert_eq!(c.len(), 3, "{c:?}");
    let logs = format!("{rd}/slurm/logs");
    assert_eq!(c[0], vec!["--parsable".to_string(), format!("--output={logs}/%x_%j.out"), format!("{rd}/slurm/prepare.sbatch")]);
    assert_eq!(
        c[1],
        vec![
            "--parsable".to_string(),
            "--kill-on-invalid-dep=yes".into(),
            "--dependency=afterok:1001".into(),
            format!("--output={logs}/%x_%A_%a.out"),
            format!("{rd}/slurm/samples.sbatch")
        ]
    );
    assert_eq!(
        c[2],
        vec![
            "--parsable".to_string(),
            "--kill-on-invalid-dep=yes".into(),
            "--dependency=afterok:1002".into(),
            format!("--output={logs}/%x_%j.out"),
            format!("{rd}/slurm/final.sbatch")
        ]
    );

    // Run the array script itself as task 1 with the stub binary: the hostile path reaches the
    // program as one argument, unexpanded, and the task id and scratch are passed through.
    fs::remove_file(&log).unwrap();
    let scratch = tmp.path().join("node scratch");
    fs::create_dir_all(&scratch).unwrap();
    let o = Command::new("bash")
        .arg(sdir.join("samples.sbatch"))
        .env("STUB_LOG", &log)
        .env("SCRATCH", &scratch)
        .env("SLURM_ARRAY_TASK_ID", "1")
        .env("HOME", "/home/EXPANDED")
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let c = calls(&log);
    assert_eq!(c.len(), 1);
    assert_eq!(c[0][0], work_dir.to_str().unwrap(), "the script works in the planned directory");
    assert_eq!(
        c[0][1..],
        ["quant", "--out", &rd, "--task-id", "1", "--samples-per-task", "3", "--tmpdir", scratch.to_str().unwrap()]
    );

    // Without $SCRATCH the job fails loudly before the program runs.
    fs::remove_file(&log).unwrap();
    let o = Command::new("bash")
        .arg(sdir.join("samples.sbatch"))
        .env("STUB_LOG", &log)
        .env_remove("SCRATCH")
        .env("SLURM_ARRAY_TASK_ID", "1")
        .output()
        .unwrap();
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("SCRATCH"));
    assert!(!log.exists(), "the program must not run without $SCRATCH");
    // No log pattern is written into a #SBATCH line for a path with spaces or quotes.
    let samples = fs::read_to_string(sdir.join("samples.sbatch")).unwrap();
    assert!(!samples.contains("#SBATCH --output"), "{samples}");
}

// ---------------------------------------------------------------------------------------------
// `qrnamap run --slurm` through the binary

/// Run the binary with `PATH` set to `empty_path` (no tools). With `stubs`, `--tool-path stubs`
/// is added so script generation finds the (stub) tools there (ruling R28).
fn qrnamap_no_tools(args: &[&str], empty_path: &Path) -> std::process::Output {
    let stubs = empty_path.parent().unwrap().join("stub_tools");
    stub_tools(&stubs, ALL_TOOLS, &empty_path.parent().unwrap().join("stub_tools.log"));
    Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args(args)
        .args(["--tool-path", stubs.to_str().unwrap()])
        .env("PATH", empty_path)
        .output()
        .unwrap()
}

/// Every tool name the generated scripts can require.
const ALL_TOOLS: &[&str] = &[
    "kallisto", "salmon", "hisat2-build", "hisat2-build-s", "hisat2-align-s", "python",
    "featureCounts", "Rscript",
];

/// Stub executables that append `<name>` to `log` when started.
fn stub_tools(dir: &Path, names: &[&str], log: &Path) {
    fs::create_dir_all(dir).unwrap();
    for n in names {
        write_exe(&dir.join(n), &format!("#!/bin/sh\necho {n} >> {}\n", sh_quote(log.to_str().unwrap()).unwrap()));
    }
}

/// A sheet with distinctive sample names and FASTQ paths (the FASTQ files need not exist:
/// `--slurm` reads no reads).
fn distinct_sheet(d: &Path) -> (PathBuf, Vec<String>) {
    let mut t = String::from("sample\tfastq1\tfastq2\tcondition\n");
    let mut secrets = Vec::new();
    for i in 0..7 {
        let name = format!("ZZsampleUnique{i}");
        let r1 = format!("/data/QQfastqUnique{i}_R1.fq.gz");
        let r2 = format!("/data/QQfastqUnique{i}_R2.fq.gz");
        let cond = if i < 4 { "ctl" } else { "trt" };
        t.push_str(&format!("{name}\t{r1}\t{r2}\t{cond}\n"));
        secrets.extend([name, r1, r2]);
    }
    let p = d.join("sheet.tsv");
    fs::write(&p, t).unwrap();
    (p, secrets)
}

#[test]
fn slurm_flag_runs_no_tool() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let empty = d.join("empty_path");
    fs::create_dir_all(&empty).unwrap();
    let (sheet, _) = distinct_sheet(d);
    let out = d.join("run");
    let fx = common::mini_fixture();
    let o = qrnamap_no_tools(
        &[
            "run", "--samples", sheet.to_str().unwrap(), "--reference", fx.join("genome.fa").to_str().unwrap(),
            "--gff", fx.join("genes.gff3").to_str().unwrap(), "--mode", "align", "--out", out.to_str().unwrap(),
            "--slurm", "--samples-per-task", "3", "--contrast", "condition,trt,ctl", "--partition", "short",
        ],
        &empty,
    );
    let stdout = String::from_utf8_lossy(&o.stdout).into_owned();
    assert!(o.status.success(), "{stdout}\n{}", String::from_utf8_lossy(&o.stderr));
    let sdir = out.join("slurm");
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch", "submit.sh"] {
        assert!(sdir.join(f).is_file(), "{f}");
    }
    assert!(sdir.join("logs").is_dir());
    assert!(stdout.contains(&format!("bash {}", sh_quote(sdir.join("submit.sh").to_str().unwrap()).unwrap())), "{stdout}");
    // Only the scripts were written: no index, no counts, no marker.
    let names: Vec<String> = fs::read_dir(&out).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().into_owned()).collect();
    assert_eq!(names, ["slurm"]);
    assert!(!d.join("stub_tools.log").exists(), "no tool may be started by --slurm");
    let samples = fs::read_to_string(sdir.join("samples.sbatch")).unwrap();
    assert!(samples.contains("#SBATCH --array=0-2\n") && samples.contains("'align'"), "{samples}");
    assert!(samples.contains("#SBATCH --partition=short\n"));
    let fin = fs::read_to_string(sdir.join("final.sbatch")).unwrap();
    assert!(fin.contains("'--contrast' 'condition,trt,ctl'"), "{fin}");
}

#[test]
fn slurm_scripts_never_contain_sample_names_or_fastq_paths() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let empty = d.join("empty_path");
    fs::create_dir_all(&empty).unwrap();
    let (sheet, secrets) = distinct_sheet(d);
    let fx = common::mini_fixture();
    for mode in ["quant", "align"] {
        let out = d.join(format!("run_{mode}"));
        let o = qrnamap_no_tools(
            &[
                "run", "--samples", sheet.to_str().unwrap(), "--reference", fx.join("genome.fa").to_str().unwrap(),
                "--gff", fx.join("genes.gff3").to_str().unwrap(), "--mode", mode, "--out", out.to_str().unwrap(), "--slurm",
            ],
            &empty,
        );
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        let mut n = 0;
        for e in fs::read_dir(out.join("slurm")).unwrap() {
            let p = e.unwrap().path();
            if p.is_file() {
                n += 1;
                let text = fs::read_to_string(&p).unwrap();
                for sec in &secrets {
                    assert!(!text.contains(sec.as_str()), "{} contains {sec}", p.display());
                }
                for part in ["ZZsample", "QQfastq"] {
                    assert!(!text.contains(part), "{} contains {part}", p.display());
                }
            }
        }
        assert_eq!(n, 4);
        let stdout = String::from_utf8_lossy(&o.stdout);
        assert!(!stdout.contains("ZZsample") && !stdout.contains("QQfastq"), "{stdout}");
    }
}

#[test]
fn plain_run_dir_gets_sbatch_output_directive() {
    let s = slurm_scripts(&plan(7, 3, Path::new("/bigdata/lab/run1"))).unwrap();
    assert!(s.samples.contains("#SBATCH --output=/bigdata/lab/run1/slurm/logs/%x_%A_%a.out\n"), "{}", s.samples);
    assert!(s.prepare.contains("#SBATCH --output=/bigdata/lab/run1/slurm/logs/%x_%j.out\n"), "{}", s.prepare);
}

// ---------------------------------------------------------------------------------------------
// Ruling R28: tool paths baked into the scripts, tool check, kill-on-invalid-dep, recovery.

/// A sheet of `n` paired samples with tiny real FASTQ files, condition ctl/trt alternating.
fn small_sheet(d: &Path, n: usize) -> PathBuf {
    let mut t = String::from("sample\tfastq1\tfastq2\tcondition\n");
    for i in 0..n {
        let (r1, r2) = (d.join(format!("s{i}_1.fq")), d.join(format!("s{i}_2.fq")));
        fs::write(&r1, "@r\nACGT\n+\nIIII\n").unwrap();
        fs::write(&r2, "@r\nACGT\n+\nIIII\n").unwrap();
        let cond = if i % 2 == 0 { "ctl" } else { "trt" };
        t.push_str(&format!("s{i}\t{}\t{}\t{cond}\n", r1.display(), r2.display()));
    }
    let p = d.join("sheet.tsv");
    fs::write(&p, t).unwrap();
    p
}

/// `qrnamap run --slurm` with the given PATH and extra arguments.
fn generate(d: &Path, mode: &str, path_env: &str, extra: &[&str]) -> std::process::Output {
    let sheet = small_sheet(d, 4);
    let fx = common::mini_fixture();
    let out = d.join("run");
    Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args([
            "run", "--samples", sheet.to_str().unwrap(), "--reference", fx.join("genome.fa").to_str().unwrap(), "--gff",
            fx.join("genes.gff3").to_str().unwrap(), "--mode", mode, "--out", out.to_str().unwrap(), "--slurm",
        ])
        .args(extra)
        .env("PATH", path_env)
        .output()
        .unwrap()
}

fn read_script(d: &Path, name: &str) -> String {
    fs::read_to_string(d.join("run").join("slurm").join(name)).unwrap()
}

fn export_line(text: &str) -> String {
    text.lines().find(|l| l.starts_with("export PATH=")).unwrap_or_else(|| panic!("no export PATH line:\n{text}")).to_string()
}

#[test]
fn sbatch_files_export_baked_tool_paths_first() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let explicit = d.join("explicit_tools");
    fs::create_dir_all(&explicit).unwrap();
    let on_path = d.join("path_tools");
    stub_tools(&on_path, &["kallisto"], &d.join("log"));
    let o = generate(d, "quant", on_path.to_str().unwrap(), &["--tool-path", explicit.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let want = format!(
        "export PATH={}:{}:\"$PATH\"",
        sh_quote(fs::canonicalize(&explicit).unwrap().to_str().unwrap()).unwrap(),
        sh_quote(on_path.to_str().unwrap()).unwrap()
    );
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch"] {
        let text = read_script(d, f);
        assert_eq!(export_line(&text), want, "{f}");
        let pos = |needle: &str| text.find(needle).unwrap_or_else(|| panic!("{needle} in {f}"));
        assert!(pos("export PATH=") < pos("exec "), "{f}: PATH is set before the program starts");
    }
    // A directory given twice (explicitly and found on PATH) is written once.
    let o = generate(d, "quant", on_path.to_str().unwrap(), &["--tool-path", on_path.to_str().unwrap()]);
    assert!(o.status.success());
    let canon = sh_quote(fs::canonicalize(&on_path).unwrap().to_str().unwrap()).unwrap();
    assert_eq!(export_line(&read_script(d, "prepare.sbatch")), format!("export PATH={canon}:\"$PATH\""));
}

#[test]
fn generation_fails_when_required_tool_missing() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let empty = d.join("empty");
    fs::create_dir_all(&empty).unwrap();
    let o = generate(d, "quant", empty.to_str().unwrap(), &[]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(err.contains("required tools not found") && err.contains("kallisto"), "{err}");
    assert!(!d.join("run").exists(), "no scripts when a tool is missing");
    // Align mode with a contrast lists every missing tool.
    let o = generate(d, "align", empty.to_str().unwrap(), &["--contrast", "condition,trt,ctl"]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    for t in ["hisat2-build", "hisat2-build-s", "hisat2-align-s", "python", "featureCounts", "Rscript"] {
        assert!(err.contains(t), "{t}: {err}");
    }
    // A --tool-path directory with kallisto is enough.
    let tools = d.join("tools");
    stub_tools(&tools, &["kallisto"], &d.join("log"));
    let o = generate(d, "quant", empty.to_str().unwrap(), &["--tool-path", tools.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    // A --tool-path directory that does not exist is a named error.
    let o = generate(d, "quant", empty.to_str().unwrap(), &["--tool-path", d.join("nope").to_str().unwrap()]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success() && err.contains("--tool-path") && err.contains("nope"), "{err}");
}

/// Stub `qrnamap` binary: records PATH and starts `kallisto` by name.
fn stub_qrnamap(path: &Path, log: &Path) {
    let l = sh_quote(log.to_str().unwrap()).unwrap();
    write_exe(path, &format!("#!/bin/sh\nprintf 'PATH=%s\\n' \"$PATH\" >> {l}\nkallisto\n"));
}

#[test]
fn tool_path_hostile_directory_round_trips() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let tools = d.join("tools $HOME it's `echo pwned`");
    let log = d.join("calls.log");
    stub_tools(&tools, &["kallisto"], &log);
    let bin = d.join("qrnamap_stub");
    stub_qrnamap(&bin, &log);
    let empty = d.join("empty");
    fs::create_dir_all(&empty).unwrap();
    let o = generate(d, "quant", empty.to_str().unwrap(), &["--tool-path", tools.to_str().unwrap(), "--qrnamap-bin", bin.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch", "submit.sh"] {
        let o = Command::new("bash").arg("-n").arg(d.join("run/slurm").join(f)).output().unwrap();
        assert!(o.status.success(), "bash -n {f}: {}", String::from_utf8_lossy(&o.stderr));
    }
    let scratch = d.join("scratch");
    fs::create_dir_all(&scratch).unwrap();
    let o = Command::new("env")
        .args(["-i", "PATH=/usr/bin:/bin", "HOME=/home/EXPANDED"])
        .arg(format!("SCRATCH={}", scratch.display()))
        .arg("bash")
        .arg(d.join("run/slurm/prepare.sbatch"))
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let text = fs::read_to_string(&log).unwrap();
    let canon = fs::canonicalize(&tools).unwrap();
    assert_eq!(text, format!("PATH={}:/usr/bin:/bin\nkallisto\n", canon.display()), "the directory arrives literally");
}

#[test]
fn scripts_work_without_submitter_path() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let scratch = d.join("scratch");
    fs::create_dir_all(&scratch).unwrap();
    // Stub variant: the tool is found on the generating PATH only (no --tool-path).
    let tools = d.join("path_tools");
    let log = d.join("calls.log");
    stub_tools(&tools, &["kallisto"], &log);
    let bin = d.join("qrnamap_stub");
    stub_qrnamap(&bin, &log);
    let o = generate(d, "quant", tools.to_str().unwrap(), &["--qrnamap-bin", bin.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let run_prepare = || {
        Command::new("env")
            .args(["-i", "PATH=/usr/bin:/bin"])
            .arg(format!("SCRATCH={}", scratch.display()))
            .arg("bash")
            .arg(d.join("run/slurm/prepare.sbatch"))
            .output()
            .unwrap()
    };
    let o = run_prepare();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(fs::read_to_string(&log).unwrap().ends_with("kallisto\n"), "the stub tool was started through the baked PATH");

    // Real variant: real kallisto and the real qrnamap, prepare only, under `env -i`.
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let real = d.join("real");
    fs::create_dir_all(&real).unwrap();
    let o = generate(&real, "quant", &std::env::var("PATH").unwrap(), &["--qrnamap-bin", env!("CARGO_BIN_EXE_qrnamap")]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let o = Command::new("env")
        .args(["-i", "PATH=/usr/bin:/bin"])
        .arg(format!("SCRATCH={}", scratch.display()))
        .arg("bash")
        .arg(real.join("run/slurm/prepare.sbatch"))
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(real.join("run/.prepared").is_file() && real.join("run/transcripts.kidx").is_file());
}

#[test]
fn prepare_tool_check_lists_missing_tools() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let tools = d.join("tools");
    let log = d.join("calls.log");
    stub_tools(&tools, ALL_TOOLS, &log);
    let bin = d.join("qrnamap_stub");
    stub_qrnamap(&bin, &log);
    let empty = d.join("empty");
    fs::create_dir_all(&empty).unwrap();
    let o = generate(
        d,
        "align",
        empty.to_str().unwrap(),
        &["--tool-path", tools.to_str().unwrap(), "--qrnamap-bin", bin.to_str().unwrap(), "--contrast", "condition,trt,ctl"],
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    // The tools disappear after generation (for example a removed environment).
    for t in ["featureCounts", "Rscript"] {
        fs::remove_file(tools.join(t)).unwrap();
    }
    let scratch = d.join("scratch");
    fs::create_dir_all(&scratch).unwrap();
    let o = Command::new("env")
        .args(["-i", "PATH=/nonexistent"])
        .arg(format!("SCRATCH={}", scratch.display()))
        .arg("/bin/bash")
        .arg(d.join("run/slurm/prepare.sbatch"))
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(err.contains("required tools not found on PATH: featureCounts Rscript"), "{err}");
    assert!(!log.exists(), "the program must not start when a tool is missing");
}

fn run_submit(sub: &Path, d: &Path) -> Vec<Vec<String>> {
    let stub = d.join("stub_sbatch");
    fs::create_dir_all(&stub).unwrap();
    let log = d.join("sbatch.log");
    write_exe(
        &stub.join("sbatch"),
        "#!/bin/sh\nn=$(grep -c '^CALL$' \"$STUB_LOG\" 2>/dev/null || true)\nn=$((${n:-0}+1))\n{ echo CALL; printf '%s\\n' \"$@\"; } >> \"$STUB_LOG\"\necho \"100$n\"\n",
    );
    let path = format!("{}:{}", stub.display(), std::env::var("PATH").unwrap());
    let o = Command::new("bash").arg(sub).env("PATH", path).env("STUB_LOG", &log).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    calls(&log)
}

#[test]
fn submit_sh_uses_kill_on_invalid_dep_for_dependents() {
    let tmp = tempfile::tempdir().unwrap();
    let p = plan(7, 3, &tmp.path().join("run"));
    let sub = write_slurm_scripts(&p).unwrap();
    let c = run_submit(&sub, tmp.path());
    assert_eq!(c.len(), 3);
    let flag = "--kill-on-invalid-dep=yes".to_string();
    assert!(!c[0].contains(&flag), "prepare has no dependency: {:?}", c[0]);
    assert!(c[1].contains(&flag) && c[1].contains(&"--dependency=afterok:1001".to_string()), "{:?}", c[1]);
    assert!(c[2].contains(&flag) && c[2].contains(&"--dependency=afterok:1002".to_string()), "{:?}", c[2]);
}

#[test]
fn submit_sh_comments_explain_recovery() {
    let run_dir = Path::new("/bigdata/lab/run 1");
    let s = slurm_scripts(&plan(7, 3, run_dir)).unwrap();
    let q = |name: &str| sh_quote(&format!("/bigdata/lab/run 1/slurm/{name}")).unwrap();
    let logs = |pat: &str| sh_quote(&format!("/bigdata/lab/run 1/slurm/logs/{pat}")).unwrap();
    let comments: String = s.submit.lines().filter(|l| l.starts_with('#')).map(|l| format!("{l}\n")).collect();
    for want in [
        "scancel FINAL_JOB_ID".to_string(),
        format!("sbatch --array=TASK_IDS --output={} {}", logs("%x_%A_%a.out"), q("samples.sbatch")),
        format!("sbatch --kill-on-invalid-dep=yes --dependency=afterok:NEW_ARRAY_JOB_ID --output={} {}", logs("%x_%j.out"), q("final.sbatch")),
    ] {
        assert!(comments.contains(&want), "missing {want:?} in\n{comments}");
    }
    // The recovery text is only comments: the script runs three sbatch calls.
    assert_eq!(s.submit.lines().filter(|l| !l.starts_with('#') && l.contains("sbatch")).count(), 3, "{}", s.submit);
}

// ---------------------------------------------------------------------------------------------
// Task 5: 8 GB default, --tmpdir for `run --slurm`, SCRATCH message

const SCRATCH_MESSAGE: &str = "SCRATCH is not set: this job needs node-local scratch. Submit from a login shell where $SCRATCH is set, or regenerate the scripts with --tmpdir DIR.";

#[test]
fn default_mem_is_8_gb_in_scripts() {
    assert_eq!(DEFAULT_MEM_GB, 8);
    let s = slurm_scripts(&plan(7, 3, Path::new("/bigdata/lab/run1"))).unwrap();
    for (name, text) in [("prepare", &s.prepare), ("samples", &s.samples), ("final", &s.final_)] {
        assert!(text.contains("#SBATCH --mem=8G\n"), "{name}:\n{text}");
    }
}

#[test]
fn mem_gb_zero_is_still_an_error() {
    let mut p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    p.mem_gb = 0;
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::BadValue { what: "--mem-gb", .. })));
}

#[test]
fn tmpdir_option_replaces_scratch_check() {
    let mut p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    p.tmpdir = Some(PathBuf::from("/bigdata/lab/tmp dir"));
    for (name, text) in all_scripts(&p) {
        assert!(!text.contains("SCRATCH"), "{name}:\n{text}");
        if name.ends_with(".sbatch") {
            assert!(text.contains("--tmpdir '/bigdata/lab/tmp dir'\n"), "{name}:\n{text}");
            assert!(text.contains("set -euo pipefail"), "{name}");
        }
    }
}

#[test]
fn tmpdir_with_newline_is_refused() {
    let mut p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    p.tmpdir = Some(PathBuf::from("/bigdata/lab/t\nmp"));
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::NewlineOrNul(_))));
    let tmp = tempfile::tempdir().unwrap();
    let run_dir = tmp.path().join("run");
    let mut q = plan(7, 3, &run_dir);
    q.tmpdir = Some(PathBuf::from("/bigdata/lab/t\nmp"));
    assert!(matches!(write_slurm_scripts(&q), Err(SlurmError::NewlineOrNul(_))));
    assert!(!run_dir.exists(), "no run directory is written when --tmpdir is refused");
    p.tmpdir = Some(PathBuf::from("relative/tmp"));
    assert!(matches!(slurm_scripts(&p), Err(SlurmError::RelativePath(_))));
}

/// Stub `qrnamap`: records every argument on its own line after a `CALL` marker, then starts
/// `kallisto` so that the prepare tool check has something to find.
fn stub_qrnamap_args(path: &Path, log: &Path) {
    let l = sh_quote(log.to_str().unwrap()).unwrap();
    write_exe(path, &format!("#!/bin/sh\n{{ echo CALL; printf '%s\\n' \"$@\"; }} >> {l}\n"));
}

#[test]
fn without_tmpdir_scratch_check_remains_with_new_message() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let tools = d.join("tools");
    let log = d.join("calls.log");
    stub_tools(&tools, &["kallisto"], &log);
    let bin = d.join("qrnamap_stub");
    stub_qrnamap_args(&bin, &log);
    let o = generate(d, "quant", tools.to_str().unwrap(), &["--qrnamap-bin", bin.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch"] {
        let o = Command::new("env")
            .args(["-i", &format!("PATH={}:/usr/bin:/bin", tools.display())])
            .arg("SLURM_ARRAY_TASK_ID=0")
            .arg("bash")
            .arg(d.join("run/slurm").join(f))
            .output()
            .unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success(), "{f}");
        let line = err.lines().next().unwrap_or("");
        assert!(line.ends_with(&format!("SCRATCH: {SCRATCH_MESSAGE}")), "{f}: {err:?}");
        assert_eq!(err.lines().count(), 1, "{f}: {err:?}");
        assert!(!log.exists(), "the program must not start without scratch");
    }
}

#[test]
fn hostile_tmpdir_round_trips() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let hostile = d.join("tmp dir $HOME it's `echo pwned` 100% back\\slash; a&b *");
    let tools = d.join("tools");
    let log = d.join("calls.log");
    stub_tools(&tools, &["kallisto"], &d.join("tool.log"));
    let bin = d.join("qrnamap_stub");
    stub_qrnamap_args(&bin, &log);
    let o = generate(d, "quant", tools.to_str().unwrap(), &["--qrnamap-bin", bin.to_str().unwrap(), "--tmpdir", hostile.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch", "submit.sh"] {
        let o = Command::new("bash").arg("-n").arg(d.join("run/slurm").join(f)).output().unwrap();
        assert!(o.status.success(), "bash -n {f}: {}", String::from_utf8_lossy(&o.stderr));
    }
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch"] {
        assert!(!read_script(d, f).contains("SCRATCH"), "{f}");
        // No SCRATCH in the environment: the job must still run.
        let o = Command::new("env")
            .args(["-i", "PATH=/usr/bin:/bin", "HOME=/home/EXPANDED", "SLURM_ARRAY_TASK_ID=0"])
            .arg("bash")
            .arg(d.join("run/slurm").join(f))
            .output()
            .unwrap();
        assert!(o.status.success(), "{f}: {}", String::from_utf8_lossy(&o.stderr));
    }
    let want = hostile.to_str().unwrap();
    let all = calls(&log);
    assert_eq!(all.len(), 3);
    for c in &all {
        let i = c.iter().position(|a| a == "--tmpdir").expect("--tmpdir present");
        assert_eq!(c[i + 1], want, "one literal argument");
        assert_eq!(c.iter().filter(|a| *a == want).count(), 1);
    }
}

#[test]
fn scripts_never_contain_sample_names_or_fastq_paths_with_tmpdir() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let empty = d.join("empty_path");
    fs::create_dir_all(&empty).unwrap();
    let (sheet, secrets) = distinct_sheet(d);
    let fx = common::mini_fixture();
    let out = d.join("run_tmpdir");
    let o = qrnamap_no_tools(
        &[
            "run", "--samples", sheet.to_str().unwrap(), "--reference", fx.join("genome.fa").to_str().unwrap(),
            "--gff", fx.join("genes.gff3").to_str().unwrap(), "--mode", "quant", "--out", out.to_str().unwrap(), "--slurm",
            "--tmpdir", d.join("job tmp").to_str().unwrap(),
        ],
        &empty,
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for e in fs::read_dir(out.join("slurm")).unwrap() {
        let p = e.unwrap().path();
        if p.is_file() {
            let text = fs::read_to_string(&p).unwrap();
            for sec in &secrets {
                assert!(!text.contains(sec.as_str()), "{} contains {sec}", p.display());
            }
            assert!(!text.contains("ZZsample") && !text.contains("QQfastq"), "{}", p.display());
        }
    }
}

#[test]
fn tmpdir_under_scratch_warns() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let tools = d.join("tools");
    stub_tools(&tools, &["kallisto"], &d.join("log"));
    let o = generate(d, "quant", tools.to_str().unwrap(), &["--tmpdir", "/scratch/someone/12345/jobtmp"]);
    assert!(o.status.success(), "a warning, not an error: {}", String::from_utf8_lossy(&o.stderr));
    let err = String::from_utf8_lossy(&o.stderr);
    let w = err.lines().find(|l| l.contains("--tmpdir")).unwrap_or_else(|| panic!("no --tmpdir warning: {err}"));
    assert!(w.starts_with("warning:") && w.contains("job-specific") && w.contains("other nodes"), "{w}");
    // A directory under $SCRATCH warns too.
    let d2 = d.join("two");
    fs::create_dir_all(&d2).unwrap();
    let sheet = small_sheet(&d2, 2);
    let fx = common::mini_fixture();
    let o = Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args([
            "run", "--samples", sheet.to_str().unwrap(), "--reference", fx.join("genome.fa").to_str().unwrap(), "--gff",
            fx.join("genes.gff3").to_str().unwrap(), "--mode", "quant", "--out", d2.join("run").to_str().unwrap(), "--slurm",
            "--tmpdir", "/node/local/x/jobtmp",
        ])
        .env("PATH", tools.to_str().unwrap())
        .env("SCRATCH", "/node/local/x")
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).lines().any(|l| l.contains("--tmpdir")), "{}", String::from_utf8_lossy(&o.stderr));
}

#[test]
fn tmpdir_under_scratchy_does_not_warn() {
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let tools = d.join("tools");
    stub_tools(&tools, &["kallisto"], &d.join("log"));
    let o = generate(d, "quant", tools.to_str().unwrap(), &["--tmpdir", "/scratchy/jobtmp"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!err.contains("--tmpdir"), "{err}");
}

// ---------------------------------------------------------------------------------------------
// Final review fix wave (R35-6): timing pointer and the conditional --mem comment

const MEM_DEFAULT_LINE: &str = "# The --mem default is small:";

#[test]
fn header_points_to_amb_timing_not_to_missing_timing() {
    let s = slurm_scripts(&plan(7, 3, Path::new("/bigdata/lab/run1"))).unwrap();
    for (name, text) in [("prepare", &s.prepare), ("samples", &s.samples), ("final", &s.final_)] {
        assert!(!text.contains("timing exists yet"), "{name}:\n{text}");
        assert!(
            text.contains("quant 22.7 to 26.8 s and align 35.2 to 52.1 s per million read pairs"),
            "{name}:\n{text}"
        );
        assert!(text.contains("validation/amb/README.md"), "{name}:\n{text}");
        assert!(text.contains("#SBATCH --time=02:00:00\n"), "{name}:\n{text}");
        assert!(text.contains("The --time default is a placeholder"), "{name}:\n{text}");
    }
}

#[test]
fn mem_default_comment_only_with_the_default() {
    let s = slurm_scripts(&plan(7, 3, Path::new("/bigdata/lab/run1"))).unwrap();
    for (name, text) in [("prepare", &s.prepare), ("samples", &s.samples), ("final", &s.final_)] {
        assert!(text.contains(MEM_DEFAULT_LINE), "{name}:\n{text}");
    }
    let mut p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    p.mem_gb = 16;
    let s = slurm_scripts(&p).unwrap();
    for (name, text) in [("prepare", &s.prepare), ("samples", &s.samples), ("final", &s.final_)] {
        assert!(!text.contains(MEM_DEFAULT_LINE), "{name}:\n{text}");
        assert!(text.contains("#SBATCH --mem=16G\n"), "{name}:\n{text}");
    }
}

// ---------------------------------------------------------------------------------------------
// Salmon: --quantifier in the scripts, tools and the rate note

#[test]
fn slurm_forwards_quantifier() {
    use qrnamap_cli::run::required_tools;
    let tmp = tempfile::tempdir().unwrap();
    let d = tmp.path();
    let tools = d.join("tools");
    stub_tools(&tools, &["salmon"], &d.join("log"));
    let o = generate(d, "quant", tools.to_str().unwrap(), &["--quantifier", "salmon"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch"] {
        let text = read_script(d, f);
        assert!(text.contains("'--quantifier' 'salmon'"), "{f}:\n{text}");
    }
    let prep = read_script(d, "prepare.sbatch");
    assert!(prep.contains("for tool in 'salmon'; do"), "{prep}");
    // Salmon needs salmon on PATH, not kallisto.
    let o = generate(d, "quant", tools.to_str().unwrap(), &[]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success() && err.contains("kallisto"), "{err}");
    // The default is written out too, so the scripts are explicit.
    stub_tools(&tools, &["kallisto"], &d.join("log"));
    let o = generate(d, "quant", tools.to_str().unwrap(), &[]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch"] {
        assert!(read_script(d, f).contains("'--quantifier' 'kallisto'"), "{f}");
    }
    // Align scripts carry no --quantifier.
    let align_tools =
        ["hisat2-build", "hisat2-build-s", "hisat2-align-s", "python", "featureCounts"];
    stub_tools(&tools, &align_tools, &d.join("log"));
    let o = generate(d, "align", tools.to_str().unwrap(), &[]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch"] {
        assert!(!read_script(d, f).contains("--quantifier"), "{f}");
    }
    assert_eq!(required_tools(Mode::Quant, Quantifier::Salmon, false), vec!["salmon"]);
    assert_eq!(required_tools(Mode::Quant, Quantifier::Kallisto, false), vec!["kallisto"]);
    assert_eq!(required_tools(Mode::Quant, Quantifier::Salmon, true), vec!["salmon", "Rscript"]);
    let align = required_tools(Mode::Align, Quantifier::Salmon, false);
    assert!(!align.contains(&"salmon") && !align.contains(&"kallisto"), "{align:?}");
}

#[test]
fn slurm_rate_note_names_kallisto_only() {
    let mut p = plan(7, 3, Path::new("/bigdata/lab/run1"));
    let s = slurm_scripts(&p).unwrap();
    for (name, text) in [("prepare", &s.prepare), ("samples", &s.samples), ("final", &s.final_)] {
        assert!(text.contains("The quant rates were measured with kallisto."), "{name}:\n{text}");
        assert!(!text.contains("salmon"), "{name}:\n{text}");
    }
    p.quantifier = Some(Quantifier::Salmon);
    let s = slurm_scripts(&p).unwrap();
    for (name, text) in [("prepare", &s.prepare), ("samples", &s.samples), ("final", &s.final_)] {
        assert!(text.contains("The quant rates were measured with kallisto."), "{name}:\n{text}");
        assert!(text.contains("Salmon rates are not measured."), "{name}:\n{text}");
    }
}
