//! SLURM script generation. No real SLURM job is submitted: `submit.sh` runs against a stub
//! `sbatch` on a temporary PATH that records its arguments.

mod common;

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

use qrnamap_cli::slurm::{sh_quote, slurm_array_script, slurm_scripts, write_slurm_scripts, Plan, SlurmError};
use qrnamap_core::Mode;

fn plan(n: usize, k: usize, run_dir: &Path) -> Plan {
    Plan {
        n_samples: n,
        per_task: k,
        cpus: 4,
        mem_gb: 16,
        time: "02:00:00".into(),
        partition: Some("batch".into()),
        account: None,
        run_dir: run_dir.to_path_buf(),
        work_dir: PathBuf::from("/bigdata/lab/work"),
        qrnamap_bin: PathBuf::from("/opt/qrnamap/bin/qrnamap"),
        mode: Mode::Quant,
        run_args: vec!["--samples".into(), "/bigdata/lab/sheet.tsv".into(), "--out".into(), run_dir.to_str().unwrap().into()],
        task_args: vec!["--samples".into(), "/bigdata/lab/sheet.tsv".into(), "--out".into(), run_dir.to_str().unwrap().into()],
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
            for want in ["#SBATCH --cpus-per-task=4", "#SBATCH --mem=16G", "#SBATCH --time=02:00:00", "#SBATCH --partition=batch", "#SBATCH --job-name="] {
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
        vec!["--parsable".to_string(), "--dependency=afterok:1001".into(), format!("--output={logs}/%x_%A_%a.out"), format!("{rd}/slurm/samples.sbatch")]
    );
    assert_eq!(
        c[2],
        vec!["--parsable".to_string(), "--dependency=afterok:1002".into(), format!("--output={logs}/%x_%j.out"), format!("{rd}/slurm/final.sbatch")]
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

fn qrnamap_no_tools(args: &[&str], empty_path: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qrnamap")).args(args).env("PATH", empty_path).output().unwrap()
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
