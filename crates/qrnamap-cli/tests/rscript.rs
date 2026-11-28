//! Finding an Rscript that loads DESeq2 (`deseq2_version`, `run_de` with `DeArgs::rscript`), with stand-in Rscript
//! scripts. Unix only (the stand-ins are shell scripts).
#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

use qrnamap_cli::de::{deseq2_version, run_de, DeArgs};

fn script(dir: &Path, name: &str, body: &str) -> PathBuf {
    let p = dir.join(name);
    fs::write(&p, format!("#!/bin/sh\n{body}\n")).unwrap();
    fs::set_permissions(&p, fs::Permissions::from_mode(0o755)).unwrap();
    p
}

#[test]
fn deseq2_version_reads_the_version_or_names_the_failure() {
    let d = tempfile::tempdir().unwrap();
    let ok = script(d.path(), "ok", "printf '1.50.2'");
    assert_eq!(deseq2_version(&ok).unwrap(), "1.50.2");

    let missing = script(
        d.path(),
        "missing",
        "echo \"Error in packageVersion('DESeq2') : there is no package called 'DESeq2'\" >&2; echo 'Execution halted' >&2; exit 1",
    );
    let e = deseq2_version(&missing).unwrap_err();
    assert!(e.contains("cannot load DESeq2") && e.contains("Execution halted"), "{e}");

    let e = deseq2_version(&d.path().join("no-such-rscript")).unwrap_err();
    assert!(e.starts_with("cannot start"), "{e}");
}

/// `run_de` runs the Rscript named in `DeArgs::rscript`: a stand-in that records its arguments and fails shows the inputs
/// were checked and handed over (the script, counts, coldata and prefix arguments).
#[test]
fn run_de_uses_the_given_rscript() {
    let d = tempfile::tempdir().unwrap();
    let log = d.path().join("args.log");
    let fake = script(d.path(), "Rscript", &format!("echo \"$@\" > {}; echo 'stand-in R' >&2; exit 3", log.display()));
    fs::write(d.path().join("counts.tsv"), "gene\ta1\ta2\tb1\tb2\ng1\t10\t12\t30\t33\ng2\t5\t0\t7\t9\n").unwrap();
    fs::write(
        d.path().join("samples.tsv"),
        "sample\tfastq1\tfastq2\tcondition\na1\ta1.fq\t\tctl\na2\ta2.fq\t\tctl\nb1\tb1.fq\t\ttrt\nb2\tb2.fq\t\ttrt\n",
    )
    .unwrap();
    let args = DeArgs {
        counts: d.path().join("counts.tsv"),
        samples: d.path().join("samples.tsv"),
        design: "~ condition".to_string(),
        contrast: "condition,trt,ctl".to_string(),
        shrink: None,
        out: d.path().join("out").join("de"),
        tmpdir: Some(d.path().join("tmp")),
        rscript: Some(fake),
    };
    let e = run_de(&args).unwrap_err().to_string();
    assert!(e.contains("stand-in R"), "{e}");
    let called = fs::read_to_string(&log).unwrap();
    assert!(called.starts_with("--vanilla ") && called.contains("de.R") && called.contains("counts.tsv") && called.contains("coldata.tsv"), "{called}");
}
