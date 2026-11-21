//! `run_call` tests with the real tools (ruling R5: skip when missing, fail under
//! `QRNAMAP_REQUIRE_TOOLS=1`).

mod common;

use qrnamap_cli::exec::{run_call, ExecError};
use qrnamap_core::{Tool, ToolCall};

#[test]
fn native_failure_reports_stderr() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let call = ToolCall {
        tool: Tool::KallistoQuant,
        args: vec!["quant".into(), "--bogus-flag".into()],
        inputs: vec![],
        outputs: vec![],
    };
    let err = run_call(&call, dir.path()).unwrap_err();
    assert!(matches!(err, ExecError::Failed { .. }), "{err:?}");
    let text = err.to_string();
    assert!(text.contains("kallisto"), "{text}");
    // kallisto's own words, from its stderr.
    assert!(
        text.contains("unrecognized option '--bogus-flag'"),
        "{text}"
    );
}

#[test]
fn run_call_missing_output_is_error() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("never_written.idx");
    let call = ToolCall {
        tool: Tool::KallistoIndex,
        args: vec!["version".into()],
        inputs: vec![],
        outputs: vec![missing.clone()],
    };
    let err = run_call(&call, dir.path()).unwrap_err();
    match &err {
        ExecError::MissingOutput { path, .. } => assert_eq!(path, &missing),
        other => panic!("expected MissingOutput, got {other:?}"),
    }
    assert!(err.to_string().contains("never_written.idx"));
}

#[test]
fn run_call_missing_input_is_error_before_running() {
    let dir = tempfile::tempdir().unwrap();
    let call = ToolCall {
        tool: Tool::Hisat2Align,
        args: vec![],
        inputs: vec![dir.path().join("no_index")],
        outputs: vec![],
    };
    let err = run_call(&call, dir.path()).unwrap_err();
    assert!(matches!(err, ExecError::MissingInput { .. }), "{err:?}");
}

#[test]
fn stderr_tail_is_bounded() {
    let long = "x".repeat(100_000);
    let t = qrnamap_cli::exec::tail(long.as_bytes());
    assert!(t.len() <= qrnamap_cli::exec::TAIL_BYTES + 3, "{}", t.len());
}

#[test]
fn error_tail_strips_ansi_colour_codes() {
    use qrnamap_cli::exec::error_tail;
    let coloured = b"\x1b[1;31merror\x1b[0m: bad index \x1b[2K\x1b[33mwarn\x1b[m\n";
    assert_eq!(error_tail(coloured), "error: bad index warn\n");
    // Plain text, a lone ESC and a cut-off sequence at the end.
    assert_eq!(error_tail(b"plain [1m text\n"), "plain [1m text\n");
    assert_eq!(error_tail(b"a\x1bb"), "a\x1bb");
    assert_eq!(error_tail(b"end\x1b[31"), "end");
}

#[test]
fn failed_tool_error_has_no_colour_codes() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    // A stub kallisto that writes coloured stderr and fails.
    let bin = d.join("bin");
    std::fs::create_dir_all(&bin).unwrap();
    let stub = bin.join("kallisto");
    let body = "#!/bin/sh\nprintf '\\033[1;31mERROR\\033[0m: stub failure\\n' >&2\nexit 1\n";
    std::fs::write(&stub, body).unwrap();
    std::fs::set_permissions(&stub, std::fs::Permissions::from_mode(0o755)).unwrap();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("tx2gene.tsv.zst"));
    let (g, gff) = (common::mini_genome(), common::mini_gff());
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    let (r1, r2) = (d.join("r_1.fq"), d.join("r_2.fq"));
    std::fs::write(&r1, "@r\nACGT\n+\nIIII\n").unwrap();
    std::fs::write(&r2, "@r\nACGT\n+\nIIII\n").unwrap();
    let sheet = d.join("sheet.tsv");
    let text = format!(
        "sample\tfastq1\tfastq2\nS1\t{}\t{}\n",
        r1.display(),
        r2.display()
    );
    std::fs::write(&sheet, text).unwrap();
    let s = |p: &std::path::Path| p.to_str().unwrap().to_string();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_qrnamap"))
        .args([
            "quant",
            "--samples",
            &s(&sheet),
            "--transcripts",
            &s(&fa),
            "--tx2gene",
            &s(&t2g),
        ])
        .args(["--out", &s(&d.join("out")), "--tmpdir", &s(d)])
        .env("PATH", &bin)
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(err.contains("ERROR: stub failure"), "{err}");
    assert!(!err.contains('\x1b'), "{err:?}");
}
