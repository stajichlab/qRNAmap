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
    assert!(text.contains("unrecognized option '--bogus-flag'"), "{text}");
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
