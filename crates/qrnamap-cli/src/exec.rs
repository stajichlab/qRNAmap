//! Run a `ToolCall` as a native program.
//!
//! Safety contract (ruling R18): the program and its arguments go to
//! `std::process::Command` as separate argv elements. No shell is started, so a space,
//! quote or `$` in a path or sample name is never interpreted.

use std::fmt;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use qrnamap_core::{argv_native, Tool, ToolCall};

/// At most this many bytes of a tool's stderr (and stdout) are kept in an error.
pub const TAIL_BYTES: usize = 4096;

#[derive(Debug)]
pub enum ExecError {
    /// The program could not be started (for example, it is not on PATH).
    Spawn { program: String, source: std::io::Error },
    /// A declared input does not exist. Checked before the program runs.
    MissingInput { program: String, path: PathBuf },
    /// The program exited with a non-zero status or was killed by a signal.
    Failed { program: String, status: String, stderr_tail: String, stdout_tail: String },
    /// The program exited with status 0 but a declared output is missing.
    MissingOutput { program: String, path: PathBuf, stderr_tail: String },
}

impl fmt::Display for ExecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExecError::Spawn { program, source } => {
                write!(f, "cannot start '{program}': {source} (is it on PATH?)")
            }
            ExecError::MissingInput { program, path } => {
                write!(f, "'{program}': input {} does not exist", path.display())
            }
            ExecError::Failed { program, status, stderr_tail, stdout_tail } => {
                write!(f, "'{program}' failed ({status}). stderr (last {TAIL_BYTES} bytes at most):\n{stderr_tail}")?;
                if stderr_tail.trim().is_empty() && !stdout_tail.trim().is_empty() {
                    write!(f, "\nstdout (last {TAIL_BYTES} bytes at most):\n{stdout_tail}")?;
                }
                Ok(())
            }
            ExecError::MissingOutput { program, path, stderr_tail } => write!(
                f,
                "'{program}' exited with status 0 but did not write {}. stderr (last {TAIL_BYTES} bytes at most):\n{stderr_tail}",
                path.display()
            ),
        }
    }
}

impl std::error::Error for ExecError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ExecError::Spawn { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// The last `TAIL_BYTES` bytes of `bytes` as text (lossy UTF-8), with `...` in front when cut.
pub fn tail(bytes: &[u8]) -> String {
    if bytes.len() <= TAIL_BYTES {
        return String::from_utf8_lossy(bytes).into_owned();
    }
    let mut start = bytes.len() - TAIL_BYTES;
    // Do not start inside a UTF-8 sequence.
    while start < bytes.len() && (bytes[start] & 0b1100_0000) == 0b1000_0000 {
        start += 1;
    }
    format!("...{}", String::from_utf8_lossy(&bytes[start..]))
}

/// Remove ANSI CSI sequences (`ESC [`, parameter and intermediate bytes, one final byte in
/// `@`..`~`), such as colour codes. A sequence cut off at the end is removed. Any other byte,
/// a lone ESC too, is kept.
fn strip_csi(bytes: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == 0x1b && bytes.get(i + 1) == Some(&b'[') {
            let mut j = i + 2;
            while j < bytes.len() && (0x20..=0x3f).contains(&bytes[j]) {
                j += 1;
            }
            // Skip the final byte when it is there.
            i = if j < bytes.len() && (0x40..=0x7e).contains(&bytes[j]) { j + 1 } else { j };
            continue;
        }
        out.push(bytes[i]);
        i += 1;
    }
    out
}

/// The tool output tail kept in an error: `tail` after ANSI CSI sequences (colour codes) are
/// removed, so they do not reach the error message.
pub fn error_tail(bytes: &[u8]) -> String {
    tail(&strip_csi(bytes))
}

fn resolve(base: &Path, p: &Path) -> PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        base.join(p)
    }
}

/// An input exists. A HISAT2 index prefix (input of `Hisat2Align`) exists when its eight
/// `.N.ht2` files exist.
fn input_exists(tool: Tool, p: &Path) -> bool {
    if p.exists() {
        return true;
    }
    tool == Tool::Hisat2Align
        && (1..=8).all(|i| {
            let mut s = p.as_os_str().to_owned();
            s.push(format!(".{i}.ht2"));
            PathBuf::from(s).is_file()
        })
}

/// Run `argv_native(call)` with `scratch` as the working directory. No shell is used.
///
/// Before the run, every declared input must exist. After the run, the exit status must be 0
/// and every declared output must exist. Relative paths in `inputs`/`outputs` are checked
/// against `scratch`, as the program sees them; callers should pass absolute paths.
/// stdout and stderr are captured; errors carry a bounded tail of them.
pub fn run_call(call: &ToolCall, scratch: &Path) -> Result<(), ExecError> {
    run_call_stderr(call, scratch).map(|_| ())
}

/// `run_call`, returning the bounded stderr tail of a successful run (HISAT2 prints its
/// alignment summary there).
pub fn run_call_stderr(call: &ToolCall, scratch: &Path) -> Result<String, ExecError> {
    let argv = argv_native(call);
    let program = argv[0].clone();
    for p in &call.inputs {
        let p = resolve(scratch, p);
        if !input_exists(call.tool, &p) {
            return Err(ExecError::MissingInput { program, path: p });
        }
    }
    let out = Command::new(&program)
        .args(&argv[1..])
        .current_dir(scratch)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .map_err(|source| ExecError::Spawn { program: program.clone(), source })?;
    if !out.status.success() {
        return Err(ExecError::Failed {
            program,
            status: out.status.to_string(),
            stderr_tail: error_tail(&out.stderr),
            stdout_tail: tail(&out.stdout),
        });
    }
    for p in &call.outputs {
        let p = resolve(scratch, p);
        if !p.exists() {
            let stderr_tail = error_tail(&out.stderr);
            return Err(ExecError::MissingOutput { program, path: p, stderr_tail });
        }
    }
    Ok(tail(&out.stderr))
}
