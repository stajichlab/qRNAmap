//! Helpers shared by integration tests that call external tools.
//! Copied from `crates/qrnamap-core/tests/common/mod.rs` (no shared dev crate).
#![allow(dead_code)]

pub mod sim;

use std::path::PathBuf;

/// Find `name` on PATH. When it is missing: print a skip message and return `None`,
/// or panic when `QRNAMAP_REQUIRE_TOOLS=1` (ruling R5).
pub fn require_tool(name: &str) -> Option<PathBuf> {
    let found = std::env::var_os("PATH").and_then(|paths| {
        std::env::split_paths(&paths)
            .map(|d| d.join(name))
            .find(|p| p.is_file())
    });
    if found.is_none() {
        let required = std::env::var("QRNAMAP_REQUIRE_TOOLS")
            .map(|v| v == "1")
            .unwrap_or(false);
        if required {
            panic!("required tool '{name}' is not on PATH and QRNAMAP_REQUIRE_TOOLS=1");
        }
        eprintln!("SKIP: tool '{name}' is not on PATH (set QRNAMAP_REQUIRE_TOOLS=1 to make this a failure)");
    }
    found
}

/// True when every tool is on PATH (each missing tool prints a skip line or panics).
pub fn require_tools(names: &[&str]) -> bool {
    names.iter().all(|n| require_tool(n).is_some())
}

/// Path of the shared mini fixture directory at the workspace root.
pub fn mini_fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/mini")
}

pub fn mini_genome() -> PathBuf {
    mini_fixture().join("genome.fa")
}

pub fn mini_gff() -> PathBuf {
    mini_fixture().join("genes.gff3")
}
