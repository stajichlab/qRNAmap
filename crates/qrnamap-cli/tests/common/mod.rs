//! Helpers shared by integration tests that call external tools.
//! Copied from `crates/qrnamap-core/tests/common/mod.rs` (no shared dev crate).
#![allow(dead_code)]

pub mod fetch_server;
pub mod large;
pub mod sim;

use std::path::{Path, PathBuf};

use qrnamap_core::{build_transcripts, Sample, SampleSheet, TxRecord};
use qseq_io::{parse_gff3, read_fasta};

/// Pairs per transcript, in fixture order tA1, tB1, tB2, tR1, tC1, tU1. Sum 10,000.
pub const COUNTS: [u64; 6] = [3000, 1500, 1000, 2000, 1500, 1000];

/// The six mini fixture transcripts, built by the core builder.
pub fn records() -> Vec<TxRecord> {
    let open = |p: PathBuf| std::io::BufReader::new(std::fs::File::open(p).unwrap());
    let genome = read_fasta(open(mini_genome())).unwrap();
    let model = parse_gff3(open(mini_gff())).unwrap();
    let recs = build_transcripts(&genome, &model).unwrap();
    let ids: Vec<&str> = recs.iter().map(|r| r.id.as_str()).collect();
    assert_eq!(ids, ["tA1", "tB1", "tB2", "tR1", "tC1", "tU1"]);
    recs
}

/// A paired sample parsed through the real sample sheet parser.
pub fn sample(name: &str, r1: &[&Path], r2: &[&Path]) -> Sample {
    let j = |v: &[&Path]| {
        v.iter()
            .map(|p| p.to_str().unwrap())
            .collect::<Vec<_>>()
            .join(",")
    };
    let sheet = format!("sample\tfastq1\tfastq2\n{name}\t{}\t{}\n", j(r1), j(r2));
    SampleSheet::parse(&sheet).unwrap().samples.remove(0)
}

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
