//! Parity with gffread (`-w` transcript FASTA) and with HISAT2's `hisat2_extract_splice_sites.py`.
//! Each test skips with a message when a tool is missing, and fails instead when
//! `QRNAMAP_REQUIRE_TOOLS=1` (ruling R5). Every disagreement is printed.

mod common;

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{mini_fixture, require_tool};
use qrnamap_core::{build_transcripts, splice_sites};
use qseq_io::{parse_gff3, read_fasta, GeneModel};

fn work_dir(name: &str) -> PathBuf {
    let d =
        PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    // gffread writes genome.fa.fai next to the FASTA, so work on copies.
    for f in ["genome.fa", "genes.gff3"] {
        std::fs::copy(mini_fixture().join(f), d.join(f)).unwrap();
    }
    d
}

fn run(cmd: &mut Command) -> Vec<u8> {
    let out = cmd
        .output()
        .unwrap_or_else(|e| panic!("cannot run {cmd:?}: {e}"));
    assert!(
        out.status.success(),
        "{cmd:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn model(dir: &Path) -> GeneModel {
    parse_gff3(std::fs::read(dir.join("genes.gff3")).unwrap().as_slice()).unwrap()
}

#[test]
fn gffread_w_sequences_match() {
    let Some(gffread) = require_tool("gffread") else {
        return;
    };
    let d = work_dir("gffread_w");
    run(Command::new(&gffread).current_dir(&d).args([
        "-g",
        "genome.fa",
        "-w",
        "gffread.fa",
        "genes.gff3",
    ]));
    let theirs = read_fasta(std::fs::read(d.join("gffread.fa")).unwrap().as_slice()).unwrap();
    let genome = read_fasta(std::fs::read(d.join("genome.fa")).unwrap().as_slice()).unwrap();
    let ours = build_transcripts(&genome, &model(&d)).unwrap();

    let mut problems = Vec::new();
    for (id, seq) in &theirs {
        match ours.iter().find(|r| &r.id == id) {
            None => problems.push(format!("{id}: in gffread output only")),
            Some(r) if r.seq != seq.to_ascii_uppercase() => problems.push(format!(
                "{id}: sequence differs\n  gffread: {}\n  ours:    {}",
                String::from_utf8_lossy(&seq.to_ascii_uppercase()),
                String::from_utf8_lossy(&r.seq)
            )),
            Some(_) => {}
        }
    }
    for r in &ours {
        if !theirs.iter().any(|(id, _)| id == &r.id) {
            problems.push(format!("{}: in our output only", r.id));
        }
    }
    for p in &problems {
        eprintln!("gffread -w disagreement: {p}");
    }
    eprintln!(
        "gffread -w: {} records compared, {} disagreements",
        theirs.len(),
        problems.len()
    );
    assert!(
        problems.is_empty(),
        "{} disagreements with gffread -w (see stderr)",
        problems.len()
    );
}

fn hisat2_sites(gffread: &Path, script: &Path, d: &Path, extra: &[&str]) -> BTreeSet<String> {
    let mut gtf_cmd = Command::new(gffread);
    gtf_cmd
        .current_dir(d)
        .arg("-T")
        .args(extra)
        .args(["-o", "genes.gtf", "genes.gff3"]);
    run(&mut gtf_cmd);
    let out = run(Command::new(script).current_dir(d).arg("genes.gtf"));
    String::from_utf8(out)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

fn report(
    label: &str,
    ours: &BTreeSet<String>,
    theirs: &BTreeSet<String>,
) -> (Vec<String>, Vec<String>) {
    let only_ours: Vec<String> = ours.difference(theirs).cloned().collect();
    let only_theirs: Vec<String> = theirs.difference(ours).cloned().collect();
    eprintln!(
        "{label}: ours {} lines, hisat2 {} lines",
        ours.len(),
        theirs.len()
    );
    for l in &only_ours {
        eprintln!("{label}: only in ours:   {l}");
    }
    for l in &only_theirs {
        eprintln!("{label}: only in hisat2: {l}");
    }
    (only_ours, only_theirs)
}

#[test]
fn hisat2_splice_sites_match() {
    let Some(gffread) = require_tool("gffread") else {
        return;
    };
    let Some(script) = require_tool("hisat2_extract_splice_sites.py") else {
        return;
    };
    let d = work_dir("hisat2_ss");
    let ours: BTreeSet<String> = splice_sites(&model(&d))
        .lines()
        .map(str::to_string)
        .collect();

    // `gffread -T --force-exons` writes exon lines for the CDS-only gene: the sets must be equal.
    let forced = hisat2_sites(&gffread, &script, &d, &["--force-exons"]);
    let (a, b) = report("gffread -T --force-exons", &ours, &forced);
    assert!(
        a.is_empty() && b.is_empty(),
        "splice sites differ from hisat2 on the --force-exons GTF (see stderr)"
    );

    // Plain `gffread -T` writes only CDS lines for the CDS-only gene geneC, and the HISAT2 script reads
    // exon lines only, so it misses the geneC intron. Our rule (R4) keeps it. Documented in the fixture README.
    let plain = hisat2_sites(&gffread, &script, &d, &[]);
    let (a, b) = report("gffread -T", &ours, &plain);
    assert_eq!(
        a,
        vec!["chrB\t429\t500\t-".to_string()],
        "only the documented geneC intron may be ours-only"
    );
    assert!(
        b.is_empty(),
        "hisat2 found splice sites we do not (see stderr)"
    );
}
