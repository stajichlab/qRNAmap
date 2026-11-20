//! `genes.tsv.zst`: class, biotype and mitochondrial flag per gene, written next to the counts.

mod common;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use common::{records, sim, COUNTS};
use qrnamap_cli::pipeline::{read_table, write_genes_table, PipelineError};
use qrnamap_core::Strandedness;
use qseq_io::{Gene, GeneModel, Strand, Transcript};

const HEADER: &str = "gene\tcontig\tstrand\tbiotype\tclass\tmito\tn_transcripts\n";

fn qrnamap(args: &[&str], path_env: Option<&Path>) -> Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_qrnamap"));
    c.args(args).env_remove("SCRATCH").env_remove("TMPDIR");
    if let Some(p) = path_env {
        c.env("PATH", p);
    }
    c.output().unwrap()
}

fn ok(o: &Output) -> String {
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    assert!(o.status.success(), "stdout:\n{}\nstderr:\n{err}", String::from_utf8_lossy(&o.stdout));
    err
}

fn fail(o: &Output) -> String {
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    assert!(!o.status.success(), "expected failure; stderr:\n{err}");
    err
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

fn mini() -> (PathBuf, PathBuf) {
    (common::mini_genome(), common::mini_gff())
}

/// One sample S1 with simulated reads.
fn sheet(d: &Path) -> PathBuf {
    let pairs = sim::simulate(&records(), &COUNTS, Strandedness::Unstranded, 41);
    let (r1, r2) = (d.join("s_1.fq"), d.join("s_2.fq"));
    sim::write_fastq(&pairs, &r1, &r2);
    let p = d.join("sheet.tsv");
    fs::write(&p, format!("sample\tfastq1\tfastq2\nS1\t{}\t{}\n", r1.display(), r2.display())).unwrap();
    p
}

fn expected(mito_chr: &str) -> String {
    let m = |c: &str| if c == mito_chr { "yes" } else { "no" };
    format!(
        "{HEADER}geneA\tchrA\t+\tmRNA\tprotein_coding\t{}\t1\n\
         geneB\tchrA\t-\tmRNA\tprotein_coding\t{}\t2\n\
         geneR\tchrB\t+\trRNA\trRNA\t{}\t1\n\
         geneC\tchrB\t-\tmRNA\tprotein_coding\t{}\t1\n\
         geneU\tchrB\t+\tmRNA\tprotein_coding\t{}\t1\n",
        m("chrA"), m("chrA"), m("chrB"), m("chrB"), m("chrB")
    )
}

fn transcripts(d: &Path, extra: &[&str]) -> PathBuf {
    let (g, gff) = mini();
    let table = d.join("genes.tsv.zst");
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    let mut a = vec![
        "transcripts", "--reference", s(&g), "--gff", s(&gff), "--out", s(&fa), "--tx2gene", s(&t2g),
        "--genes-table", s(&table),
    ];
    a.extend(extra);
    ok(&qrnamap(&a, None));
    table
}

#[test]
fn transcripts_writes_genes_table() {
    let dir = tempfile::tempdir().unwrap();
    let table = transcripts(dir.path(), &[]);
    // read_table decompresses by suffix (zstd frame), so this also proves the file is .zst.
    assert_eq!(read_table(&table).unwrap(), expected("none"));
    assert_eq!(fs::read(&table).unwrap()[..4], [0x28, 0xb5, 0x2f, 0xfd], "zstd magic number");
}

#[test]
fn mito_contig_flag_marks_genes() {
    let dir = tempfile::tempdir().unwrap();
    let table = transcripts(dir.path(), &["--mito-contig", "chrB", "--mito-contig", "no_such_contig"]);
    assert_eq!(read_table(&table).unwrap(), expected("chrB"));
}

#[test]
fn mito_region_in_gff_marks_genes_without_flag() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let gff = d.join("g.gff3");
    let mut text = String::from("##gff-version 3\nchrB\tmini\tregion\t1\t1500\t.\t+\t.\tID=chrB:1..1500;genome=mitochondrion\n");
    text.push_str(&fs::read_to_string(common::mini_gff()).unwrap().replace("##gff-version 3\n", ""));
    fs::write(&gff, text).unwrap();
    let table = d.join("genes.tsv");
    let genome = common::mini_genome();
    let (fa, t2g) = (d.join("tx.fa"), d.join("t2g.tsv"));
    ok(&qrnamap(
        &[
            "transcripts", "--reference", s(&genome), "--gff", s(&gff), "--out", s(&fa),
            "--tx2gene", s(&t2g), "--genes-table", s(&table),
        ],
        None,
    ));
    // A plain name: no compression.
    assert_eq!(fs::read_to_string(&table).unwrap(), expected("chrB"));
}

#[test]
fn plain_quant_writes_genes_table_in_out() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (_, gff) = mini();
    transcripts(d, &[]);
    fs::remove_file(d.join("genes.tsv.zst")).unwrap();
    let out = d.join("q");
    ok(&qrnamap(
        &[
            "quant", "--samples", s(&sh), "--transcripts", s(&d.join("tx.fa.zst")), "--tx2gene", s(&d.join("t2g.tsv.zst")),
            "--out", s(&out), "--tmpdir", s(&d.join("tmp")), "--gff", s(&gff), "--mito-contig", "chrB",
        ],
        None,
    ));
    assert_eq!(read_table(&out.join("genes.tsv.zst")).unwrap(), expected("chrB"));
}

#[test]
fn plain_quant_without_gff_writes_no_table_and_says_so() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    transcripts(d, &[]);
    let out = d.join("q");
    let err = ok(&qrnamap(
        &[
            "quant", "--samples", s(&sh), "--transcripts", s(&d.join("tx.fa.zst")), "--tx2gene", s(&d.join("t2g.tsv.zst")),
            "--out", s(&out), "--tmpdir", s(&d.join("tmp")),
        ],
        None,
    ));
    assert!(err.contains("genes.tsv.zst is not written"), "{err}");
    assert!(!out.join("genes.tsv.zst").exists());
}

#[test]
fn plain_align_writes_genes_table_in_out() {
    if !common::require_tools(&["hisat2-build", "hisat2-align-s", "featureCounts"]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let out = d.join("a");
    ok(&qrnamap(
        &[
            "align", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--out", s(&out), "--tmpdir",
            s(&d.join("tmp")), "--mito-contig", "chrA",
        ],
        None,
    ));
    assert_eq!(read_table(&out.join("genes.tsv.zst")).unwrap(), expected("chrA"));
}

#[test]
fn run_prepare_writes_genes_table() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    ok(&qrnamap(
        &[
            "run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out),
            "--step", "prepare", "--tmpdir", s(&d.join("tmp")), "--mito-contig", "chrA",
        ],
        None,
    ));
    assert_eq!(read_table(&out.join("genes.tsv.zst")).unwrap(), expected("chrA"));
    assert!(out.join(".prepared").is_file());
}

#[test]
fn counts_matrix_is_unchanged() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    ok(&qrnamap(
        &[
            "run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out),
            "--tmpdir", s(&d.join("tmp")),
        ],
        None,
    ));
    let counts = read_table(&out.join("counts.tsv.zst")).unwrap();
    assert!(counts.starts_with("gene\tS1\n"), "{counts}");
    // `run` takes the row order from the GFF3 (unchanged by this task); counts as in
    // cli_quant_align_matrix_end_to_end.
    assert_eq!(counts, "gene\tS1\ngeneA\t3000\ngeneB\t2500\ngeneR\t2000\ngeneC\t1500\ngeneU\t1000\n");
    assert!(out.join("genes.tsv.zst").is_file());
}

/// Run with an empty PATH, so the index build fails after the old outputs were removed.
#[test]
fn old_genes_table_is_removed_before_a_failing_run() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let empty = d.join("empty_path");
    fs::create_dir_all(&empty).unwrap();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    let tmp = d.join("tmp");
    let cases: Vec<(&str, Vec<&str>)> = vec![
        ("run", vec!["--reference", s(&g), "--gff", s(&gff), "--mode", "quant"]),
        ("run", vec!["--reference", s(&g), "--gff", s(&gff), "--mode", "align"]),
        ("quant", vec!["--transcripts", s(&fa), "--tx2gene", s(&t2g), "--gff", s(&gff)]),
        ("align", vec!["--reference", s(&g), "--gff", s(&gff)]),
    ];
    for (i, (cmd, extra)) in cases.into_iter().enumerate() {
        let out = d.join(format!("out{i}"));
        fs::create_dir_all(&out).unwrap();
        let stale = out.join("genes.tsv.zst");
        fs::write(&stale, b"stale").unwrap();
        let stale_partial = out.join("genes.tsv.zst.partial");
        fs::write(&stale_partial, b"stale").unwrap();
        let mut a = vec![cmd, "--samples", s(&sh), "--out", s(&out), "--tmpdir", s(&tmp)];
        a.extend(extra);
        fail(&qrnamap(&a, Some(&empty)));
        assert!(!stale.exists(), "case {i} ({cmd}): the stale genes table must be gone");
        assert!(!stale_partial.exists(), "case {i} ({cmd}): the stale partial must be gone");
    }
}

fn gene(label: &str) -> Gene {
    Gene {
        label: label.into(),
        contig: "chrA".into(),
        strand: Strand::Plus,
        biotype: None,
        transcripts: vec![Transcript {
            id: "t1".into(),
            kind: "mRNA".into(),
            exons: vec![(1, 10)],
            cds: vec![],
            utr5: vec![],
            utr3: vec![],
        }],
    }
}

#[test]
fn genes_table_label_with_tab_is_refused() {
    let d = tempfile::tempdir().unwrap();
    for bad in ["a\tb", "a\nb", "a\rb"] {
        let model = GeneModel { genes: vec![gene("fine"), gene(bad)], ..GeneModel::default() };
        let path = d.path().join("genes.tsv.zst");
        let e = write_genes_table(&path, &model, &[]).unwrap_err();
        assert!(matches!(e, PipelineError::BadGeneLabel(ref l) if l == bad), "{e}");
        assert!(e.to_string().contains("tab, newline or carriage return"), "{e}");
        assert!(!path.exists(), "no table");
        assert!(!d.path().join("genes.tsv.zst.partial").exists(), "no partial file");
    }
}

#[test]
fn write_genes_table_suffix_picks_the_compression() {
    let d = tempfile::tempdir().unwrap();
    let model = GeneModel { genes: vec![gene("g1")], ..GeneModel::default() };
    let want = format!("{HEADER}g1\tchrA\t+\tmRNA\tprotein_coding\tno\t1\n");
    for name in ["g.tsv", "g.tsv.gz", "g.tsv.zst"] {
        let p = d.path().join(name);
        write_genes_table(&p, &model, &[]).unwrap();
        assert_eq!(read_table(&p).unwrap(), want, "{name}");
    }
    assert_eq!(fs::read(d.path().join("g.tsv")).unwrap(), want.as_bytes());
}

#[test]
fn run_slurm_forwards_mito_contig_quoted() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = d.join("sheet.tsv");
    fs::write(&sh, "sample\tfastq1\tfastq2\nZzSampleOne\t/x/a_1.fq\t/x/a_2.fq\nZzSampleTwo\t/x/b_1.fq\t/x/b_2.fq\n").unwrap();
    let (g, gff) = mini();
    let out = d.join("cl");
    let hostile = "my contig$HOME";
    ok(&qrnamap(
        &[
            "run", "--samples", s(&sh), "--reference", s(&g), "--gff", s(&gff), "--mode", "quant", "--out", s(&out), "--slurm",
            "--qrnamap-bin", env!("CARGO_BIN_EXE_qrnamap"), "--mito-contig", hostile, "--mito-contig", "chrM",
        ],
        None,
    ));
    let sdir = out.join("slurm");
    let prepare = fs::read_to_string(sdir.join("prepare.sbatch")).unwrap();
    assert!(prepare.contains("'--mito-contig' 'my contig$HOME' '--mito-contig' 'chrM'"), "{prepare}");
    for f in ["prepare.sbatch", "samples.sbatch", "final.sbatch", "submit.sh"] {
        let text = fs::read_to_string(sdir.join(f)).unwrap();
        assert!(!text.contains("ZzSample"), "{f} holds a sample name");
        let o = Command::new("bash").arg("-n").arg(sdir.join(f)).output().unwrap();
        assert!(o.status.success(), "bash -n {f}: {}", String::from_utf8_lossy(&o.stderr));
    }
    // The unquoted `$HOME` must not appear outside single quotes.
    assert!(!prepare.replace("'my contig$HOME'", "").contains("$HOME"), "{prepare}");
}
