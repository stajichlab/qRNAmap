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
    assert!(
        o.status.success(),
        "stdout:\n{}\nstderr:\n{err}",
        String::from_utf8_lossy(&o.stdout)
    );
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
    fs::write(
        &p,
        format!(
            "sample\tfastq1\tfastq2\nS1\t{}\t{}\n",
            r1.display(),
            r2.display()
        ),
    )
    .unwrap();
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
        m("chrA"),
        m("chrA"),
        m("chrB"),
        m("chrB"),
        m("chrB")
    )
}

fn transcripts(d: &Path, extra: &[&str]) -> PathBuf {
    let (g, gff) = mini();
    let table = d.join("genes.tsv.zst");
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    let mut a = vec![
        "transcripts",
        "--reference",
        s(&g),
        "--gff",
        s(&gff),
        "--out",
        s(&fa),
        "--tx2gene",
        s(&t2g),
        "--genes-table",
        s(&table),
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
    assert_eq!(
        fs::read(&table).unwrap()[..4],
        [0x28, 0xb5, 0x2f, 0xfd],
        "zstd magic number"
    );
}

#[test]
fn transcripts_writes_gene_info_table() {
    let dir = tempfile::tempdir().unwrap();
    let _ = transcripts(dir.path(), &[]);
    let info = dir.path().join("gene_info.tsv.zst");
    assert!(
        info.exists(),
        "gene_info.tsv.zst must exist next to genes.tsv.zst"
    );
    let content = read_table(&info).unwrap();
    assert!(content.starts_with(
        "gene\tsymbol\tgene_id\told_locus_tag\ttranscripts\tproteins\tbiotype\tclass\tcontig\tstrand\tdbxref\tdescription\n"
    ));
}

#[test]
fn mito_contig_flag_marks_genes() {
    let dir = tempfile::tempdir().unwrap();
    let table = transcripts(
        dir.path(),
        &["--mito-contig", "chrB", "--mito-contig", "no_such_contig"],
    );
    assert_eq!(read_table(&table).unwrap(), expected("chrB"));
}

#[test]
fn mito_region_in_gff_marks_genes_without_flag() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let gff = d.join("g.gff3");
    let mut text = String::from("##gff-version 3\nchrB\tmini\tregion\t1\t1500\t.\t+\t.\tID=chrB:1..1500;genome=mitochondrion\n");
    text.push_str(
        &fs::read_to_string(common::mini_gff())
            .unwrap()
            .replace("##gff-version 3\n", ""),
    );
    fs::write(&gff, text).unwrap();
    let table = d.join("genes.tsv");
    let genome = common::mini_genome();
    let (fa, t2g) = (d.join("tx.fa"), d.join("t2g.tsv"));
    ok(&qrnamap(
        &[
            "transcripts",
            "--reference",
            s(&genome),
            "--gff",
            s(&gff),
            "--out",
            s(&fa),
            "--tx2gene",
            s(&t2g),
            "--genes-table",
            s(&table),
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
            "quant",
            "--samples",
            s(&sh),
            "--transcripts",
            s(&d.join("tx.fa.zst")),
            "--tx2gene",
            s(&d.join("t2g.tsv.zst")),
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
            "--gff",
            s(&gff),
            "--mito-contig",
            "chrB",
        ],
        None,
    ));
    assert_eq!(
        read_table(&out.join("genes.tsv.zst")).unwrap(),
        expected("chrB")
    );
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
            "quant",
            "--samples",
            s(&sh),
            "--transcripts",
            s(&d.join("tx.fa.zst")),
            "--tx2gene",
            s(&d.join("t2g.tsv.zst")),
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
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
            "align",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
            "--mito-contig",
            "chrA",
        ],
        None,
    ));
    assert_eq!(
        read_table(&out.join("genes.tsv.zst")).unwrap(),
        expected("chrA")
    );
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
            "run",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--mode",
            "quant",
            "--out",
            s(&out),
            "--step",
            "prepare",
            "--tmpdir",
            s(&d.join("tmp")),
            "--mito-contig",
            "chrA",
        ],
        None,
    ));
    assert_eq!(
        read_table(&out.join("genes.tsv.zst")).unwrap(),
        expected("chrA")
    );
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
            "run",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--mode",
            "quant",
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
        ],
        None,
    ));
    let counts = read_table(&out.join("counts.tsv.zst")).unwrap();
    assert!(counts.starts_with("gene\tS1\n"), "{counts}");
    // `run` takes the row order from the GFF3 (unchanged by this task); counts as in
    // cli_quant_align_matrix_end_to_end.
    assert_eq!(
        counts,
        "gene\tS1\ngeneA\t3000\ngeneB\t2500\ngeneR\t2000\ngeneC\t1500\ngeneU\t1000\n"
    );
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
        (
            "run",
            vec!["--reference", s(&g), "--gff", s(&gff), "--mode", "quant"],
        ),
        (
            "run",
            vec!["--reference", s(&g), "--gff", s(&gff), "--mode", "align"],
        ),
        (
            "quant",
            vec![
                "--transcripts",
                s(&fa),
                "--tx2gene",
                s(&t2g),
                "--gff",
                s(&gff),
            ],
        ),
        ("align", vec!["--reference", s(&g), "--gff", s(&gff)]),
    ];
    for (i, (cmd, extra)) in cases.into_iter().enumerate() {
        let out = d.join(format!("out{i}"));
        fs::create_dir_all(&out).unwrap();
        let stale = out.join("genes.tsv.zst");
        fs::write(&stale, b"stale").unwrap();
        let stale_partial = out.join("genes.tsv.zst.partial");
        fs::write(&stale_partial, b"stale").unwrap();
        let mut a = vec![
            cmd,
            "--samples",
            s(&sh),
            "--out",
            s(&out),
            "--tmpdir",
            s(&tmp),
        ];
        a.extend(extra);
        fail(&qrnamap(&a, Some(&empty)));
        assert!(
            !stale.exists(),
            "case {i} ({cmd}): the stale genes table must be gone"
        );
        assert!(
            !stale_partial.exists(),
            "case {i} ({cmd}): the stale partial must be gone"
        );
    }
}

fn gene(label: &str) -> Gene {
    Gene {
        label: label.into(),
        contig: "chrA".into(),
        strand: Strand::Plus,
        biotype: None,
        symbol: None,
        gene_id: None,
        dbxref: vec![],
        old_locus_tag: None,
        description: None,
        transcripts: vec![Transcript {
            id: "t1".into(),
            kind: "mRNA".into(),
            exons: vec![(1, 10)],
            cds: vec![],
            utr5: vec![],
            utr3: vec![],
            proteins: vec![],
        }],
    }
}

#[test]
fn genes_table_label_with_tab_is_refused() {
    let d = tempfile::tempdir().unwrap();
    for bad in ["a\tb", "a\nb", "a\rb"] {
        let model = GeneModel {
            genes: vec![gene("fine"), gene(bad)],
            ..GeneModel::default()
        };
        let path = d.path().join("genes.tsv.zst");
        let e = write_genes_table(&path, &model, &[]).unwrap_err();
        assert!(
            matches!(e, PipelineError::BadGeneLabel(ref l) if l == bad),
            "{e}"
        );
        assert!(
            e.to_string().contains("tab, newline or carriage return"),
            "{e}"
        );
        assert!(!path.exists(), "no table");
        assert!(
            !d.path().join("genes.tsv.zst.partial").exists(),
            "no partial file"
        );
    }
}

#[test]
fn write_genes_table_suffix_picks_the_compression() {
    let d = tempfile::tempdir().unwrap();
    let model = GeneModel {
        genes: vec![gene("g1")],
        ..GeneModel::default()
    };
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
            "run",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--mode",
            "quant",
            "--out",
            s(&out),
            "--slurm",
            "--qrnamap-bin",
            env!("CARGO_BIN_EXE_qrnamap"),
            "--mito-contig",
            hostile,
            "--mito-contig",
            "chrM",
        ],
        None,
    ));
    let sdir = out.join("slurm");
    let prepare = fs::read_to_string(sdir.join("prepare.sbatch")).unwrap();
    assert!(
        prepare.contains("'--mito-contig' 'my contig$HOME' '--mito-contig' 'chrM'"),
        "{prepare}"
    );
    for f in [
        "prepare.sbatch",
        "samples.sbatch",
        "final.sbatch",
        "submit.sh",
    ] {
        let text = fs::read_to_string(sdir.join(f)).unwrap();
        assert!(!text.contains("ZzSample"), "{f} holds a sample name");
        let o = Command::new("bash")
            .arg("-n")
            .arg(sdir.join(f))
            .output()
            .unwrap();
        assert!(
            o.status.success(),
            "bash -n {f}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
    }
    // The unquoted `$HOME` must not appear outside single quotes.
    assert!(
        !prepare.replace("'my contig$HOME'", "").contains("$HOME"),
        "{prepare}"
    );
}

// ---------------------------------------------------------------------------------------------
// Final review fix wave (R35)

/// The mini GFF3 without geneU (and its transcript tU1).
fn gff_without_gene_u(d: &Path) -> PathBuf {
    let text: String = fs::read_to_string(common::mini_gff())
        .unwrap()
        .lines()
        .filter(|l| !l.contains("geneU") && !l.contains("tU1"))
        .map(|l| format!("{l}\n"))
        .collect();
    let p = d.join("no_u.gff3");
    fs::write(&p, text).unwrap();
    p
}

/// C1 (R35-1), the reviewer's scenario: a finished quant run with `--mito-contig chrB`, then
/// `strand-check` into the same directory. The gene table must keep its bytes.
#[test]
fn strand_check_leaves_existing_genes_table_untouched() {
    if !common::require_tools(&[
        "kallisto",
        "hisat2-build",
        "hisat2-align-s",
        "featureCounts",
    ]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let out = d.join("R");
    let tmp = d.join("tmp");
    ok(&qrnamap(
        &[
            "run",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--mode",
            "quant",
            "--out",
            s(&out),
            "--tmpdir",
            s(&tmp),
            "--mito-contig",
            "chrB",
        ],
        None,
    ));
    let table = out.join("genes.tsv.zst");
    let before = fs::read(&table).unwrap();
    assert_eq!(read_table(&table).unwrap(), expected("chrB"));
    ok(&qrnamap(
        &[
            "strand-check",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--out",
            s(&out),
            "--tmpdir",
            s(&tmp),
        ],
        None,
    ));
    assert!(out.join("strand_check.tsv").is_file());
    assert_eq!(
        fs::read(&table).unwrap(),
        before,
        "strand-check changed genes.tsv.zst"
    );
    assert!(!out.join("genes.tsv.zst.partial").exists());
}

#[test]
fn strand_check_does_not_create_genes_table_in_a_fresh_dir() {
    if !common::require_tools(&["hisat2-build", "hisat2-align-s", "featureCounts"]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let out = d.join("sc");
    ok(&qrnamap(
        &[
            "strand-check",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
        ],
        None,
    ));
    assert!(out.join("strand_check.tsv").is_file());
    assert!(
        !out.join("genes.tsv.zst").exists(),
        "strand-check wrote a gene table"
    );
    assert!(!out.join("genes.tsv.zst.partial").exists());
}

fn tx2gene_map(t2g: &Path) -> std::collections::HashMap<String, String> {
    qrnamap_cli::pipeline::read_tx2gene(t2g).unwrap()
}

fn model_of(gff: &Path) -> GeneModel {
    qseq_io::parse_gff3(std::io::BufReader::new(fs::File::open(gff).unwrap())).unwrap()
}

/// I4 (R35-3), tool-free: the check itself, and the CLI with an empty PATH (the error comes
/// before any tool and before the table is written).
#[test]
fn quant_gff_not_matching_tx2gene_is_named_error() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (g, gff) = mini();
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    qrnamap_cli::pipeline::run_transcripts(&g, &gff, &fa, &t2g, None, None).unwrap();
    let no_u = gff_without_gene_u(d);
    let e = qrnamap_cli::pipeline::check_tx2gene_in_model(&tx2gene_map(&t2g), &model_of(&no_u))
        .unwrap_err();
    assert!(
        matches!(&e, PipelineError::GenesTableMismatch { missing, n_missing: 1 } if missing == &["geneU".to_string()]),
        "{e:?}"
    );
    assert!(e.to_string().contains("geneU"), "{e}");

    // A GFF3 with none of the genes: the list is bounded, the count is the full one.
    let other = d.join("other.gff3");
    let mut text = String::from("##gff-version 3\n");
    for i in 0..30 {
        text.push_str(&format!("chrA\tx\tgene\t1\t10\t.\t+\t.\tID=zz{i:02}\n"));
    }
    fs::write(&other, text).unwrap();
    let e = qrnamap_cli::pipeline::check_tx2gene_in_model(&tx2gene_map(&t2g), &model_of(&other))
        .unwrap_err();
    match &e {
        PipelineError::GenesTableMismatch { missing, n_missing } => {
            assert_eq!(*n_missing, 5);
            assert_eq!(missing, &["geneA", "geneB", "geneC", "geneR", "geneU"]);
        }
        _ => panic!("{e:?}"),
    }

    // CLI: plain quant with the GFF3 that lacks geneU; no tool on PATH.
    let sh = sheet(d);
    let out = d.join("q");
    fs::create_dir_all(&out).unwrap();
    fs::write(out.join("genes.tsv.zst"), b"stale").unwrap();
    let empty = d.join("empty_path");
    fs::create_dir_all(&empty).unwrap();
    let err = fail(&qrnamap(
        &[
            "quant",
            "--samples",
            s(&sh),
            "--transcripts",
            s(&fa),
            "--tx2gene",
            s(&t2g),
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
            "--gff",
            s(&no_u),
        ],
        Some(&empty),
    ));
    assert!(
        err.contains("1 gene(s) of the tx2gene table are not genes of the --gff GFF3: geneU"),
        "{err}"
    );
    assert!(
        !err.contains("kallisto"),
        "the check must come before any tool: {err}"
    );
    assert!(
        !out.join("genes.tsv.zst").exists(),
        "no gene table after the error"
    );
    assert!(!out.join("transcripts.kidx").exists());
}

#[test]
fn quant_gff_superset_of_tx2gene_is_ok() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (g, gff) = mini();
    let no_u = gff_without_gene_u(d);
    let (fa, t2g) = (d.join("tx.fa.zst"), d.join("t2g.tsv.zst"));
    // tx2gene without geneU; the full GFF3 has one gene more.
    qrnamap_cli::pipeline::run_transcripts(&g, &no_u, &fa, &t2g, None, None).unwrap();
    qrnamap_cli::pipeline::check_tx2gene_in_model(&tx2gene_map(&t2g), &model_of(&gff)).unwrap();
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let sh = sheet(d);
    let out = d.join("q");
    ok(&qrnamap(
        &[
            "quant",
            "--samples",
            s(&sh),
            "--transcripts",
            s(&fa),
            "--tx2gene",
            s(&t2g),
            "--out",
            s(&out),
            "--tmpdir",
            s(&d.join("tmp")),
            "--gff",
            s(&gff),
        ],
        None,
    ));
    assert_eq!(
        read_table(&out.join("genes.tsv.zst")).unwrap(),
        expected("none")
    );
}

/// R35-7: a `--mito-contig` name that matches no contig of the GFF3 gives a note, not an error.
#[test]
fn mito_contig_matching_no_contig_prints_note() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (g, gff) = mini();
    let table = d.join("genes.tsv");
    let o = qrnamap(
        &[
            "transcripts",
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--out",
            s(&d.join("tx.fa")),
            "--tx2gene",
            s(&d.join("t2g.tsv")),
            "--genes-table",
            s(&table),
            "--mito-contig",
            "chrB",
            "--mito-contig",
            "chrb",
        ],
        None,
    );
    let err = ok(&o);
    assert!(
        err.contains("note: --mito-contig 'chrb' matches no contig of the GFF3"),
        "{err}"
    );
    assert!(!err.contains("'chrB' matches"), "chrB exists: {err}");
    assert_eq!(fs::read_to_string(&table).unwrap(), expected("chrB"));
}

#[test]
fn transcripts_mito_contig_without_genes_table_prints_note() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (g, gff) = mini();
    let err = ok(&qrnamap(
        &[
            "transcripts",
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--out",
            s(&d.join("tx.fa")),
            "--tx2gene",
            s(&d.join("t2g.tsv")),
            "--mito-contig",
            "chrB",
        ],
        None,
    ));
    assert!(
        err.contains("note: --mito-contig has no effect without --genes-table"),
        "{err}"
    );
}

#[test]
fn quant_task_id_mito_contig_prints_note() {
    if common::require_tool("kallisto").is_none() {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sh = sheet(d);
    let (g, gff) = mini();
    let out = d.join("run");
    let tmp = d.join("tmp");
    ok(&qrnamap(
        &[
            "run",
            "--samples",
            s(&sh),
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--mode",
            "quant",
            "--out",
            s(&out),
            "--step",
            "prepare",
            "--tmpdir",
            s(&tmp),
        ],
        None,
    ));
    let before = fs::read(out.join("genes.tsv.zst")).unwrap();
    let err = ok(&qrnamap(
        &[
            "quant",
            "--samples",
            s(&sh),
            "--transcripts",
            s(&out.join("transcripts.fa.zst")),
            "--tx2gene",
            s(&out.join("tx2gene.tsv.zst")),
            "--out",
            s(&out),
            "--tmpdir",
            s(&tmp),
            "--task-id",
            "0",
            "--reference",
            s(&g),
            "--gff",
            s(&gff),
            "--mito-contig",
            "chrB",
        ],
        None,
    ));
    assert!(
        err.contains("note: --mito-contig has no effect with --task-id"),
        "{err}"
    );
    assert_eq!(fs::read(out.join("genes.tsv.zst")).unwrap(), before);
}

#[test]
fn real_genomes_gene_info_validation() {
    use flate2::read::GzDecoder;
    use qrnamap_core::{build_transcripts, gene_info_tsv, tx2gene_tsv};
    use qseq_io::{parse_gff3, read_fasta};
    use std::collections::HashSet;
    use std::io::BufReader;

    let genomes = [
        (
            "yeast",
            "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures/genomes/yeast/GCF_000146045.2_R64_genomic.fna.gz",
            "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures/genomes/yeast/GCF_000146045.2_R64_genomic.gff.gz",
        ),
        (
            "afum",
            "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures/genomes/afum/GCF_000002655.1_ASM265v1_genomic.fna.gz",
            "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures/genomes/afum/GCF_000002655.1_ASM265v1_genomic.gff.gz",
        ),
        (
            "h99",
            "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/page-v1/dl_h99/GCF_000149245.1_CNA3_genomic.fna.gz",
            "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/page-v1/dl_h99/GCF_000149245.1_CNA3_genomic.gff.gz",
        ),
    ];

    for (name, fna_path, gff_path) in genomes {
        if !Path::new(fna_path).exists() || !Path::new(gff_path).exists() {
            eprintln!("skip {name}: files not found");
            continue;
        }
        let gff_file = fs::File::open(gff_path).unwrap();
        let model = parse_gff3(BufReader::new(GzDecoder::new(gff_file))).unwrap();
        let fna_file = fs::File::open(fna_path).unwrap();
        let genome = read_fasta(BufReader::new(GzDecoder::new(fna_file))).unwrap();
        let records = build_transcripts(&genome, &model).unwrap();
        let tx2gene = tx2gene_tsv(&records);

        let tsv = gene_info_tsv(&model, &[]);
        let lines: Vec<&str> = tsv.lines().collect();
        assert!(!lines.is_empty(), "{name}: gene_info.tsv must not be empty");
        assert_eq!(
            lines[0],
            "gene\tsymbol\tgene_id\told_locus_tag\ttranscripts\tproteins\tbiotype\tclass\tcontig\tstrand\tdbxref"
        );

        let mut n_symbol = 0;
        let mut n_gene_id = 0;
        let mut n_old_locus_tag = 0;
        let mut n_proteins = 0;
        let mut n_dbxref = 0;
        let mut n_fungidb = 0;
        let mut all_txs_in_info = HashSet::new();

        for line in &lines[1..] {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 11, "11 columns required on {line}");
            let _gene = cols[0];
            let symbol = cols[1];
            let gene_id = cols[2];
            let old_locus_tag = cols[3];
            let transcripts = cols[4];
            let proteins = cols[5];
            let dbxref = cols[10];

            if !symbol.is_empty() {
                n_symbol += 1;
            }
            if !gene_id.is_empty() {
                n_gene_id += 1;
            }
            if !old_locus_tag.is_empty() {
                n_old_locus_tag += 1;
            }
            if !proteins.is_empty() {
                n_proteins += 1;
            }
            if !dbxref.is_empty() {
                n_dbxref += 1;
            }
            if dbxref.split(',').any(|x| x.starts_with("FungiDB:")) {
                n_fungidb += 1;
            }
            if !transcripts.is_empty() {
                for tx in transcripts.split(',') {
                    all_txs_in_info.insert(tx.to_string());
                }
            }
        }

        eprintln!(
            "{name}: total_genes={} symbol={} gene_id={} old_locus_tag={} proteins={} dbxref={}",
            lines.len() - 1,
            n_symbol,
            n_gene_id,
            n_old_locus_tag,
            n_proteins,
            n_dbxref
        );

        for line in tx2gene.lines() {
            let tx = line.split('\t').next().unwrap();
            assert!(
                all_txs_in_info.contains(tx),
                "{name}: transcript '{tx}' from tx2gene must appear in gene_info.tsv"
            );
        }

        assert!(n_symbol > 0, "{name}: must have symbols");
        assert!(n_gene_id > 0, "{name}: must have gene_ids");
        assert!(n_proteins > 0, "{name}: must have proteins");
        assert!(n_dbxref > 0, "{name}: must have dbxref entries");
        if name == "afum" {
            assert!(n_old_locus_tag > 0, "afum must have old_locus_tags");
            // Af293 (GCF_000002655.1) carries Dbxref=FungiDB:AFUA_1G00100,GeneID:... on gene features;
            // the results-page link-outs are built from this column.
            assert!(
                n_fungidb > 0,
                "afum (Af293) must have FungiDB dbxref entries"
            );
        }
    }

    // Every locus name in counts.tsv (expected_quant.tsv for the fixture) appears in gene_info.tsv
    let mito_gff = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/browser-page/fixtures/mito.gff3");
    let exp_quant = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/browser-page/fixtures/expected_quant.tsv");
    if mito_gff.exists() && exp_quant.exists() {
        let model = parse_gff3(BufReader::new(fs::File::open(&mito_gff).unwrap())).unwrap();
        let tsv = gene_info_tsv(&model, &[]);
        let info_genes: HashSet<&str> = tsv
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        let quant_text = fs::read_to_string(&exp_quant).unwrap();
        for line in quant_text.lines().skip(1) {
            let gene = line.split('\t').nth(1).unwrap();
            assert!(
                info_genes.contains(gene),
                "locus '{gene}' from counts must appear in gene_info.tsv"
            );
        }
    }
}
