//! The wrapper functions, called natively (no WebAssembly engine needed).

use qrnamap_core::{kallisto_index, kallisto_quant, RunConfig, Sample, Strandedness};
use qrnamap_wasm::*;
use serde_json::Value;
use std::collections::BTreeMap;

fn json(s: &str) -> Value {
    serde_json::from_str(s).unwrap()
}

fn args(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|a| a.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn catalog_lists_genomes_with_a_gff() {
    let v = json(&catalog_json());
    let a = v.as_array().unwrap();
    assert_eq!(a.len(), genome_catalog::builtin().usable_for_rna().len());
    assert!(!a.is_empty());
    let g = &a[0];
    for k in [
        "accession",
        "label",
        "organism",
        "genome_size",
        "fasta",
        "gff",
    ] {
        assert!(!g[k].is_null(), "{k}");
    }
    assert!(g["fasta"]["url"].as_str().unwrap().starts_with("https://"));
}

#[test]
fn kallisto_index_args_equal_toolcall() {
    let got = kallisto_args_native(r#"{"kind":"index","fasta":"t.fa","index":"t.idx"}"#).unwrap();
    assert_eq!(args(&json(&got)), kallisto_index("t.fa", "t.idx").args);
}

#[test]
fn kallisto_quant_args_equal_toolcall() {
    let sample = Sample {
        name: "S1".into(),
        fastq1: vec!["a_1.fq".into(), "b_1.fq".into()],
        fastq2: vec!["a_2.fq".into(), "b_2.fq".into()],
        covariates: BTreeMap::new(),
    };
    let cfg = RunConfig {
        threads: 1,
        strandedness: Strandedness::Reverse,
        ..RunConfig::default()
    };
    let want = kallisto_quant(&sample, "t.idx", "out", &cfg).unwrap().args;
    let req = r#"{"kind":"quant","sample":{"name":"S1","fastq1":["a_1.fq","b_1.fq"],"fastq2":["a_2.fq","b_2.fq"]},
        "index":"t.idx","out":"out","threads":1,"strandedness":"reverse"}"#;
    assert_eq!(args(&json(&kallisto_args_native(req).unwrap())), want);

    // Single-end needs a fragment length, as in toolcall.
    let se =
        r#"{"kind":"quant","sample":{"name":"S1","fastq1":["a.fq"]},"index":"t.idx","out":"out"}"#;
    assert!(kallisto_args_native(se).is_err());
    let se = r#"{"kind":"quant","sample":{"name":"S1","fastq1":["a.fq"]},"index":"t.idx","out":"out","fragment_len":[200,20]}"#;
    let sample = Sample {
        name: "S1".into(),
        fastq1: vec!["a.fq".into()],
        fastq2: vec![],
        covariates: BTreeMap::new(),
    };
    let cfg = RunConfig {
        fragment_len: Some((200.0, 20.0)),
        ..RunConfig::default()
    };
    let want = kallisto_quant(&sample, "t.idx", "out", &cfg).unwrap().args;
    assert_eq!(args(&json(&kallisto_args_native(se).unwrap())), want);
}

#[test]
fn kallisto_args_reject_bad_requests() {
    assert!(kallisto_args_native("not json").is_err());
    assert!(kallisto_args_native(r#"{"kind":"salmon"}"#).is_err());
    assert!(kallisto_args_native(r#"{"kind":"index"}"#).is_err());
}

#[test]
fn de_command_quotes_words() {
    assert_eq!(
        de_command_native("samples.tsv", "counts_preview.tsv", "~ condition", "condition,AmpB,control"),
        "qrnamap de --counts counts_preview.tsv --samples samples.tsv --design '~ condition' --contrast condition,AmpB,control --out de"
    );
    let c = de_command_native("a b.tsv", "c.tsv", "~ x", "x,it's,ref");
    assert!(c.contains("'a b.tsv'"), "{c}");
    assert!(c.contains(r"'x,it'\''s,ref'"), "{c}");
}

#[test]
fn validate_sheet_reports_errors_and_samples() {
    let ok = json(&validate_sheet_native(
        "sample\tfastq1\tfastq2\tcondition\nA\ta_1.fq\ta_2.fq\tctl\nB\tb_1.fq\tb_2.fq\ttrt\n",
    ));
    assert_eq!(ok["ok"], true);
    assert_eq!(ok["samples"], serde_json::json!(["A", "B"]));
    assert_eq!(ok["columns"], serde_json::json!([["condition", "factor"]]));
    let bad = json(&validate_sheet_native("sample\tfastq1\nA\ta.fq\nA\tb.fq\n"));
    assert_eq!(bad["ok"], false);
    assert!(bad["errors"][0]
        .as_str()
        .unwrap()
        .contains("duplicate sample name"));
    let empty = json(&validate_sheet_native(""));
    assert_eq!(empty["ok"], false);
}

#[test]
fn validate_sheet_gives_rows_and_refuses_reserved_names() {
    let ok = json(&validate_sheet_native(
        "# comment\r\nsample\tfastq1\tfastq2\tcondition\tdose:numeric\r\n\r\nA\tdir/a 1.fq.gz\t\tctl\t1.5\r\nB\tb1.fq,b2.fq\tc1.fq,c2.fq\ttrt\t2\r\n",
    ));
    assert_eq!(ok["ok"], true, "{ok}");
    assert_eq!(
        ok["rows"],
        serde_json::json!([
            {"sample": "A", "fastq1": ["dir/a 1.fq.gz"], "fastq2": [], "covariates": {"condition": "ctl", "dose": "1.5"}},
            {"sample": "B", "fastq1": ["b1.fq", "b2.fq"], "fastq2": ["c1.fq", "c2.fq"], "covariates": {"condition": "trt", "dose": "2"}},
        ])
    );
    assert_eq!(
        ok["columns"],
        serde_json::json!([["condition", "factor"], ["dose", "numeric"]])
    );
    for name in ["index", "counts.tsv.zst", "slurm.partial"] {
        let bad = json(&validate_sheet_native(&format!(
            "sample\tfastq1\nA\ta.fq\n{name}\tb.fq\n"
        )));
        assert_eq!(bad["ok"], false, "{name}");
        assert_eq!(bad["rows"], serde_json::json!([]));
        assert_eq!(
            bad["errors"][0],
            format!("sample '{name}': the name is reserved: qrnamap writes '{name}' in the output directory. Rename the sample in the sheet")
        );
    }
}

#[test]
fn matrix_and_gene_counts_through_the_wrapper() {
    let t2g = "tx1\tgA\ntx2\tgA\ntx3\tgB\n";
    let ab = "target_id\tlength\teff_length\test_counts\ttpm\ntx1\t1\t1\t1.4\t0\ntx2\t1\t1\t1.4\t0\ntx3\t1\t1\t0\t0\n";
    let gc = json(&gene_counts_native(ab, t2g).unwrap());
    assert_eq!(gc["counts"], serde_json::json!({"gA": 3, "gB": 0}));
    assert_eq!(gc["total"], 3);
    let input = serde_json::json!({
        "tx2gene": t2g,
        "samples": [{"name": "S1", "counts": gc["counts"]}, {"name": "S2", "counts": {"gB": 5}}],
    });
    let m = json(&assemble_matrix_native(&input.to_string()).unwrap());
    assert_eq!(m["tsv"], "gene\tS1\tS2\ngA\t3\t0\ngB\t0\t5\n");
    // Duplicate sample names are refused.
    let dup = serde_json::json!({"tx2gene": t2g, "samples": [{"name": "S1", "counts": {}}, {"name": "S1", "counts": {}}]});
    assert!(assemble_matrix_native(&dup.to_string())
        .unwrap_err()
        .contains("duplicate"));
    // An unknown transcript is an error naming it.
    let e = gene_counts_native("target_id\test_counts\nzz\t1\n", t2g).unwrap_err();
    assert!(e.contains("zz"), "{e}");
}

#[test]
fn qc_small_designs_give_messages_not_failures() {
    let one = json(&qc_native("gene\tS1\ng1\t5\ng2\t7\n").unwrap());
    assert!(one["correlation"].is_null());
    assert!(one["correlation_error"]
        .as_str()
        .unwrap()
        .contains("at least 2"));
    assert!(one["pca"].is_null());
    let two = json(&qc_native("gene\tS1\tS2\ng1\t5\t9\ng2\t7\t1\n").unwrap());
    assert!(two["correlation"].is_array());
    assert!(two["pca"].is_null());
    assert!(two["pca_error"].as_str().unwrap().contains("at least 3"));
    let e = qc_native("gene\tS1\tS2\ng1\t5\t0\n").unwrap_err();
    assert!(e.contains("S2"), "{e}");
    assert!(qc_native("").is_err());
    assert!(qc_native("gene\tS1\ng1\tx\n").is_err());
    let three = json(
        &qc_native("gene\tA\tB\tC\ng1\t10\t200\t5\ng2\t500\t20\t90\ng3\t40\t40\t40\n").unwrap(),
    );
    assert_eq!(three["pca"]["scores"].as_array().unwrap().len(), 3);
}

const FASTA: &str = ">chr1\nACGTACGTACGTACGTACGTAAAACCCCGGGGTTTTACGT\n";
const GFF: &str = "##gff-version 3\nchr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=gene1;Name=g1\nchr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=gene1\nchr1\tsrc\texon\t1\t20\t.\t+\t.\tID=e1;Parent=tx1\n";

#[test]
fn transcriptome_builds_and_reports_bad_input() {
    let t = build_transcriptome_native(FASTA.as_bytes(), GFF).unwrap();
    assert_eq!(t.n_transcripts, 1);
    assert_eq!(t.n_genes, 1);
    assert!(String::from_utf8(t.fasta.clone())
        .unwrap()
        .starts_with(">tx1\nACGTACGTACGTACGTACGT\n"));
    assert_eq!(t.tx2gene, "tx1\tg1\n");
    assert!(t.gene_info_tsv.starts_with(
        "gene\tsymbol\tgene_id\told_locus_tag\ttranscripts\tproteins\tbiotype\tclass\tcontig\tstrand\tdbxref\n"
    ));
    // FASTA names that do not match the GFF3.
    let e = build_transcriptome_native(">other\nACGT\n".as_bytes(), GFF).unwrap_err();
    assert!(e.contains("chr1") && e.contains("names match"), "{e}");
    // GFF3 with no transcripts.
    let e = build_transcriptome_native(FASTA.as_bytes(), "##gff-version 3\n").unwrap_err();
    assert!(!e.is_empty());
    // Empty inputs.
    assert!(build_transcriptome_native(b"", GFF)
        .unwrap_err()
        .contains("empty"));
    assert!(build_transcriptome_native(FASTA.as_bytes(), "")
        .unwrap_err()
        .contains("empty"));
}

#[test]
fn kallisto_quant_args_refuse_other_thread_counts() {
    for t in ["2", "0", "4294967297"] {
        let req = format!(
            r#"{{"kind":"quant","sample":{{"name":"S","fastq1":["a_1.fq"],"fastq2":["a_2.fq"]}},"index":"i","out":"o","threads":{t}}}"#
        );
        let e = kallisto_args_native(&req).unwrap_err();
        assert!(e.contains("threads"), "{e}");
    }
}

#[test]
fn qc_nan_correlation_becomes_null() {
    // Sample C has the same count in every gene, so it has no variance.
    let q =
        json(&qc_native("gene\tA\tB\tC\ng1\t10\t200\t5\ng2\t500\t20\t5\ng3\t40\t40\t5\n").unwrap());
    let c = &q["correlation"];
    for k in 0..3 {
        assert!(c[2][k].is_null(), "row C, col {k}");
        assert!(c[k][2].is_null(), "row {k}, col C");
    }
    assert!(c[0][1].is_number());
    assert!(c[0][0].is_number());
}

#[test]
fn group_fastq_names_returns_rows_and_skips() {
    let v = json(&group_fastq_names_native(
        r#"{"names": ["a_R2.fq.gz", "a_R1.fq.gz", "notes.txt"], "numericMates": false}"#,
    ));
    assert_eq!(v["ok"], true);
    assert_eq!(v["samples"][0]["name"], "a");
    assert_eq!(v["samples"][0]["fastq1"][0], "a_R1.fq.gz");
    assert_eq!(v["samples"][0]["fastq2"][0], "a_R2.fq.gz");
    assert_eq!(v["samples"][0]["rule"], "read-tag");
    assert_eq!(v["skipped"][0]["path"], "notes.txt");
    let bad = json(&group_fastq_names_native("{}"));
    assert_eq!(bad["ok"], false);
    assert!(bad["error"].as_str().unwrap().contains("names"));
}
