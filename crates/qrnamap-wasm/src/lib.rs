//! WebAssembly wrapper for the browser page. It exposes `qrnamap-core`, `qseq-io` and
//! `genome-catalog` to JavaScript. The logic lives in those crates. This crate converts
//! between JavaScript values and Rust values (JSON strings, byte arrays) and does nothing else.
//!
//! Every exported function is a thin layer over a plain Rust function in this file, so native
//! tests can call it. Errors are strings. JSON numbers that are not finite (NaN) become `null`.

use std::collections::BTreeMap;
use std::io::Cursor;

use qrnamap_core::{
    build_transcripts, count_matrix, gene_counts_from_abundance, gene_info_tsv, gene_table,
    gene_table_tsv, genes_from_tx2gene, kallisto_index, kallisto_quant, log_cpm, parse_tx2gene,
    pca, reserved_output_name, sample_correlation, tx2gene_tsv, write_fasta, ColType, CountMatrix,
    QcError, RunConfig, Sample, SampleSheet, Strandedness,
};
use qseq_io::{parse_gff3, read_fasta};
use serde_json::{json, Value};
use wasm_bindgen::prelude::*;

/// The number of top-variance genes the page uses for the PCA.
pub const PCA_TOP_GENES: usize = 500;

// ---------------------------------------------------------------------------------------------
// Catalog

/// Genomes that have a GFF3 (the ones RNA-seq quantification can use), as a JSON array.
pub fn catalog_json() -> String {
    let cat = genome_catalog::builtin();
    let file =
        |f: &genome_catalog::FileRef| json!({"name": f.name, "url": f.url, "bytes": f.bytes});
    let genomes: Vec<Value> = cat
        .usable_for_rna()
        .into_iter()
        .map(|g| {
            json!({
                "accession": g.accession,
                "source": g.source,
                "label": g.label,
                "organism": g.organism,
                "strain": g.strain,
                "level": g.level,
                "assembly": g.assembly,
                "genome_size": g.genome_size,
                "genes": g.genes,
                "fasta": file(&g.fasta),
                "gff": g.gff.as_ref().map(file),
                "size_warning": genome_catalog::size_warning(g),
            })
        })
        .collect();
    Value::Array(genomes).to_string()
}

// ---------------------------------------------------------------------------------------------
// Transcriptome

/// Result of `build_transcriptome`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptomeOut {
    pub fasta: Vec<u8>,
    pub tx2gene: String,
    pub genes_tsv: String,
    pub gene_info_tsv: String,
    pub n_transcripts: usize,
    pub n_genes: usize,
}

fn check_label(l: &str) -> Result<(), String> {
    if l.contains(['\t', '\n', '\r']) {
        return Err(format!(
            "gene or transcript name {l:?} has a tab, newline or carriage return and cannot go in a table"
        ));
    }
    Ok(())
}

/// Transcript FASTA, `tx2gene` and gene table from a genome FASTA (plain bytes) and a GFF3
/// (plain text). The same steps as `qrnamap transcripts`.
pub fn build_transcriptome_native(fasta: &[u8], gff3: &str) -> Result<TranscriptomeOut, String> {
    if fasta.iter().all(|b| b.is_ascii_whitespace()) {
        return Err("the genome FASTA is empty".into());
    }
    if gff3.trim().is_empty() {
        return Err("the GFF3 is empty".into());
    }
    let genome =
        read_fasta(Cursor::new(fasta)).map_err(|e| format!("cannot read the genome FASTA: {e}"))?;
    if genome.is_empty() {
        return Err("the genome FASTA has no sequences".into());
    }
    let model = parse_gff3(Cursor::new(gff3.as_bytes()))
        .map_err(|e| format!("cannot read the GFF3: {e}"))?;
    for g in &model.genes {
        check_label(&g.label)?;
        for t in &g.transcripts {
            check_label(&t.id)?;
        }
    }
    let records = build_transcripts(&genome, &model).map_err(|e| {
        format!("{e}. Check that the FASTA sequence names match the first column of the GFF3.")
    })?;
    if records.is_empty() {
        return Err(
            "the GFF3 gives no usable transcripts (no gene has an exon, CDS or UTR interval)"
                .into(),
        );
    }
    let mut fasta_out = Vec::new();
    write_fasta(&records, &mut fasta_out).map_err(|e| e.to_string())?;
    let tx2gene = tx2gene_tsv(&records);
    let genes_tsv = gene_table_tsv(&gene_table(&model, &[]));
    let gene_info_tsv = gene_info_tsv(&model, &[]);
    let n_genes = genes_from_tx2gene(&tx2gene)
        .map_err(|e| e.to_string())?
        .len();
    Ok(TranscriptomeOut {
        n_transcripts: records.len(),
        n_genes,
        fasta: fasta_out,
        tx2gene,
        genes_tsv,
        gene_info_tsv,
    })
}

/// The transcriptome made from a genome and a gene model.
#[wasm_bindgen]
pub struct Transcriptome {
    inner: TranscriptomeOut,
}

#[wasm_bindgen]
impl Transcriptome {
    /// Transcript FASTA as bytes (same bytes as `qrnamap transcripts --out x.fa`).
    #[wasm_bindgen(getter)]
    pub fn fasta(&self) -> Vec<u8> {
        self.inner.fasta.clone()
    }
    /// `transcript<TAB>gene` rows (same bytes as `qrnamap transcripts --tx2gene x.tsv`).
    #[wasm_bindgen(getter)]
    pub fn tx2gene(&self) -> String {
        self.inner.tx2gene.clone()
    }
    /// The gene table (same bytes as `qrnamap transcripts --genes-table x.tsv`).
    #[wasm_bindgen(getter)]
    pub fn genes_tsv(&self) -> String {
        self.inner.genes_tsv.clone()
    }
    /// The gene alias table (same bytes as `qrnamap transcripts` writes next to `genes.tsv`).
    #[wasm_bindgen(getter)]
    pub fn gene_info_tsv(&self) -> String {
        self.inner.gene_info_tsv.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn n_transcripts(&self) -> usize {
        self.inner.n_transcripts
    }
    #[wasm_bindgen(getter)]
    pub fn n_genes(&self) -> usize {
        self.inner.n_genes
    }
}

// ---------------------------------------------------------------------------------------------
// Gene counts and the matrix

/// Gene counts of one kallisto `abundance.tsv`. JSON: `{counts: {gene: n}, total, n_genes}`.
pub fn gene_counts_native(abundance_tsv: &str, tx2gene_text: &str) -> Result<String, String> {
    let tx2gene = parse_tx2gene(tx2gene_text).map_err(|e| e.to_string())?;
    let counts = gene_counts_from_abundance(abundance_tsv, &tx2gene).map_err(|e| e.to_string())?;
    let total: u64 = counts.values().sum();
    Ok(json!({"counts": counts, "total": total, "n_genes": counts.len()}).to_string())
}

/// The gene-by-sample matrix. Input JSON: `{"tx2gene": "<text>", "samples": [{"name": "S1",
/// "counts": {"gene": n}}]}`. The rows are the genes of `tx2gene` in first-appearance order,
/// as in the CLI. Output JSON: `{tsv, warnings, n_genes, samples}`.
pub fn assemble_matrix_native(samples_json: &str) -> Result<String, String> {
    let v: Value = serde_json::from_str(samples_json).map_err(|e| format!("bad JSON: {e}"))?;
    let tx2gene = v
        .get("tx2gene")
        .and_then(Value::as_str)
        .ok_or("missing 'tx2gene' text")?;
    let genes = genes_from_tx2gene(tx2gene).map_err(|e| e.to_string())?;
    let samples = v
        .get("samples")
        .and_then(Value::as_array)
        .ok_or("missing 'samples' array")?;
    let mut per_sample: Vec<(String, BTreeMap<String, u64>)> = Vec::with_capacity(samples.len());
    for s in samples {
        let name = s
            .get("name")
            .and_then(Value::as_str)
            .ok_or("a sample has no 'name'")?;
        if name.is_empty() || name.contains(['\t', '\n', '\r']) {
            return Err(format!(
                "sample name {name:?} is empty or has a tab or newline"
            ));
        }
        if per_sample.iter().any(|(n, _)| n == name) {
            return Err(format!("duplicate sample name '{name}'"));
        }
        let obj = s
            .get("counts")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("sample '{name}' has no 'counts' object"))?;
        let mut counts = BTreeMap::new();
        for (g, n) in obj {
            let n = n.as_u64().ok_or_else(|| {
                format!("sample '{name}': count of gene '{g}' is not a non-negative integer")
            })?;
            counts.insert(g.clone(), n);
        }
        per_sample.push((name.to_string(), counts));
    }
    let (matrix, warnings) = count_matrix(&per_sample, Some(&genes)).map_err(|e| e.to_string())?;
    Ok(json!({
        "tsv": matrix.to_tsv(),
        "warnings": warnings,
        "n_genes": matrix.genes.len(),
        "samples": matrix.samples,
    })
    .to_string())
}

// ---------------------------------------------------------------------------------------------
// QC

/// Read a matrix TSV (`gene<TAB>sample...`, then one integer row per gene).
pub fn parse_matrix_tsv(text: &str) -> Result<CountMatrix, String> {
    let mut lines = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.trim_end_matches('\r')))
        .filter(|(_, l)| !l.trim().is_empty());
    let (_, header) = lines.next().ok_or("the count matrix is empty")?;
    let cols: Vec<&str> = header.split('\t').collect();
    if cols.len() < 2 || cols[0] != "gene" {
        return Err(
            "the count matrix header must start with 'gene' and name at least one sample".into(),
        );
    }
    let samples: Vec<String> = cols[1..].iter().map(|s| s.to_string()).collect();
    let mut genes = Vec::new();
    let mut counts = Vec::new();
    for (line, row) in lines {
        let f: Vec<&str> = row.split('\t').collect();
        if f.len() != cols.len() {
            return Err(format!(
                "count matrix line {line}: expected {} columns",
                cols.len()
            ));
        }
        let nums: Result<Vec<u64>, _> = f[1..].iter().map(|c| c.trim().parse::<u64>()).collect();
        counts.push(nums.map_err(|_| {
            format!("count matrix line {line}: a count is not a non-negative integer")
        })?);
        genes.push(f[0].to_string());
    }
    Ok(CountMatrix {
        genes,
        samples,
        counts,
    })
}

/// QC numbers for a matrix TSV. JSON: `{samples, n_genes, correlation, correlation_error, pca,
/// pca_error, pca_method}`. Correlation is a samples-by-samples array or null. `pca` is
/// `{scores: [{sample, pc1, pc2}], variance_explained: [pc1, pc2], top_genes}` or null. An
/// error string says why a part is null (too few samples). A zero-library sample or a bad
/// matrix is an `Err`.
pub fn qc_native(matrix_tsv: &str) -> Result<String, String> {
    let m = parse_matrix_tsv(matrix_tsv)?;
    if let Err(QcError::ZeroLibrary(s)) = log_cpm(&m) {
        return Err(QcError::ZeroLibrary(s).to_string());
    }
    let (correlation, correlation_error) = match sample_correlation(&m) {
        Ok(c) => (json!(c), Value::Null),
        Err(e) => (Value::Null, json!(e.to_string())),
    };
    let (pca_v, pca_error) = match pca(&m, PCA_TOP_GENES) {
        Ok(p) => {
            let scores: Vec<Value> = m
                .samples
                .iter()
                .zip(&p.scores)
                .map(|(s, xy)| json!({"sample": s, "pc1": xy[0], "pc2": xy[1]}))
                .collect();
            (
                json!({
                    "scores": scores,
                    "variance_explained": p.variance_explained,
                    "top_genes": PCA_TOP_GENES.min(m.genes.len()),
                }),
                Value::Null,
            )
        }
        Err(e) => (Value::Null, json!(e.to_string())),
    };
    Ok(json!({
        "samples": m.samples,
        "n_genes": m.genes.len(),
        "correlation": correlation,
        "correlation_error": correlation_error,
        "pca": pca_v,
        "pca_error": pca_error,
        "pca_method": "PCA of log2(CPM + 1) on the highest-variance genes (not DESeq2 variance-stabilised)",
    })
    .to_string())
}

// ---------------------------------------------------------------------------------------------
// Commands

fn strandedness(s: &str) -> Result<Strandedness, String> {
    match s {
        "" | "unstranded" => Ok(Strandedness::Unstranded),
        "forward" => Ok(Strandedness::Forward),
        "reverse" => Ok(Strandedness::Reverse),
        other => Err(format!(
            "strandedness must be unstranded, forward or reverse, not '{other}'"
        )),
    }
}

fn str_list(v: Option<&Value>) -> Result<Vec<String>, String> {
    match v {
        None | Some(Value::Null) => Ok(vec![]),
        Some(Value::Array(a)) => a
            .iter()
            .map(|x| {
                x.as_str()
                    .map(String::from)
                    .ok_or_else(|| "a file name is not a string".to_string())
            })
            .collect(),
        Some(_) => Err("a file list must be an array of strings".into()),
    }
}

/// The kallisto argument list (after the program name) from `toolcall`, as a JSON array.
/// Request JSON, index: `{"kind":"index","fasta":"transcripts.fa","index":"transcripts.idx"}`.
/// Quant: `{"kind":"quant","sample":{"name","fastq1":[..],"fastq2":[..]},"index","out",
/// "threads":1,"strandedness":"unstranded","fragment_len":[mean,sd]}`.
pub fn kallisto_args_native(request_json: &str) -> Result<String, String> {
    let v: Value = serde_json::from_str(request_json).map_err(|e| format!("bad JSON: {e}"))?;
    let text = |k: &str| -> Result<&str, String> {
        v.get(k)
            .and_then(Value::as_str)
            .ok_or_else(|| format!("missing '{k}'"))
    };
    match text("kind")? {
        "index" => {
            let call = kallisto_index(text("fasta")?, text("index")?);
            Ok(json!(call.args).to_string())
        }
        "quant" => {
            let s = v.get("sample").ok_or("missing 'sample'")?;
            let sample = Sample {
                name: s
                    .get("name")
                    .and_then(Value::as_str)
                    .unwrap_or("sample")
                    .to_string(),
                fastq1: str_list(s.get("fastq1"))?,
                fastq2: str_list(s.get("fastq2"))?,
                covariates: BTreeMap::new(),
            };
            // The vendored kallisto WASM requires `-t 1` when a thread count is supplied (upstream
            // README), so the browser export accepts only 1 and refuses any other value.
            let threads = u32::try_from(v.get("threads").and_then(Value::as_u64).unwrap_or(1))
                .map_err(|_| "'threads' is too large".to_string())?;
            if threads != 1 {
                return Err(format!(
                    "threads must be 1 in the browser build (kallisto WASM requires -t 1), not {threads}"
                ));
            }
            let mut cfg = RunConfig {
                threads,
                strandedness: strandedness(
                    v.get("strandedness").and_then(Value::as_str).unwrap_or(""),
                )?,
                ..RunConfig::default()
            };
            if let Some(Value::Array(fl)) = v.get("fragment_len") {
                if let [Some(mean), Some(sd)] = [
                    fl.first().and_then(Value::as_f64),
                    fl.get(1).and_then(Value::as_f64),
                ] {
                    cfg.fragment_len = Some((mean, sd));
                }
            }
            let call = kallisto_quant(&sample, text("index")?, text("out")?, &cfg)
                .map_err(|e| e.to_string())?;
            Ok(json!(call.args).to_string())
        }
        other => Err(format!("kind must be index or quant, not '{other}'")),
    }
}

/// Quote one word for a POSIX shell: single quotes, unless the word is plain.
fn shell_word(s: &str) -> String {
    if !s.is_empty()
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_-./,:=+@%".contains(&b))
    {
        s.to_string()
    } else {
        format!("'{}'", s.replace('\'', r"'\''"))
    }
}

/// The `qrnamap de` command line for the sample sheet and matrix the user saved.
pub fn de_command_native(
    sheet_name: &str,
    counts_name: &str,
    design: &str,
    contrast: &str,
) -> String {
    format!(
        "qrnamap de --counts {} --samples {} --design {} --contrast {} --out de",
        shell_word(counts_name),
        shell_word(sheet_name),
        shell_word(design),
        shell_word(contrast)
    )
}

// ---------------------------------------------------------------------------------------------
// Sample sheet

/// Check a sample sheet with the core parser, and refuse the sample names that the command-line
/// tool refuses (`RESERVED_OUTPUT_NAMES`, with the same message). JSON: `{ok, errors: [text],
/// samples: [name], columns: [[name, "factor"|"numeric"]], rows: [{sample, fastq1: [path],
/// fastq2: [path], covariates: {column: value}}]}`. On an error, `samples`, `columns` and `rows`
/// are empty. This function does not fail.
pub fn validate_sheet_native(sheet_tsv: &str) -> String {
    let failed = |errors: Vec<String>| json!({"ok": false, "errors": errors, "samples": [], "columns": [], "rows": []});
    let sheet = match SampleSheet::parse(sheet_tsv) {
        Ok(sheet) => sheet,
        Err(errs) => return failed(errs.iter().map(|e| e.to_string()).collect()).to_string(),
    };
    // Same text as the command-line error PipelineError::ReservedSampleName.
    let reserved: Vec<String> = sheet
        .samples
        .iter()
        .filter_map(|s| reserved_output_name(&s.name))
        .map(|r| {
            format!(
                "sample '{r}': the name is reserved: qrnamap writes '{r}' in the output \
                 directory. Rename the sample in the sheet"
            )
        })
        .collect();
    if !reserved.is_empty() {
        return failed(reserved).to_string();
    }
    json!({
        "ok": true,
        "errors": [],
        "samples": sheet.samples.iter().map(|s| s.name.clone()).collect::<Vec<_>>(),
        "columns": sheet.columns.iter().map(|(n, t)| json!([n, match t {
            ColType::Factor => "factor",
            ColType::Numeric => "numeric",
        }])).collect::<Vec<_>>(),
        "rows": sheet.samples.iter().map(|s| json!({
            "sample": s.name,
            "fastq1": s.fastq1,
            "fastq2": s.fastq2,
            "covariates": s.covariates,
        })).collect::<Vec<_>>(),
    })
    .to_string()
}

/// Sample sheet rows from FASTQ file names (see `qrnamap_core::fastq_names`). The request is
/// `{"names": [path, ...], "numericMates": bool}` (numericMates defaults to true). Returns JSON
/// `{"ok", "error", "samples": [{name, fromName, fastq1, fastq2, rule, problems, notes}],
/// "skipped": [{path, reason}]}` and does not throw.
pub fn group_fastq_names_native(request_json: &str) -> String {
    use qrnamap_core::fastq_names::{group_fastq_names, GroupOptions};
    let fail = |error: String| {
        json!({"ok": false, "error": error, "samples": [], "skipped": []}).to_string()
    };
    let request: Value = match serde_json::from_str(request_json) {
        Ok(v) => v,
        Err(x) => return fail(format!("bad request: {x}")),
    };
    let Some(names) = request["names"].as_array() else {
        return fail("bad request: names must be an array of strings".into());
    };
    let names: Vec<&str> = names.iter().filter_map(Value::as_str).collect();
    let opts = GroupOptions {
        numeric_mates: request["numericMates"].as_bool().unwrap_or(true),
    };
    let g = group_fastq_names(&names, &opts);
    json!({
        "ok": true,
        "error": null,
        "samples": g.samples.iter().map(|r| json!({
            "name": r.name,
            "fromName": r.from_name,
            "fastq1": r.fastq1,
            "fastq2": r.fastq2,
            "rule": r.rule.as_str(),
            "problems": r.problems,
            "notes": r.notes,
        })).collect::<Vec<_>>(),
        "skipped": g.skipped.iter().map(|s| json!({"path": s.path, "reason": s.reason})).collect::<Vec<_>>(),
    })
    .to_string()
}

// ---------------------------------------------------------------------------------------------
// wasm-bindgen exports

/// Genomes for the menu, as a JSON array. See `catalog_json`.
#[wasm_bindgen]
pub fn catalog() -> String {
    catalog_json()
}

/// Build the transcriptome from the genome FASTA bytes (plain, not gzip) and the GFF3 text.
#[wasm_bindgen]
pub fn build_transcriptome(fasta_bytes: &[u8], gff3_text: &str) -> Result<Transcriptome, String> {
    build_transcriptome_native(fasta_bytes, gff3_text).map(|inner| Transcriptome { inner })
}

/// Gene counts of one kallisto `abundance.tsv` as JSON.
#[wasm_bindgen]
pub fn gene_counts(abundance_tsv: &str, tx2gene_tsv: &str) -> Result<String, String> {
    gene_counts_native(abundance_tsv, tx2gene_tsv)
}

/// The gene-by-sample matrix as JSON (`tsv`, `warnings`).
#[wasm_bindgen]
pub fn assemble_matrix(samples_json: &str) -> Result<String, String> {
    assemble_matrix_native(samples_json)
}

/// Correlation and PCA as JSON.
#[wasm_bindgen]
pub fn qc(matrix_tsv: &str) -> Result<String, String> {
    qc_native(matrix_tsv)
}

/// The kallisto argument list as JSON (see `kallisto_args_native`).
#[wasm_bindgen]
pub fn kallisto_args(request_json: &str) -> Result<String, String> {
    kallisto_args_native(request_json)
}

/// The `qrnamap de` command for the saved files.
#[wasm_bindgen]
pub fn de_command(sheet_name: &str, counts_name: &str, design: &str, contrast: &str) -> String {
    de_command_native(sheet_name, counts_name, design, contrast)
}

/// Validate a sample sheet. Returns JSON and does not throw.
#[wasm_bindgen]
pub fn validate_sheet(sheet_tsv: &str) -> String {
    validate_sheet_native(sheet_tsv)
}

/// Sample sheet rows from FASTQ file names. Returns JSON and does not throw. See
/// `group_fastq_names_native`.
#[wasm_bindgen]
pub fn group_fastq_names(request_json: &str) -> String {
    group_fastq_names_native(request_json)
}
