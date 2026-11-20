//! DESeq2 step: design checks and the generated R script. No I/O.
//!
//! The script holds no user text. Paths (counts file, coldata file, output prefix) arrive through
//! `commandArgs(trailingOnly = TRUE)`. The only values written into the script are validated
//! column names, the contrast level and reference, and the design formula from
//! [`Design::to_r`]. Names and levels are written as R string literals with `\` and `"`
//! escaped, as a second guard after the sheet parser's character checks.

use crate::counts::CountMatrix;
use crate::design::{Contrast, Design};
use crate::sheet::{ColType, SampleSheet};
use std::collections::BTreeSet;
use std::fmt;

/// Method of `lfcShrink`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shrink {
    Normal,
    Apeglm,
    Ashr,
}

impl Shrink {
    pub fn name(self) -> &'static str {
        match self {
            Shrink::Normal => "normal",
            Shrink::Apeglm => "apeglm",
            Shrink::Ashr => "ashr",
        }
    }
}

/// A DE input problem. Every variant names the cause.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DeError {
    /// A tested group has fewer than 2 samples.
    TooFewReplicates { group: String, n: usize },
    /// All counts of this sample are zero.
    AllZeroSample(String),
    /// R reported that the design is not full rank. The field is the tail of R's stderr.
    NotFullRank(String),
    /// The contrast level or reference occurs in no sample.
    LevelNotInData { factor: String, level: String },
    /// The sample sets of the matrix and the sheet differ, or the matrix is malformed.
    SampleMismatch(String),
    /// The contrast level and reference are equal.
    SameLevelAndReference(String),
    /// The contrast factor is not a factor column of the sheet.
    ContrastNotFactor(String),
    /// The contrast factor is not a term of the design.
    ContrastNotInDesign(String),
}

impl fmt::Display for DeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DeError::TooFewReplicates { group, n } => write!(
                f,
                "group '{group}' has {n} sample(s); DESeq2 needs at least 2 per tested group"
            ),
            DeError::AllZeroSample(s) => write!(f, "sample '{s}' has all counts equal to zero"),
            DeError::NotFullRank(m) => write!(f, "the design is not full rank. R said:\n{m}"),
            DeError::LevelNotInData { factor, level } => {
                write!(f, "level '{level}' of factor '{factor}' occurs in no sample")
            }
            DeError::SampleMismatch(m) => write!(f, "sample mismatch between the count matrix and the sheet: {m}"),
            DeError::SameLevelAndReference(l) => {
                write!(f, "contrast level and reference are both '{l}'")
            }
            DeError::ContrastNotFactor(c) => write!(f, "contrast column '{c}' is not a factor column of the sheet"),
            DeError::ContrastNotInDesign(c) => write!(f, "contrast factor '{c}' is not a term of the design"),
        }
    }
}

impl std::error::Error for DeError {}

/// Check the inputs of a DE run before R starts.
///
/// Order: contrast factor is a factor column; level differs from reference; the matrix sample
/// set equals the sheet sample set; level and reference occur in the data; each has at least 2
/// samples; no sample has all counts zero.
pub fn check_design(
    sheet: &SampleSheet,
    contrast: &Contrast,
    matrix: &CountMatrix,
) -> Result<(), DeError> {
    match sheet.columns.iter().find(|(c, _)| *c == contrast.factor) {
        Some((_, ColType::Factor)) => {}
        _ => return Err(DeError::ContrastNotFactor(contrast.factor.clone())),
    }
    if contrast.level == contrast.reference {
        return Err(DeError::SameLevelAndReference(contrast.level.clone()));
    }
    // Matrix shape and sample set.
    if matrix.counts.len() != matrix.genes.len() || matrix.counts.iter().any(|r| r.len() != matrix.samples.len()) {
        return Err(DeError::SampleMismatch("the count matrix rows do not match its header".into()));
    }
    let in_matrix: BTreeSet<&str> = matrix.samples.iter().map(String::as_str).collect();
    if in_matrix.len() != matrix.samples.len() {
        return Err(DeError::SampleMismatch("a sample name occurs twice in the count matrix".into()));
    }
    let in_sheet: BTreeSet<&str> = sheet.samples.iter().map(|s| s.name.as_str()).collect();
    if let Some(s) = in_sheet.difference(&in_matrix).next() {
        return Err(DeError::SampleMismatch(format!("sample '{s}' is in the sheet but not in the matrix")));
    }
    if let Some(s) = in_matrix.difference(&in_sheet).next() {
        return Err(DeError::SampleMismatch(format!("sample '{s}' is in the matrix but not in the sheet")));
    }
    // Replicates per tested group.
    let n_in = |level: &str| {
        sheet
            .samples
            .iter()
            .filter(|s| s.covariates.get(&contrast.factor).map(String::as_str) == Some(level))
            .count()
    };
    for group in [&contrast.level, &contrast.reference] {
        let n = n_in(group);
        if n == 0 {
            return Err(DeError::LevelNotInData { factor: contrast.factor.clone(), level: group.clone() });
        }
    }
    for group in [&contrast.level, &contrast.reference] {
        let n = n_in(group);
        if n < 2 {
            return Err(DeError::TooFewReplicates { group: group.clone(), n });
        }
    }
    // All-zero samples.
    for (j, name) in matrix.samples.iter().enumerate() {
        if matrix.counts.iter().all(|row| row[j] == 0) {
            return Err(DeError::AllZeroSample(name.clone()));
        }
    }
    Ok(())
}

/// Columns the DE step needs from the sheet: the design terms, once each, in formula order.
fn used_columns(design: &Design) -> Vec<&str> {
    let mut seen = Vec::new();
    for t in design.terms() {
        if !seen.contains(&t.as_str()) {
            seen.push(t.as_str());
        }
    }
    seen
}

/// Tab-separated `coldata`: `sample` plus the design columns. Fastq paths and other sheet
/// columns are never written. One row per sample, in sheet order.
pub fn coldata_tsv(sheet: &SampleSheet, design: &Design) -> String {
    let cols = used_columns(design);
    let mut s = String::from("sample");
    for c in &cols {
        s.push('\t');
        s.push_str(c);
    }
    s.push('\n');
    for smp in &sheet.samples {
        s.push_str(&smp.name);
        for c in &cols {
            s.push('\t');
            s.push_str(smp.covariates.get(*c).map(String::as_str).unwrap_or(""));
        }
        s.push('\n');
    }
    s
}

/// An R string literal: `"` and `\` escaped; control characters written as escapes.
pub fn r_string(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '\\' => o.push_str("\\\\"),
            '"' => o.push_str("\\\""),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 || c as u32 == 0x7f => o.push_str(&format!("\\x{:02x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

/// The R script. Arguments: `counts_file coldata_file out_prefix`.
///
/// `columns` lists the coldata columns (name and type) as written by [`coldata_tsv`]. Factor
/// columns are converted with `factor()`. The contrast factor gets the reference level first.
/// Numeric columns stay numeric. `shrink` adds an `lfcShrink` step after the unshrunk results
/// are written.
pub fn deseq2_script(
    design: &Design,
    contrast: &Contrast,
    columns: &[(String, ColType)],
    shrink: Option<Shrink>,
) -> String {
    let q = r_string;
    let mut s = String::new();
    let mut line = |t: &str| {
        s.push_str(t);
        s.push('\n');
    };
    line("# Generated by qrnamap. Paths arrive as arguments; this script holds no paths.");
    line("args <- commandArgs(trailingOnly = TRUE)");
    line("if (length(args) != 3) stop(\"usage: script.R counts coldata prefix\")");
    line("suppressPackageStartupMessages(library(DESeq2))");
    line("counts <- read.delim(args[1], row.names = 1, check.names = FALSE, quote = \"\",");
    line("  comment.char = \"\", na.strings = character(0), colClasses = c(\"gene\" = \"character\"))");
    line("counts <- as.matrix(counts)");
    line("storage.mode(counts) <- \"integer\"");
    let cc: Vec<String> = std::iter::once(format!("{} = \"character\"", q("sample")))
        .chain(columns.iter().map(|(n, t)| {
            format!("{} = \"{}\"", q(n), if *t == ColType::Numeric { "numeric" } else { "character" })
        }))
        .collect();
    line("coldata <- read.delim(args[2], row.names = 1, check.names = FALSE, quote = \"\",");
    line(&format!("  comment.char = \"\", na.strings = character(0), colClasses = c({}))", cc.join(", ")));
    for (n, t) in columns {
        if *t != ColType::Factor {
            continue;
        }
        if *n == contrast.factor {
            let r = q(&contrast.reference);
            line(&format!(
                "coldata[[{n}]] <- factor(coldata[[{n}]], levels = c({r}, setdiff(sort(unique(coldata[[{n}]])), {r})))",
                n = q(n)
            ));
        } else {
            line(&format!("coldata[[{n}]] <- factor(coldata[[{n}]])", n = q(n)));
        }
    }
    line("if (!setequal(colnames(counts), rownames(coldata))) {");
    line("  stop(\"sample names differ between the count matrix and the coldata\")");
    line("}");
    line("coldata <- coldata[colnames(counts), , drop = FALSE]");
    line(&format!("dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = {})", design.to_r()));
    line("dds <- DESeq(dds, test = \"Wald\", quiet = TRUE)");
    line(&format!(
        "contrast_spec <- c({}, {}, {})",
        q(&contrast.factor),
        q(&contrast.level),
        q(&contrast.reference)
    ));
    line("res <- results(dds, contrast = contrast_spec)");
    line("write_table <- function(x, path) {");
    line("  out <- data.frame(gene = rownames(x), as.data.frame(x), check.names = FALSE)");
    line("  write.csv(out, path, row.names = FALSE)");
    line("}");
    line("write_table(res, paste0(args[3], \".results.csv\"))");
    line("write_table(counts(dds, normalized = TRUE), paste0(args[3], \".normalized.csv\"))");
    line("vsd <- varianceStabilizingTransformation(dds, blind = FALSE)");
    line("write_table(assay(vsd), paste0(args[3], \".vst.csv\"))");
    line("write_versions <- function(shrink_label) {");
    line("  writeLines(c(");
    line("    paste(\"DESeq2\", as.character(packageVersion(\"DESeq2\"))),");
    line("    paste(\"R\", R.version.string),");
    line(&format!("    paste(\"design:\", {}),", q(&design.to_r())));
    line(&format!(
        "    paste(\"contrast:\", {}),",
        q(&format!("{},{},{}", contrast.factor, contrast.level, contrast.reference))
    ));
    line("    paste(\"test: Wald\"),");
    line("    paste(\"shrink:\", shrink_label)");
    line("  ), paste0(args[3], \".versions.txt\"))");
    line("}");
    line("write_versions(\"none\")");
    line("cat(\"DESeq2\", as.character(packageVersion(\"DESeq2\")), \"\\n\")");
    if let Some(m) = shrink {
        let name = m.name();
        if m != Shrink::Normal {
            line(&format!("if (!requireNamespace(\"{name}\", quietly = TRUE)) {{"));
            line(&format!(
                "  stop(\"MISSING_PACKAGE:{name} is not installed; the unshrunk results were written\")"
            ));
            line("}");
        }
        match m {
            Shrink::Apeglm => {
                line(&format!(
                    "coef_name <- paste0({}, \"_\", {}, \"_vs_\", {})",
                    q(&contrast.factor),
                    q(&contrast.level),
                    q(&contrast.reference)
                ));
                line("if (!(coef_name %in% resultsNames(dds))) stop(\"apeglm needs the coefficient named in resultsNames(dds)\")");
                line("shrunk <- lfcShrink(dds, coef = coef_name, type = \"apeglm\", quiet = TRUE)");
            }
            Shrink::Normal => {
                line("shrunk <- lfcShrink(dds, contrast = contrast_spec, type = \"normal\", quiet = TRUE)");
            }
            Shrink::Ashr => {
                line("shrunk <- lfcShrink(dds, contrast = contrast_spec, res = res, type = \"ashr\", quiet = TRUE)");
            }
        }
        line("write_table(shrunk, paste0(args[3], \".results.shrunk.csv\"))");
        line(&format!("write_versions(\"{name}\")"));
    }
    s
}
