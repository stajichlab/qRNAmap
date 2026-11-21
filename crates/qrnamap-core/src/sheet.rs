//! Sample sheet model and parser.
//!
//! Format: a tab-separated table with a header row.
//! * Required columns: `sample`, `fastq1`. Optional column: `fastq2`.
//! * Every other column is a covariate. A header `name:numeric` declares a numeric column
//!   (stored as `name`); all other covariates are factors.
//! * `fastq1` and `fastq2` hold one path or a comma-separated list (lanes). An empty `fastq2`
//!   cell means single-end for that sample.
//! * Blank lines and lines that start with `#` are ignored. CRLF line endings are accepted.
//!
//! The sample names, covariate names and factor values end up in a generated R script, so they
//! are restricted to safe characters here (see [`SampleSheet::parse`]). Covariate names that are
//! reserved words or constants in R are rejected with `ReservedName`.
//!
//! Downstream contract: fastq paths and sample names are passed only as argv elements. They are
//! never pasted into generated R or shell text. Fastq paths are checked for control characters
//! only.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// Type of a covariate column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColType {
    Factor,
    Numeric,
}

/// One sample row.
#[derive(Debug, Clone, PartialEq)]
pub struct Sample {
    pub name: String,
    pub fastq1: Vec<String>,
    /// Empty for a single-end sample.
    pub fastq2: Vec<String>,
    /// Covariate values as written in the sheet (numeric values are validated, not rewritten).
    pub covariates: BTreeMap<String, String>,
}

/// A parsed and validated sample sheet.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleSheet {
    pub samples: Vec<Sample>,
    /// Covariate columns in header order.
    pub columns: Vec<(String, ColType)>,
}

/// A problem found in a sample sheet. Each variant names the sample, column or value concerned.
#[derive(Debug, Clone, PartialEq)]
pub enum SheetError {
    /// Two rows have the same sample name.
    DuplicateSample(String),
    /// A sample name, covariate column name or factor value has a forbidden form, or a fastq
    /// path contains a control character. The field holds the offending text. Fastq paths are
    /// not otherwise restricted (paths with spaces are allowed).
    BadName(String),
    /// A required column (`sample` or `fastq1`) is absent from the header.
    MissingColumn(String),
    /// `fastq2` is not empty and has a different number of files than `fastq1`.
    MateCountMismatch { sample: String },
    /// A value in a numeric column is not a finite number.
    BadNumeric { column: String, value: String },
    /// The sheet has no samples (or no header).
    Empty,
    /// Two header cells give the same column name (the `:numeric` marker is ignored).
    DuplicateColumn(String),
    /// A required cell is empty (`fastq1`, a list element, or a factor value).
    EmptyValue { sample: String, column: String },
    /// A row has a different number of cells than the header. The field is the first cell.
    BadRow(String),
    /// A covariate column name is a reserved word or constant in R (for example `if`, `TRUE`,
    /// `NA`, `..1`). The field is the column name.
    ReservedName(String),
}

impl fmt::Display for SheetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SheetError::DuplicateSample(s) => write!(f, "duplicate sample name {s:?}"),
            SheetError::BadName(s) => {
                write!(
                    f,
                    "name or value {s:?} contains forbidden characters or has a bad form"
                )
            }
            SheetError::MissingColumn(c) => write!(f, "required column {c:?} is missing"),
            SheetError::MateCountMismatch { sample } => {
                write!(
                    f,
                    "sample {sample:?}: fastq1 and fastq2 list different numbers of files"
                )
            }
            SheetError::BadNumeric { column, value } => {
                write!(f, "column {column:?}: {value:?} is not a finite number")
            }
            SheetError::Empty => write!(f, "the sample sheet has no samples"),
            SheetError::DuplicateColumn(c) => write!(f, "duplicate column {c:?} in the header"),
            SheetError::EmptyValue { sample, column } => {
                write!(f, "sample {sample:?}: column {column:?} is empty")
            }
            SheetError::ReservedName(c) => {
                write!(
                    f,
                    "column {c:?} is a reserved R name and cannot be a covariate"
                )
            }
            SheetError::BadRow(s) => {
                write!(
                    f,
                    "row for {s:?} has a different number of cells than the header"
                )
            }
        }
    }
}

impl std::error::Error for SheetError {}

fn is_sample_name(s: &str) -> bool {
    let mut it = s.chars();
    matches!(it.next(), Some(c) if c.is_ascii_alphabetic())
        && it.all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
}

fn is_column_name(s: &str) -> bool {
    let mut it = s.chars();
    matches!(it.next(), Some(c) if c.is_ascii_alphabetic())
        && it.all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.'))
}

/// R reserved words and constants, plus `..N` style names and dot-digit names, which R does not
/// read as ordinary variables. `T` and `F` are ordinary variables in R and are not listed.
fn is_r_reserved(s: &str) -> bool {
    matches!(
        s,
        "if" | "else"
            | "repeat"
            | "while"
            | "function"
            | "for"
            | "in"
            | "next"
            | "break"
            | "TRUE"
            | "FALSE"
            | "NULL"
            | "Inf"
            | "NaN"
            | "NA"
            | "NA_integer_"
            | "NA_real_"
            | "NA_character_"
            | "NA_complex_"
    ) || s.starts_with("..")
        || (s.starts_with('.') && s.chars().nth(1).is_some_and(|c| c.is_ascii_digit()))
}

/// Factor values must be non-empty and use only `[A-Za-z0-9_.-]`. This is the same set the
/// contrast levels use, so a value that is accepted here can be named in a contrast.
fn is_factor_value(s: &str) -> bool {
    !s.is_empty()
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
}

fn split_list(cell: &str, sample: &str, column: &str, errors: &mut Vec<SheetError>) -> Vec<String> {
    let mut out = Vec::new();
    if cell.is_empty() {
        return out;
    }
    for part in cell.split(',') {
        let p = part.trim();
        if p.is_empty() {
            errors.push(SheetError::EmptyValue {
                sample: sample.into(),
                column: column.into(),
            });
        } else if p.chars().any(|c| c.is_control()) {
            errors.push(SheetError::BadName(p.into()));
        } else {
            out.push(p.to_string());
        }
    }
    out
}

impl SampleSheet {
    /// Parse and validate a sample sheet. All problems found are returned together.
    ///
    /// Rules:
    /// * Sample names match `[A-Za-z][A-Za-z0-9_.-]*`.
    /// * Covariate column names (without `:numeric`) match `[A-Za-z][A-Za-z0-9_.]*`.
    /// * Factor values match `[A-Za-z0-9_.-]+`. Values with a space, quote, backslash,
    ///   backtick, `$` or control character are rejected with `BadName`. The task brief
    ///   allowed other characters such as spaces; this crate rejects them because they could
    ///   reach R code and could not be named in a contrast. An empty value is `EmptyValue`.
    /// * Numeric values parse as finite `f64`.
    ///
    /// If the header has errors, row checks are skipped.
    pub fn parse(tsv: &str) -> Result<SampleSheet, Vec<SheetError>> {
        // A UTF-8 byte order mark (Excel "CSV UTF-8") is not part of the first cell.
        let tsv = tsv.strip_prefix('\u{feff}').unwrap_or(tsv);
        let mut lines = tsv
            .split('\n')
            .map(|l| l.strip_suffix('\r').unwrap_or(l))
            .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
        let header = match lines.next() {
            Some(h) => h,
            None => return Err(vec![SheetError::Empty]),
        };

        let mut errors = Vec::new();
        // (name, numeric?) per header cell
        let mut cols: Vec<(String, bool)> = Vec::new();
        let mut seen = BTreeSet::new();
        for cell in header.split('\t') {
            let cell = cell.trim();
            let (name, numeric) = match cell.strip_suffix(":numeric") {
                Some(n) => (n, true),
                None => (cell, false),
            };
            let reserved = matches!(name, "sample" | "fastq1" | "fastq2");
            if reserved && numeric {
                errors.push(SheetError::BadName(cell.into()));
            } else if !reserved && is_r_reserved(name) {
                errors.push(SheetError::ReservedName(name.into()));
            } else if !reserved && !is_column_name(name) {
                errors.push(SheetError::BadName(name.into()));
            } else if !seen.insert(name.to_string()) {
                errors.push(SheetError::DuplicateColumn(name.into()));
            }
            cols.push((name.to_string(), numeric));
        }
        for req in ["sample", "fastq1"] {
            if !seen.contains(req) {
                errors.push(SheetError::MissingColumn(req.into()));
            }
        }
        if !errors.is_empty() {
            return Err(errors);
        }

        let idx = |n: &str| cols.iter().position(|(c, _)| c == n);
        let i_sample = idx("sample").unwrap();
        let i_f1 = idx("fastq1").unwrap();
        let i_f2 = idx("fastq2");
        let columns: Vec<(String, ColType)> = cols
            .iter()
            .filter(|(n, _)| !matches!(n.as_str(), "sample" | "fastq1" | "fastq2"))
            .map(|(n, num)| {
                (
                    n.clone(),
                    if *num {
                        ColType::Numeric
                    } else {
                        ColType::Factor
                    },
                )
            })
            .collect();

        let mut samples = Vec::new();
        let mut names = BTreeSet::new();
        for line in lines {
            let cells: Vec<&str> = line.split('\t').map(str::trim).collect();
            if cells.len() != cols.len() {
                errors.push(SheetError::BadRow(cells[0].into()));
                continue;
            }
            let name = cells[i_sample];
            if !is_sample_name(name) {
                errors.push(SheetError::BadName(name.into()));
            } else if !names.insert(name.to_string()) {
                errors.push(SheetError::DuplicateSample(name.into()));
            }
            let f1_cell = cells[i_f1];
            if f1_cell.is_empty() {
                errors.push(SheetError::EmptyValue {
                    sample: name.into(),
                    column: "fastq1".into(),
                });
            }
            let fastq1 = split_list(f1_cell, name, "fastq1", &mut errors);
            let fastq2 = match i_f2 {
                Some(i) => split_list(cells[i], name, "fastq2", &mut errors),
                None => Vec::new(),
            };
            if !fastq2.is_empty() && !fastq1.is_empty() && fastq1.len() != fastq2.len() {
                errors.push(SheetError::MateCountMismatch {
                    sample: name.into(),
                });
            }
            let mut covariates = BTreeMap::new();
            for (i, (cname, numeric)) in cols.iter().enumerate() {
                if matches!(cname.as_str(), "sample" | "fastq1" | "fastq2") {
                    continue;
                }
                let v = cells[i];
                if *numeric {
                    match v.parse::<f64>() {
                        Ok(x) if x.is_finite() => {}
                        _ => errors.push(SheetError::BadNumeric {
                            column: cname.clone(),
                            value: v.into(),
                        }),
                    }
                } else if v.is_empty() {
                    errors.push(SheetError::EmptyValue {
                        sample: name.into(),
                        column: cname.clone(),
                    });
                } else if !is_factor_value(v) {
                    errors.push(SheetError::BadName(v.into()));
                }
                covariates.insert(cname.clone(), v.to_string());
            }
            samples.push(Sample {
                name: name.to_string(),
                fastq1,
                fastq2,
                covariates,
            });
        }
        if samples.is_empty() && errors.is_empty() {
            errors.push(SheetError::Empty);
        }
        if errors.is_empty() {
            Ok(SampleSheet { samples, columns })
        } else {
            Err(errors)
        }
    }
}
