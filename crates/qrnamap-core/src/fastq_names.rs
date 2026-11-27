//! Sample sheet rows from FASTQ file names: mate pairing, lane grouping and safe sample names.
//!
//! Only the file names are read (never the file contents), and the rules are fixed, so the same
//! names always give the same rows. The rules, tried in this order on the name without its
//! extension (`.fastq`, `.fq`, each optionally `.gz` or `.zst`):
//!
//! 1. **Illumina** (bcl2fastq / BCL Convert): `<sample>_S<n>[_L<lane>]_R<1|2>_<chunk>`, for
//!    example `WT1_S1_L001_R1_001`. Index reads (`_I1_`, `_I2_`) and `Undetermined_*` files are
//!    skipped. The `_S<n>` number is dropped from the sample name unless two samples would then
//!    get the same name.
//! 2. **Read tag**: a name ending in `_R1`/`_R2`, `.R1`/`.R2`, `_r1`/`_r2` or `_read1`/`_read2`,
//!    optionally after a lane (`_L001_R1`).
//! 3. **Number** (optional, [`GroupOptions::numeric_mates`]): a name ending in `_1`/`_2` or
//!    `.1`/`.2` (SRA and ENA files, `SRR123_1`). Because `ctrl_1` and `ctrl_2` can also be two
//!    single-end replicates, a number only counts as a mate when both `_1` and `_2` of the same
//!    prefix are present; a lone `_1` stays a single-end sample with its full name.
//! 4. Anything else is a single-end sample named after the file.
//!
//! Files of one sample and mate in several lanes (or chunks) become a comma-separated list in
//! lane order, the sheet format of [`crate::sheet`]. Mate 1 and mate 2 lists must cover the same
//! lanes; when they do not, the sample is kept with a problem message.
//!
//! Sample names follow the sheet rules (a letter first, then letters, digits, `_`, `.`, `-`):
//! other characters become `_`, a name that does not start with a letter gets an `X` in front
//! (as R does), a reserved output name gets `_sample` appended, and repeated names get `-2`,
//! `-3`, ... The original prefix is kept in [`SampleRow::from_name`] so the change can be shown.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use crate::sheet::reserved_output_name;

/// Options of [`group_fastq_names`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GroupOptions {
    /// Treat a trailing `_1`/`_2` (or `.1`/`.2`) as mate 1 and mate 2 (rule 3). Default: true.
    pub numeric_mates: bool,
}

impl Default for GroupOptions {
    fn default() -> Self {
        GroupOptions {
            numeric_mates: true,
        }
    }
}

/// Which rule paired a sample's files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rule {
    Illumina,
    ReadTag,
    Number,
    Single,
}

impl Rule {
    /// Short name used in JSON and messages.
    pub fn as_str(self) -> &'static str {
        match self {
            Rule::Illumina => "illumina",
            Rule::ReadTag => "read-tag",
            Rule::Number => "number",
            Rule::Single => "single",
        }
    }
}

/// One sample built from file names.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleRow {
    /// Sample name that passes the sheet rules.
    pub name: String,
    /// The name part of the files before it was made safe (and before the Illumina `_S<n>` was dropped).
    pub from_name: String,
    /// Mate 1 files (or the single-end files) in lane order, as given.
    pub fastq1: Vec<String>,
    /// Mate 2 files in lane order; empty for a single-end sample.
    pub fastq2: Vec<String>,
    pub rule: Rule,
    /// Problems that stop the row from being a valid sheet row as it is (for example mate 1 and
    /// mate 2 cover different lanes). Empty when the row is fine.
    pub problems: Vec<String>,
    /// Things done to the row that the user may want to check (renames, a mate 2 not found).
    pub notes: Vec<String>,
}

/// A file that was left out, and why.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Skipped {
    pub path: String,
    pub reason: String,
}

/// The result of [`group_fastq_names`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Grouping {
    /// Samples in natural order of their names (`S2` before `S10`, case ignored).
    pub samples: Vec<SampleRow>,
    pub skipped: Vec<Skipped>,
}

/// A file of a group: (lane, chunk, path).
type LaneFile = (Option<u32>, Option<u32>, String);

/// A file name read by the rules.
#[derive(Debug, Clone, PartialEq)]
struct Parsed {
    /// Sample part (the group key), for Illumina including `_S<n>`.
    prefix: String,
    /// Illumina sample part without `_S<n>` (the preferred name).
    short: Option<String>,
    mate: Option<u8>,
    lane: Option<u32>,
    chunk: Option<u32>,
    rule: Rule,
}

/// The file name of a path: the part after the last `/` or `\`.
pub fn base_name(path: &str) -> &str {
    path.rsplit(['/', '\\']).next().unwrap_or(path)
}

/// The name without a FASTQ extension (`.fastq`, `.fq`, optionally `.gz` or `.zst`; any case), or `None`.
pub fn fastq_stem(name: &str) -> Option<&str> {
    let lower = name.to_ascii_lowercase();
    for ext in [
        ".fastq.gz",
        ".fq.gz",
        ".fastq.zst",
        ".fq.zst",
        ".fastq",
        ".fq",
    ] {
        if lower.ends_with(ext) && name.len() > ext.len() {
            return Some(&name[..name.len() - ext.len()]);
        }
    }
    None
}

/// `s` ends with `sep` + `tag` + digits (exactly `width` digits when `width` > 0); returns (prefix, number).
fn strip_number_suffix<'a>(
    s: &'a str,
    seps: &[char],
    tag: &str,
    width: usize,
) -> Option<(&'a str, u32)> {
    let digits = s.len() - s.trim_end_matches(|c: char| c.is_ascii_digit()).len();
    if digits == 0 || (width > 0 && digits != width) || digits > 9 {
        return None;
    }
    let (head, num) = s.split_at(s.len() - digits);
    let head = head.strip_suffix(tag)?;
    let sep = head.chars().last()?;
    if !seps.contains(&sep) {
        return None;
    }
    let prefix = &head[..head.len() - sep.len_utf8()];
    if prefix.is_empty() {
        return None;
    }
    Some((prefix, num.parse().ok()?))
}

/// Read one file name stem with rules 1 to 3. `Err` is a skip reason.
fn parse_stem(stem: &str, opts: &GroupOptions) -> Result<Parsed, String> {
    // Rule 1: Illumina <sample>_S<n>[_L<lane>]_<R|I><m>_<chunk>
    if let Some((rest, chunk)) = strip_number_suffix(stem, &['_'], "", 3) {
        let read = ["R1", "R2", "I1", "I2"].iter().find_map(|t| {
            rest.strip_suffix(t)
                .and_then(|h| h.strip_suffix('_'))
                .map(|h| (h, *t))
        });
        if let Some((head, tag)) = read {
            let (head, lane) = match strip_number_suffix(head, &['_'], "L", 3) {
                Some((h, l)) => (h, Some(l)),
                None => (head, None),
            };
            if tag.starts_with('I') {
                return Err(format!("index read ({tag}), not a sample read"));
            }
            // bcl2fastq / BCL Convert names have _S<n>; older CASAVA names have the index sequence instead.
            let short = strip_number_suffix(head, &['_'], "S", 0).map(|(s, _)| s);
            let sample = short.unwrap_or(head);
            if sample.eq_ignore_ascii_case("Undetermined")
                || sample.to_ascii_lowercase().starts_with("undetermined_")
            {
                return Err("Undetermined: reads with no sample barcode".into());
            }
            return Ok(Parsed {
                prefix: head.to_string(),
                short: short.map(str::to_string),
                mate: Some(if tag == "R1" { 1 } else { 2 }),
                lane,
                chunk: Some(chunk),
                rule: Rule::Illumina,
            });
        }
    }
    // Rule 2: read tag, optionally after a lane.
    for tag in ["R", "r", "read", "Read", "READ"] {
        if let Some((head, m)) = strip_number_suffix(stem, &['_', '.', '-'], tag, 1) {
            if m == 1 || m == 2 {
                let (head, lane) = match strip_number_suffix(head, &['_'], "L", 3) {
                    Some((h, l)) => (h, Some(l)),
                    None => (head, None),
                };
                return Ok(Parsed {
                    prefix: head.to_string(),
                    short: None,
                    mate: Some(m as u8),
                    lane,
                    chunk: None,
                    rule: Rule::ReadTag,
                });
            }
        }
    }
    // Rule 3: number.
    if opts.numeric_mates {
        if let Some((head, m)) = strip_number_suffix(stem, &['_', '.'], "", 1) {
            if m == 1 || m == 2 {
                return Ok(Parsed {
                    prefix: head.to_string(),
                    short: None,
                    mate: Some(m as u8),
                    lane: None,
                    chunk: None,
                    rule: Rule::Number,
                });
            }
        }
    }
    Ok(single(stem))
}

fn single(stem: &str) -> Parsed {
    Parsed {
        prefix: stem.to_string(),
        short: None,
        mate: None,
        lane: None,
        chunk: None,
        rule: Rule::Single,
    }
}

/// Compare strings with digit runs compared as numbers (`S2` < `S10`) and letters compared
/// without case (`dcat` < `WT`), then as plain text so the order is total.
pub fn natural_cmp(a: &str, b: &str) -> Ordering {
    let (mut x, mut y) = (a.as_bytes(), b.as_bytes());
    loop {
        match (x.first(), y.first()) {
            (None, None) => return a.cmp(b),
            (None, _) => return Ordering::Less,
            (_, None) => return Ordering::Greater,
            (Some(c), Some(d)) if c.is_ascii_digit() && d.is_ascii_digit() => {
                let n = x.iter().take_while(|c| c.is_ascii_digit()).count();
                let m = y.iter().take_while(|c| c.is_ascii_digit()).count();
                let (p, q) = (trim_zeros(&x[..n]), trim_zeros(&y[..m]));
                let ord = p.len().cmp(&q.len()).then_with(|| p.cmp(q));
                if ord != Ordering::Equal {
                    return ord;
                }
                x = &x[n..];
                y = &y[m..];
            }
            (Some(c), Some(d)) => {
                let ord = c.to_ascii_lowercase().cmp(&d.to_ascii_lowercase());
                if ord != Ordering::Equal {
                    return ord;
                }
                x = &x[1..];
                y = &y[1..];
            }
        }
    }
}

fn trim_zeros(s: &[u8]) -> &[u8] {
    let z = s.iter().take_while(|&&c| c == b'0').count();
    &s[z.min(s.len().saturating_sub(1))..]
}

/// A sample name that passes the sheet rules, made from any text (see the module notes).
pub fn safe_sample_name(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len() + 1);
    for c in raw.chars() {
        let c = if c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-') {
            c
        } else {
            '_'
        };
        if c == '_' && out.ends_with('_') {
            continue;
        }
        out.push(c);
    }
    let trimmed = out.trim_matches(|c| matches!(c, '_' | '.' | '-'));
    let mut name = if trimmed.is_empty() {
        "sample".to_string()
    } else {
        trimmed.to_string()
    };
    if !name.starts_with(|c: char| c.is_ascii_alphabetic()) {
        name.insert(0, 'X');
    }
    if reserved_output_name(&name).is_some() {
        name.push_str("_sample");
    }
    name
}

/// Problems of a fastq path that the sheet format cannot hold.
fn path_problem(path: &str) -> Option<String> {
    if path.chars().any(|c| c.is_control()) {
        return Some("the path has a tab, newline or other control character".into());
    }
    if path.contains(',') {
        return Some("the path has a comma (the sheet uses commas to list lanes)".into());
    }
    None
}

/// Group FASTQ paths (or bare file names) into sample rows. Non-FASTQ files, index reads,
/// `Undetermined` files and repeated paths are listed in [`Grouping::skipped`].
pub fn group_fastq_names<S: AsRef<str>>(paths: &[S], opts: &GroupOptions) -> Grouping {
    let mut skipped = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut files: Vec<(String, Parsed)> = Vec::new();
    for path in paths {
        let path = path.as_ref().trim_end_matches(['\r', '\n']);
        if path.trim().is_empty() {
            continue;
        }
        if !seen.insert(path.to_string()) {
            skipped.push(Skipped {
                path: path.into(),
                reason: "listed twice".into(),
            });
            continue;
        }
        if let Some(reason) = path_problem(path) {
            skipped.push(Skipped {
                path: path.into(),
                reason,
            });
            continue;
        }
        let Some(stem) = fastq_stem(base_name(path)) else {
            skipped.push(Skipped {
                path: path.into(),
                reason: "not a FASTQ file name (.fastq or .fq, optionally .gz or .zst)".into(),
            });
            continue;
        };
        match parse_stem(stem, opts) {
            Ok(p) => files.push((path.to_string(), p)),
            Err(reason) => skipped.push(Skipped {
                path: path.into(),
                reason,
            }),
        }
    }

    // Rule 3 needs both mates: a lone _1 or _2 is a single-end sample with its full name.
    let mut number_mates: BTreeMap<String, [bool; 2]> = BTreeMap::new();
    for (_, p) in files.iter().filter(|(_, p)| p.rule == Rule::Number) {
        number_mates.entry(p.prefix.clone()).or_default()[p.mate.unwrap() as usize - 1] = true;
    }
    for (path, p) in files.iter_mut() {
        if p.rule == Rule::Number && number_mates[&p.prefix] != [true, true] {
            *p = single(fastq_stem(base_name(path)).unwrap());
        }
    }

    // Group by (prefix, paired-or-not). A single-end file never joins a paired group.
    struct Group {
        prefix: String,
        short: Option<String>,
        rule: Rule,
        mates: [Vec<LaneFile>; 3], // [single, mate 1, mate 2]
    }
    let mut groups: BTreeMap<(String, bool), Group> = BTreeMap::new();
    for (path, p) in files {
        let key = (p.prefix.clone(), p.mate.is_some());
        let g = groups.entry(key).or_insert_with(|| Group {
            prefix: p.prefix.clone(),
            short: p.short.clone(),
            rule: p.rule,
            mates: Default::default(),
        });
        // Mixed rules in one group (e.g. x_R1 and x_1): keep the first rule's name.
        g.mates[p.mate.map_or(0, usize::from)].push((p.lane, p.chunk, path));
    }

    let mut rows: Vec<(String, SampleRow)> = Vec::new(); // (full prefix, row)
    for g in groups.into_values() {
        let order = |v: &mut Vec<LaneFile>| {
            v.sort_by(|a, b| {
                (a.0, a.1)
                    .cmp(&(b.0, b.1))
                    .then_with(|| natural_cmp(&a.2, &b.2))
            });
        };
        let [mut s, mut m1, mut m2] = g.mates;
        order(&mut s);
        order(&mut m1);
        order(&mut m2);
        let mut problems = Vec::new();
        let mut notes = Vec::new();
        let (fastq1, fastq2, rule) = if !s.is_empty() {
            if s.len() > 1 {
                notes.push(format!("{} files with the same name in different folders are listed as lanes; check them", s.len()));
            }
            (s, Vec::new(), Rule::Single)
        } else if m1.is_empty() {
            problems.push("only mate 2 files were found: add the mate 1 files, or change the mate 2 cell to the mate 1 cell for single-end reads".into());
            (m1, m2, g.rule)
        } else {
            if m2.is_empty() {
                notes.push("no mate 2 file was found: single-end".into());
            } else {
                let lanes = |v: &[LaneFile]| v.iter().map(|f| (f.0, f.1)).collect::<Vec<_>>();
                if m1.len() != m2.len() {
                    problems.push(format!(
                        "{} mate 1 files but {} mate 2 files: a file is missing",
                        m1.len(),
                        m2.len()
                    ));
                } else if lanes(&m1) != lanes(&m2) {
                    problems
                        .push("mate 1 and mate 2 files are of different lanes or chunks".into());
                }
            }
            (m1, m2, g.rule)
        };
        if fastq1.len() + fastq2.len() > 2 && problems.is_empty() && rule != Rule::Single {
            notes.push(format!(
                "{} lanes or chunks per mate, listed in lane order",
                fastq1.len().max(fastq2.len())
            ));
        }
        rows.push((
            g.prefix.clone(),
            SampleRow {
                name: String::new(),
                from_name: g.short.clone().unwrap_or_else(|| g.prefix.clone()),
                fastq1: fastq1.into_iter().map(|f| f.2).collect(),
                fastq2: fastq2.into_iter().map(|f| f.2).collect(),
                rule,
                problems,
                notes,
            },
        ));
    }

    // Names: the Illumina short name unless two samples would share it, made safe and unique.
    let mut short_count: BTreeMap<String, usize> = BTreeMap::new();
    for (_, r) in &rows {
        *short_count.entry(r.from_name.clone()).or_default() += 1;
    }
    let mut used: BTreeMap<String, usize> = BTreeMap::new();
    for (full_prefix, r) in rows.iter_mut() {
        let raw = if short_count[&r.from_name] > 1 && r.rule == Rule::Illumina {
            r.from_name = full_prefix.clone();
            full_prefix.clone()
        } else {
            r.from_name.clone()
        };
        let base = safe_sample_name(&raw);
        let n = used.entry(base.to_ascii_lowercase()).or_default();
        *n += 1;
        r.name = if *n == 1 { base } else { format!("{base}-{n}") };
        if r.name != raw {
            r.notes.push(format!(
                "sample name changed from \"{raw}\" to \"{}\"",
                r.name
            ));
        }
    }
    let mut samples: Vec<SampleRow> = rows.into_iter().map(|(_, r)| r).collect();
    samples.sort_by(|a, b| natural_cmp(&a.name, &b.name));
    Grouping { samples, skipped }
}

/// A sample sheet (the format of [`crate::sheet`]) of the rows, with empty covariate columns
/// `columns` added for the user to fill.
pub fn grouping_sheet_tsv(rows: &[SampleRow], columns: &[String]) -> String {
    let mut out = String::from("sample\tfastq1\tfastq2");
    for c in columns {
        out.push('\t');
        out.push_str(c);
    }
    out.push('\n');
    for r in rows {
        out.push_str(&r.name);
        out.push('\t');
        out.push_str(&r.fastq1.join(","));
        out.push('\t');
        out.push_str(&r.fastq2.join(","));
        for _ in columns {
            out.push('\t');
        }
        out.push('\n');
    }
    out
}
