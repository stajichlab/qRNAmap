//! Design formula grammar and contrast parsing.
//!
//! Grammar: `~` term (op term)*, where op is `+`, `:` or `*` and term is a column name
//! (`[A-Za-z][A-Za-z0-9_.]*`) present in the sheet columns. Space and tab around tokens are
//! ignored. Everything else is a `BadSyntax` error. [`Design::to_r`] rebuilds the formula from
//! the validated names and operators only, so the text is safe to paste into R code.

use crate::sheet::ColType;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub enum DesignError {
    /// A term is a valid name but is not a column of the sheet.
    UnknownColumn(String),
    /// The text breaks the grammar. The message describes the problem and echoes at most one
    /// escaped character.
    BadSyntax(String),
}

impl fmt::Display for DesignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DesignError::UnknownColumn(c) => write!(f, "unknown column {c:?} in design"),
            DesignError::BadSyntax(m) => write!(f, "bad syntax: {m}"),
        }
    }
}

impl std::error::Error for DesignError {}

fn bad(msg: &str) -> DesignError {
    DesignError::BadSyntax(msg.to_string())
}

fn bad_char(c: char) -> DesignError {
    DesignError::BadSyntax(format!("unexpected character '{}'", c.escape_debug()))
}

/// A validated design formula.
#[derive(Debug, Clone, PartialEq)]
pub struct Design {
    terms: Vec<String>,
    /// `ops[i]` joins `terms[i]` and `terms[i + 1]`.
    ops: Vec<char>,
}

impl Design {
    pub fn parse(formula: &str, columns: &[(String, ColType)]) -> Result<Design, DesignError> {
        let chars: Vec<char> = formula.chars().collect();
        let mut i = 0;
        let skip_ws = |i: &mut usize| {
            while *i < chars.len() && matches!(chars[*i], ' ' | '\t') {
                *i += 1;
            }
        };
        skip_ws(&mut i);
        match chars.get(i) {
            Some('~') => i += 1,
            Some(_) => return Err(bad("formula must start with '~'")),
            None => return Err(bad("empty formula")),
        }
        let mut terms: Vec<String> = Vec::new();
        let mut ops: Vec<char> = Vec::new();
        loop {
            skip_ws(&mut i);
            match chars.get(i) {
                Some(c) if c.is_ascii_alphabetic() => {
                    let start = i;
                    while i < chars.len()
                        && (chars[i].is_ascii_alphanumeric() || matches!(chars[i], '_' | '.'))
                    {
                        i += 1;
                    }
                    terms.push(chars[start..i].iter().collect());
                }
                Some('+' | ':' | '*' | '~') => return Err(bad("expected a column name")),
                Some(c) => return Err(bad_char(*c)),
                None => return Err(bad("formula ends where a column name is expected")),
            }
            skip_ws(&mut i);
            match chars.get(i) {
                None => break,
                Some(&c @ ('+' | ':' | '*')) => {
                    ops.push(c);
                    i += 1;
                }
                Some(c) if c.is_ascii_alphanumeric() || matches!(c, '_' | '.') => {
                    return Err(bad("two column names without an operator"))
                }
                Some(c) => return Err(bad_char(*c)),
            }
        }
        for t in &terms {
            if !columns.iter().any(|(c, _)| c == t) {
                return Err(DesignError::UnknownColumn(t.clone()));
            }
        }
        Ok(Design { terms, ops })
    }

    /// Column names in formula order.
    pub fn terms(&self) -> &[String] {
        &self.terms
    }

    /// R formula text. `+` is written with spaces (`~ a + b`); `:` and `*` are written
    /// without (`a:b`, `a*b`). Only validated names and `~ + : *` appear.
    pub fn to_r(&self) -> String {
        let mut s = String::from("~ ");
        for (k, t) in self.terms.iter().enumerate() {
            if k > 0 {
                match self.ops[k - 1] {
                    '+' => s.push_str(" + "),
                    op => s.push(op),
                }
            }
            s.push_str(t);
        }
        s
    }
}

/// A DE contrast: `level` versus `reference` within `factor`.
#[derive(Debug, Clone, PartialEq)]
pub struct Contrast {
    pub factor: String,
    pub level: String,
    pub reference: String,
}

fn is_level(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
}

impl Contrast {
    /// Parse `factor,level,reference`. The factor must be a `Factor` column. Level and
    /// reference match `[A-Za-z0-9_.-]+` and must differ.
    pub fn parse(text: &str, columns: &[(String, ColType)]) -> Result<Contrast, DesignError> {
        let parts: Vec<&str> = text.split(',').map(str::trim).collect();
        if parts.len() != 3 {
            return Err(bad("contrast must be 'factor,level,reference'"));
        }
        let (factor, level, reference) = (parts[0], parts[1], parts[2]);
        let valid_name = {
            let mut it = factor.chars();
            matches!(it.next(), Some(c) if c.is_ascii_alphabetic())
                && it.all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.'))
        };
        if !valid_name {
            return Err(bad("contrast factor is not a valid column name"));
        }
        match columns.iter().find(|(c, _)| c == factor) {
            None => return Err(DesignError::UnknownColumn(factor.to_string())),
            Some((_, ColType::Numeric)) => return Err(bad("contrast factor is a numeric column")),
            Some((_, ColType::Factor)) => {}
        }
        if !is_level(level) {
            return Err(bad("contrast level is empty or has a character outside [A-Za-z0-9_.-]"));
        }
        if !is_level(reference) {
            return Err(bad(
                "contrast reference is empty or has a character outside [A-Za-z0-9_.-]",
            ));
        }
        if level == reference {
            return Err(bad("contrast level and reference are the same"));
        }
        Ok(Contrast { factor: factor.into(), level: level.into(), reference: reference.into() })
    }
}
