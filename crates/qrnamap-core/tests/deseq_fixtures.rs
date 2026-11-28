//! The DESeq2 port's fixtures (`tests/data/deseq/<name>/`, made on the HPCC by `validation/de/make_deseq2_fixtures.sh`)
//! and the m - p <= 3 prior table (`src/deseq/prior_hist.tsv`) are complete and agree with each other: the same genes in
//! the same order in every per-gene file, the same samples everywhere, the spec and the run record agree. The port's
//! stage tests (phases 2 to 4) read these files; this test makes a broken or partly copied fixture fail early.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Rows of a TSV file, split on tabs (header included).
fn tsv(path: &Path) -> Vec<Vec<String>> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines().filter(|l| !l.is_empty()).map(|l| l.split('\t').map(str::to_string).collect()).collect()
}

fn key_values(path: &Path) -> HashMap<String, String> {
    tsv(path).into_iter().map(|r| (r[0].clone(), r.get(1).cloned().unwrap_or_default())).collect()
}

fn column(rows: &[Vec<String>], name: &str) -> Vec<String> {
    let i = rows[0].iter().position(|h| h == name).unwrap_or_else(|| panic!("no column {name}"));
    rows[1..].iter().map(|r| r[i].clone()).collect()
}

fn fixtures() -> Vec<PathBuf> {
    let dir = root().join("tests/data/deseq");
    let mut v: Vec<PathBuf> = fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .map(|e| e.unwrap().path())
        .filter(|p| p.is_dir())
        .collect();
    v.sort();
    v
}

#[test]
fn every_fixture_is_complete_and_consistent() {
    let all = fixtures();
    assert!(all.len() >= 16, "expected the 16 fixtures of phase 1, found {}", all.len());
    for d in &all {
        let name = d.file_name().unwrap().to_string_lossy().to_string();
        let spec = key_values(&d.join("spec.tsv"));
        let run = key_values(&d.join("run.tsv"));
        assert_eq!(run["deseq2_version"], "1.50.2", "{name}");
        for k in ["design", "contrast", "ref"] {
            assert_eq!(spec[k], run[k], "{name}: {k}");
        }

        let counts = tsv(&d.join("counts.tsv"));
        let samples: Vec<String> = counts[0][1..].to_vec();
        let genes: Vec<String> = counts[1..].iter().map(|r| r[0].clone()).collect();
        for r in &counts[1..] {
            assert_eq!(r.len(), samples.len() + 1, "{name}: ragged counts row {}", r[0]);
            assert!(r[1..].iter().all(|c| c.parse::<u64>().is_ok()), "{name}: non-integer count in {}", r[0]);
        }
        assert_eq!(run["samples"].split(',').collect::<Vec<_>>(), samples, "{name}: run samples");
        assert_eq!(run["m"].parse::<usize>().unwrap(), samples.len(), "{name}: m");
        assert_eq!(run["genes"].parse::<usize>().unwrap(), genes.len(), "{name}: genes");
        assert_eq!(column(&tsv(&d.join("coldata.tsv")), "sample"), samples, "{name}: coldata samples");
        assert_eq!(column(&tsv(&d.join("sizefactors.tsv")), "sample"), samples, "{name}: size factor samples");

        let g = tsv(&d.join("genes.tsv"));
        assert_eq!(column(&g, "gene"), genes, "{name}: genes.tsv genes");
        for c in ["baseMean", "dispGeneEst", "dispFit", "dispMAP", "dispersion", "betaIter", "res_log2FoldChange", "res_padj", "optimRow"] {
            column(&g, c);
        }
        let p = run["p"].parse::<usize>().unwrap();
        let coefs: Vec<&str> = run["coefficients"].split(',').collect();
        assert_eq!(coefs.len(), p, "{name}: coefficients");
        for c in &coefs {
            column(&g, c);
            column(&g, &format!("SE_{c}"));
        }
        for f in ["hat.tsv", "cooks.tsv"] {
            let m = tsv(&d.join(f));
            assert_eq!(m[0][1..].to_vec(), samples, "{name}: {f} samples");
            assert_eq!(column(&m, "gene"), genes, "{name}: {f} genes");
        }
        let filter = tsv(&d.join("filter.tsv"));
        assert_eq!(filter[0], ["theta", "numRej", "lowess"], "{name}: filter.tsv header");
    }
}

#[test]
fn the_prior_table_has_three_complete_blocks() {
    let t = tsv(&root().join("src/deseq/prior_hist.tsv"));
    assert_eq!(t[0].len(), 44);
    assert_eq!(t.len(), 1 + 3 * 200);
    for (k, r) in t[1..].iter().enumerate() {
        assert_eq!(r[0].parse::<usize>().unwrap(), 1 + k / 200);
        assert_eq!(r[1].parse::<usize>().unwrap(), 1 + k % 200);
        let n: u64 = r[3].parse().unwrap();
        let sum: u64 = r[4..].iter().map(|c| c.parse::<u64>().unwrap()).sum();
        assert_eq!(sum, n, "row {k}: bin counts add up to n");
        assert!(n > 9_000 && n <= 10_000, "row {k}: n {n}");
    }
}
