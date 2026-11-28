# DE checks against DESeq2

The page's differential expression (`site/analyze/de-stats.js`) is a quick screen, not DESeq2. These scripts compare
it, and later a Rust port, with R DESeq2 (version pinned in `pixi.toml`) on the same count matrix. R comes from the pixi environment `de` (Linux and macOS arm64).

- `deseq2_reference.R COUNTS COLDATA DESIGN FACTOR LEVEL REF OUT_PREFIX`: DESeq2 defaults; writes `OUT_PREFIX.results.tsv`
  and `OUT_PREFIX.sizefactors.tsv`.
- `page_de.mjs COUNTS COLDATA FACTOR LEVEL REF OUT.tsv`: the page's test on the same input, same columns.
- `compare_de.py REFERENCE.tsv OTHER.tsv`: calls at padj < 0.05 in both or one, and correlations of the log2 fold
  change and -log10 p.
- `run_compare.sh COUNTS COLDATA FACTOR LEVEL REF OUTDIR [DESIGN...]`: all three, one DESeq2 run per design.

COLDATA is a TSV with a `sample` column (the count matrix's column names) and one column per covariate.

## Reference data for the Rust port (phase 1)

- `make_deseq2_fixtures.sh WORK_DIR` (HPCC, in a job): fetches the data, builds the inputs (`build_deseq2_fixtures.R`),
  runs DESeq2 on the full data, adds the genes those runs flagged to the fixtures (`add_flagged_genes.R`), runs DESeq2 on
  every fixture with all intermediates saved (`deseq2_dump.R`), makes the m − p ≤ 3 prior table
  (`deseq2_prior_table.R`) and checks it against DESeq2 (`check_prior_table.R`). Writes
  `crates/qrnamap-core/tests/data/deseq/` (format and data sources in its README) and
  `crates/qrnamap-core/src/deseq/prior_hist.tsv`.
- `test_deseq2_dump.sh`: fixture test of the dumper, the prior table and its check on simulated counts (needs R).
- `data/afum_azole_gmm.*`: the *A. fumigatus* counts from the page (2026-11-26), the only input not downloadable.

## Results

### 2026-11-27: phase 1 reference runs (DESeq2 1.50.2, R 4.5.3, HPCC job 32520679)
Per-run cases in `deseq2_fixtures_summary.tsv` (16 fixtures, 6 full-data runs: yeast Schurch 3 vs 3 and all 42 vs 44
good replicates, pasilla, *A. fumigatus* with and without strain).
- The prior table reproduces DESeq2's `dispPriorVar` exactly on all six fixtures with m − p ≤ 3 (four at the 0.25 floor,
  0.36036… and 0.26426… above it).
- The parametric dispersion trend fitted on all 21 real-data runs; only the simulation built to fail it switched to the
  local fit.
- Cases covered by the fixtures: GLM optim fallback rows (*A. fumigatus* 7 and 2, pasilla 5), Cook's-flagged genes (18
  and 1), Cook's replacement (5 planted outliers in 7 vs 7), no filtering (null WT vs WT, 1 vs 2, 1 vs 3), dispersion
  grid genes in every fixture and dispersion outliers in most, all-zero genes, genes zero in one or both contrasted levels.
- Full-data *A. fumigatus* runs give the same calls as on 2026-11-26 (1,261 for `~ condition`, 1,390 with strain).


### 2026-11-26: *A. fumigatus*, full-depth streamed counts, 7 samples
Samples: AF293 Azole ×2, AF293 GMM ×3, CEA10 Azole ×2 (no CEA10 GMM). Azole vs GMM; DESeq2 1.50.2 (pixi `de`, osx-arm64).
The page's old test (linear-scale SE) reproduces the downloaded page output (|stat| differences ≤ 5e-5) and calls 0 genes.

| method | padj < 0.05 | agreement with DESeq2 `~ condition` (7,543 genes with a DESeq2 padj) |
|---|---|---|
| page, old test | 0 | none |
| page, `de-stats.js` (moderated t on log2) | 164 | 113 shared, 1 page-only, 1,148 DESeq2-only; LFC Spearman 0.98, -log10 p Spearman 0.92 |
| DESeq2 `~ condition` | 1,261 | |
| DESeq2 `~ strain + condition` | 1,390 | the page's test agrees less with this one (LFC Spearman 0.73) |

The page's test ranks genes like DESeq2 and makes almost no calls DESeq2 does not, but it is far less powerful (about
one call in eight), and it cannot model strain. It is a screen; the DESeq2 port is the fix.
