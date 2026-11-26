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

## Results

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
