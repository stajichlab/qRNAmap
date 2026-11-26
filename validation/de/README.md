# DE checks against DESeq2

The page's differential expression (`site/analyze/de-stats.js`) is a quick screen, not DESeq2. These scripts compare
it, and later a Rust port, with R DESeq2 (version pinned in `pixi.toml`) on the same count matrix. Run on the HPCC.

- `deseq2_reference.R COUNTS COLDATA DESIGN FACTOR LEVEL REF OUT_PREFIX`: DESeq2 defaults; writes `OUT_PREFIX.results.tsv`
  and `OUT_PREFIX.sizefactors.tsv`.
- `page_de.mjs COUNTS COLDATA FACTOR LEVEL REF OUT.tsv`: the page's test on the same input, same columns.
- `compare_de.py REFERENCE.tsv OTHER.tsv`: calls at padj < 0.05 in both or one, and correlations of the log2 fold
  change and -log10 p.
- `run_compare.sh COUNTS COLDATA FACTOR LEVEL REF OUTDIR [DESIGN...]`: all three, one DESeq2 run per design.

COLDATA is a TSV with a `sample` column (the count matrix's column names) and one column per covariate.
