# Notes on parity_results.tsv (run 44944424373)

- DESeq2 row: 1 of 38,203 numeric values is out of tolerance. The table has 6,477 x 6 = 38,862 cells; 659 are NA
  on both sides. Native DESeq2 1.50.2 with R 4.5.3, WASM DESeq2 1.52.0 with R 4.6.0. No same-version native
  reference was run, so the failure cannot be attributed to WASM or to the version difference.
