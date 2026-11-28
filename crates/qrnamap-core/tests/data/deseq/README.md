# DESeq2 port fixtures

Reference data for the Rust port of DESeq2 (`docs/superpowers/specs/2026-11-27-deseq2-rust-design.md`). Each directory
holds one DESeq2 run (R DESeq2 1.50.2 with default arguments) on a small count matrix: the inputs and every
intermediate, so `cargo test` can check each stage of the port without R. `tests/deseq_fixtures.rs` checks that every
fixture is complete and consistent.

Made on the HPCC by `validation/de/make_deseq2_fixtures.sh` (inputs from `validation/de/build_deseq2_fixtures.R`, dumps
from `validation/de/deseq2_dump.R`). Do not edit them by hand: re-run the script.

## Files of a fixture

| file | content |
|---|---|
| `counts.tsv` | `gene`, then one integer column per sample (the matrix DESeq2 got) |
| `coldata.tsv` | `sample` plus covariates |
| `spec.tsv` | design, contrast (`factor,numerator,denominator`), ref (the factor's reference level), purpose |
| `run.tsv` | key, value: factor levels as used (C-locale order, ref first), m, p, coefficient names, trend type and coefficients (dispersion = asymptDisp + extraPois / mean), varLogDispEsts, dispPriorVar and its path (`trigamma` or `simulation`), Cook's cutoff, filter threshold and theta, case counts |
| `genes.tsv` | per gene: `mcols(dds)` (baseMean, baseVar, allZero, dispGeneEst, dispGeneIter, dispFit, dispersion, dispIter, dispOutlier, dispMAP, coefficients and SE_, WaldStatistic_, WaldPvalue_ on the log2 scale, betaConv, betaIter, deviance, maxCooks, replace), `res_*` (`results()` columns) and `optimRow` (IRLS reached maxit, so DESeq2 used its optim fallback) |
| `sizefactors.tsv` | sample, sizeFactor |
| `hat.tsv`, `cooks.tsv` | gene × sample hat diagonals and Cook's distances |
| `filter.tsv` | independent filtering: theta, numRej, lowess fit |
| `messages.txt` | DESeq2's messages and warnings |

Numbers have 17 significant digits; NA is `NA`.

## Data

- **yeast_\***: Schurch et al. 2016, *RNA* 22:839, *S. cerevisiae* WT vs *snf2*, 48 replicates each (htseq counts and
  the bad-replicate list from github.com/bartongroup/profDGE48, MIT licence). Replicates chosen in a seeded order.
  The batch, numeric-covariate and 3-level groupings are synthetic labels on real counts; the 7 vs 7 fixture has
  planted outliers (named in its `spec.tsv`).
- **pasilla_\***: Brooks et al. 2011, *Genome Res* 21:193, *D. melanogaster* pasilla knock-down (Bioconductor data
  package pasilla 1.42.0, LGPL). Library type renamed `single`/`paired`.
- **afum_\***: *A. fumigatus* AF293 and CEA10, Azole vs GMM, kallisto counts streamed by the qRNAmap page
  (`validation/de/data/`, see `validation/de/README.md`).
- **sim_trend_fails**: simulated (seed 5), dispersion rising with the mean, so the parametric trend fails and DESeq2
  substitutes its local fit; the port must stop.

Each fixture is 300 seeded random genes plus all-zero genes, genes zero in one contrasted group, genes zero in both
contrasted levels (when present), and the genes that the full-data run of the same case sent to the optim fallback or
flagged by Cook's distance. DESeq2 ran on exactly this subset, so a fixture is self-consistent. The full-data runs stay
on the HPCC (`bench/2026-11-27_deseq2-fixtures/dumps/full/`) for the final validation.

## The prior table

`src/deseq/prior_hist.tsv` (made by `validation/de/deseq2_prior_table.R`) holds the simulated histograms of
`estimateDispersionsPriorVar` for m − p = 1, 2, 3: `df`, `i`, `x` (`seq(0, 8, length = 200)`), `n`, `c1..c40` (counts in
the bins of `hist(breaks = -20:20/2)`; density = c / (n · 0.5)). `validation/de/check_prior_table.R` recomputes
DESeq2's `dispPriorVar` from it on every fixture that takes this path and requires DESeq2's value.
