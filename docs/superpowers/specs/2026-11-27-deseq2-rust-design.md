# DESeq2 in Rust: design

Date: 2026-11-27. Status: draft for approval (revised after an independent review the same day, see "Review").
Parent: `2026-11-19-qrnamap-design.md`. Evidence: `validation/de/README.md` (the page's screen against DESeq2 on
*A. fumigatus*, 2026-11-26). Reference source: DESeq2 1.50.2 (Bioconductor RELEASE_3_22: `R/core.R`,
`R/fitNbinomGLMs.R`, `R/results.R`, `R/methods.R`, `src/DESeq2.cpp`); function names below are DESeq2's.

## Intent

The browser page and the desktop app give differential expression from `site/analyze/de-stats.js`: median-of-ratios size
factors, then a limma-style moderated t-test on log2 counts. On *A. fumigatus* (7 samples, Azole vs GMM) it ranks genes like
DESeq2 (log2 fold change Spearman 0.98) but calls 164 genes where DESeq2 calls 1,261, and it cannot model a second factor.
The command-line tool runs real DESeq2 through `Rscript`, which needs R and Bioconductor.

Goal: a Rust implementation of the DESeq2 method (1.50.2, `DESeq()` and `results(contrast = ...)` with default arguments)
in `qrnamap-core`, compiled to WebAssembly for the page and the desktop app, matching R DESeq2 stage by stage within the
tolerances below. The R route stays the reference and the default of `qrnamap de` until the port passes validation.

The port never gives a quiet approximation. Where it cannot do what DESeq2 does, it stops with a message that names the
command-line R engine, or (only where stated) keeps DESeq2's answer shape and lists the affected genes as a divergence.

Non-goals (version 1): LFC shrinkage, the likelihood ratio test, VST and rlog, `betaPrior = TRUE`, user size or
normalization factors (tximport offsets), `fitType = "local"`.

## The algorithm, step by step

`m` samples, `p` columns of the model matrix, μ floored at `minmu = 0.5` wherever DESeq2 floors it (weights, deviance,
sandwich covariance, hat values, gene-wise dispersion).

1. **Input checks** (DESeq2's errors): non-negative integer counts; no NA in design variables; no empty factor levels;
   m > p; not every gene equal; full column rank (a rank-deficient design names the confounded columns). Genes with all-zero
   counts get baseMean 0 and NA everywhere else.
2. **Size factors** (`estimateSizeFactorsForMatrix`, type `ratio`): per sample, the median over genes with no zero count of
   count / geometric mean. If every gene has a zero, stop with DESeq2's error; `poscounts` is an explicit option only.
3. **Model matrix**: from the existing `Design` parser (`+`, `:`, `*`) and the sheet's column types. Treatment coding;
   the contrast factor's reference level first; other factors' levels in C-locale (byte) order, and the R script that
   `qrnamap de` writes is changed to pass explicit `levels =` for every factor so both engines use the same order (R's
   `sort()` is locale-dependent). Numeric columns enter as they are.
4. **Gene-wise dispersions** (`estimateDispersionsGeneEst`, niter = 1):
   - α_init = min(rough, moments) clamped to [minDisp = 1e-8, maxDisp = max(10, m)]; the rough estimate
     (`roughDispEstimate`) uses linear-model μ floored at 1.
   - μ from `linearModelMu` when the number of distinct model-matrix rows equals p, else from the NB GLM (step 7) with α_init.
   - `fitDisp` (C++): maximize the Cox-Reid adjusted profile likelihood over log α by backtracking line search: κ₀ = 1,
     Armijo ε = 1e-4, log α kept in [−30, 10], κ × 1.1 after an accepted step capped at κ₀, κ halved every 5 accepted steps;
     stop when the gain is below 1e-6 or log α < log(minDisp/10); maxit 100. If the objective did not rise by |lp|/1e6,
     return α_init.
   - Not converged = iterations = maxit **or** iterations = 1; those genes with estimate > 1e-7 get `fitDispGrid`. Clamp.
5. **Trend** (`estimateDispersionsFit`, `parametricDispersionFit`): genes with gene-wise estimate strictly > 100·minDisp;
   start (a₀, a₁) = (0.1, 1); repeat: Gamma-family identity-link GLM (a port of R's `glm.fit` for this case, with start
   values) of dispersion on 1/baseMean, dropping residual ratios outside [1e-4, 15]; accept when Σ log(c/c_old)² < 1e-6
   **and** the GLM converged; more than 10 rounds or a non-positive coefficient is a failure. **On failure the port stops**
   ("the parametric dispersion trend did not fit; DESeq2 uses a local fit here, which this version does not have: use
   `qrnamap de` with the R engine"). It does not fall back to the mean.
6. **Prior variance and MAP** (`estimateDispersionsPriorVar`, `estimateDispersionsMAP`):
   - varLogDispEsts = MAD² of log(gene-wise) − log(trend) over the genes in the fit.
   - m − p > 3: dispPriorVar = max(varLogDispEsts − trigamma((m − p)/2), 0.25).
   - m − p ≤ 3 (and m > p): DESeq2 compares a 40-bin histogram (`brks = -20:20/2`) of the residuals with simulated
     `log(rchisq(1e4, m−p)) + rnorm(1e4, 0, √x) − log(m−p)` for x on 200 points in [0, 8], under `set.seed(2)`
     (deterministic); smooths the KL divergences with `loess(span = 0.2)`; takes the argmin on a 1000-point grid; floors at
     0.25; trigamma is not subtracted. The simulated histograms depend only on m − p, so they are **precomputed in R** for
     m − p = 1, 2, 3 (3 × 200 × 40 values, a data file in the crate, with the script that makes it); the port computes the
     data histogram, the KL divergences, and a port of R's `loess` (degree 2) for the argmin.
   - MAP: start from the gene-wise estimate if it is > 0.1·trend, else the trend; maximize the Cox-Reid likelihood plus a
     normal prior on log α (mean log trend, variance dispPriorVar); converged = iterations < maxit; else the grid with the
     prior. Dispersion outliers: log(gene-wise) > log(trend) + 2·√varLogDispEsts keep their gene-wise value. Clamp.
7. **NB GLM** (`fitNbinomGLMs`, C++ `fitBeta`): IRLS on the natural log scale; start β = least squares (QR) on
   log(normalized count + 0.1); ridge λ = 1e-6 per coefficient **on the log2 scale** (1e-6 / ln(2)² on the natural scale);
   each step solves the QR of [√W·X; √Λ]; converged when |dev − dev_old| / (|dev| + 0.1) < 1e-8 (after the first
   iteration); maxit 100; a row is abandoned (iterations = maxit) when any |β| > 30 or the test is NaN. Rows not converged,
   with NA β or non-positive variance, go to DESeq2's `optim` fallback: L-BFGS-B on the **log2** scale, bounds ±30,
   objective the NB log-likelihood plus `dnorm(β, 0, √(1/λ))`, numerical gradient; its estimate is kept whether or not it
   converged. The port implements a bounded quasi-Newton optimizer on the same objective (not a port of R's `lbfgsb`); these
   rows keep their p-values and are listed as a divergence. Covariance: (X'WX + Λ)⁻¹ X'WX (X'WX + Λ)⁻¹. Coefficients
   reported on the log2 scale. Deviance with R's `dnbinom_mu` (saddle-point form).
8. **Wald test and contrast** (`results`, `getContrast`): stat = c'β / √(c'Σc), two-sided normal p (`useT = FALSE`).
   A numerator that is the reference level flips the sign of the coefficient and stat. A contrast between two non-reference
   levels recomputes Σ from the stored β (μ floored, no iterations), as `getContrast` does. If every count in the two
   contrasted levels is zero (gene not all-zero): LFC 0, stat 0, p 1, lfcSE kept.
9. **Cook's distances** (`calculateCooksDistance`, `results`):
   - dispersion from `robustMethodOfMomentsDisp` (not the final dispersions): trimmed variance per cell
     (`trimmedCellVariance`, trim 1/3, 1/4 or 1/8 by cell size, scale 2.04, 1.86 or 1.51), floored at 0.04; hat values
     from the fitted GLM;
   - cells are identical rows of the model matrix (`nOrMoreInCell`), not factor levels; maxCooks over samples in cells with
     ≥ 3 samples; a numeric covariate makes every sample its own cell, so maxCooks is NA and nothing is flagged;
   - cutoff qf(0.99, p, m − p) (port of the F quantile); a flagged gene gets **pvalue** NA (before filtering and BH);
   - with a single two-level factor, a flagged gene is unflagged when at least 3 counts exceed the outlier count.
   - Outlier **replacement and refit** (cells with ≥ 7 samples, `minReplicatesForReplace = 7`): replaced counts; the refit
     reuses the trend and dispPriorVar and re-estimates gene-wise and MAP dispersions only for replaced genes; baseMean from
     the replaced counts; maxCooks NA when every sample is replaceable; genes whose replaced counts are all zero get LFC 0,
     p 1. **Until phase 5 the port stops when any cell has ≥ 7 samples** (message names the R engine).
10. **Independent filtering** (`pvalueAdjustment`, alpha 0.1): θ = 50 points from mean(baseMean == 0) to 0.95 (to 1 if the
    start is ≥ 0.95); cutoffs by `quantile()` type 7; keep genes with baseMean ≥ cutoff; BH (`p.adjust`, NA-aware) at each
    θ; numRej = count of padj < alpha; `lowess(numRej ~ θ, f = 1/5)`; if max(numRej) ≤ 10, j = 1 (no filtering); else
    threshold = max(fit) − RMSE, where RMSE uses residuals only at θ with numRej > 0; j = first θ with numRej > threshold,
    else first above 0.9·max(fit), else 0.8·max(fit), else 1. Filtered genes get padj NA. Ports: R's `lowess`, quantile
    type 7, BH with NA handling.

Output: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj; and a diagnostics record (size factors, trend
coefficients, varLogDispEsts, dispPriorVar and which path, filter j and θ, counts and lists of dispersion outliers,
grid-search genes, `optim` rows, Cook's-flagged genes, filtered genes; notes).

## Where the code goes

- `crates/qrnamap-core/src/deseq/` (no I/O, `f64`): `checks.rs`, `size_factors.rs`, `model_matrix.rs`,
  `dispersion.rs` (gene-wise, trend, prior, MAP), `glm.rs` (IRLS, quasi-Newton fallback, Wald, contrasts), `cooks.rs`,
  `filter.rs` (independent filtering, BH, quantile type 7), `smooth.rs` (`lowess`, `loess` degree 2), `special.rs`
  (R-quality lgamma for large arguments, digamma, trigamma, normal tail to about 1e-300, `qf`, `dnbinom_mu`),
  `linalg.rs` (dense QR and Cholesky for p × p), `prior_table.rs` (the precomputed m − p ≤ 3 histograms). No new crate
  dependencies. Entry: `deseq(counts, model_matrix, contrast, opts) -> Result<DeseqResult, DeseqError>`.
- `qrnamap-wasm`: `deseq(request_json) -> response_json`; runs in the page worker.
- Page (`site/analyze/results.js`): free-text design formula (validated by `Design::parse`, default `~ <factor>`),
  contrasts on factors only, with a note for designs with interactions that the contrast is the effect at the other
  factors' reference levels; the volcano and table handle NA; the label says "DESeq2 method (Rust port of DESeq2 1.50.2
  defaults)" with the notes and divergences. `de-stats.js` is removed once the port is in.
- Desktop app: nothing extra.
- CLI: `qrnamap de --engine r|rust` (default `r`); `--engine rust` with `--shrink` is an error; the Rust engine writes
  `results.csv` and `normalized.csv`, no VST file, `versions.txt` with `engine: rust`.
- Licence: qRNAmap is GPL-3.0. DESeq2 is LGPL (≥ 3), R's `lowess`/`loess` code is GPL (≥ 2); all may be conveyed under
  GPL-3.0. Files that follow their source say so and credit it (Love, Huber, Anders 2014; R Core).

## Validation

1. **Reference dumper** (`validation/de/deseq2_reference.R`, extended; pixi environment `de`, on the HPCC): per gene
   `dispGeneEst`, `dispGeneIter`, `dispFit`, `dispersion`, `dispOutlier`, `dispMAP`, β, SE, `betaConv`, `betaIter`,
   `maxCooks`, the filter j and θ; per run size factors, trend coefficients, varLogDispEsts, dispPriorVar.
2. **Fixtures in the repo** (`crates/qrnamap-core/tests/data/deseq/`, a few hundred genes each, public data or R
   simulations with a fixed seed), so `cargo test` checks every stage without R, **natively and under wasm32**:
   - yeast amphotericin B (GSE80357), 3 vs 3; and 3 vs 3 plus a batch factor (m − p = 3);
   - 2 vs 2 (m − p = 2) and 1 vs 2 (m − p = 1); a group of 1 (1 vs 3); unbalanced 2 vs 5;
   - *A. fumigatus* `~ strain + condition`; an interaction design; a numeric covariate;
   - a 3-level factor with a contrast between two non-reference levels;
   - genes all zero in one group, and genes all zero in both contrasted levels;
   - rows that hit the dispersion grid and the GLM fallback; a dataset whose parametric trend fails (the port must stop);
   - a null dataset (max numRej ≤ 10, no filtering);
   - the two-group Cook's exception; a simulated set with ≥ 7 per group and planted outliers (phase 5).
3. **Full datasets on the HPCC**: the public ones above in full, plus *C. neoformans* H99 *lkh1* (GSE233612, 3 vs 3).
4. **Tolerances** (per stage, so a small bug cannot hide behind a loose end result):
   - size factors: relative difference < 1e-9;
   - gene-wise, trend and MAP dispersions: relative difference < 1e-6 for 99.9% of genes; the rest listed with both values;
   - β and log2 fold change on **every** non-all-zero gene: |Δ| < 1e-3·max(1, lfcSE); SE relative < 1e-4;
   - p-values: relative difference of −log10 p < 1e-4 where p > 1e-300;
   - exact match of the filter j, the Cook's-flagged set, and the NA pattern;
   - padj < 0.05 calls: every discordant gene listed, and each must have both padj within 1e-3 of 0.05;
   - documented divergences (GLM fallback rows) are listed, not counted.
   A dataset that misses a tolerance is a bug to explain, not a tolerance to loosen.
5. **Speed and memory in the browser** (laptop, a browser measurement): 6,000 and 12,000 genes × 6, 24 and 48 samples in
   Chromium, Firefox and WebKit; in `docs/browser/page-v1-measurements.md`.

## Phases

1. Reference dumper; fixtures; the m − p ≤ 3 histogram table (HPCC).
2. `special.rs`, `linalg.rs`, `smooth.rs` (`lowess` and `loess` against R on fixed inputs), quantile type 7, BH.
3. Checks, size factors, model matrix, IRLS (with its abandonment rules), gene-wise dispersions; per-stage fixture tests.
4. Trend (with the stop on failure), prior variance (both paths), MAP, final GLM with the quasi-Newton fallback, Wald and
   contrasts, Cook's flagging and the two-group exception, independent filtering; full results within tolerance on every
   fixture except the ≥ 7 replicate one (which must stop with the message).
5. Cook's replacement and refit; the ≥ 7 replicate fixture.
6. WebAssembly and the page (worker, formula box, contrast choosers, NA handling, notes); browser tests; speed.
7. `qrnamap de --engine rust`; full validation on the HPCC; results in `validation/de/README.md`. Then the owner decides
   whether the page drops the "screen" wording and whether the CLI default changes.

## Open questions

1. How often does the parametric trend fail on fungal data? Measure on the HPCC by running DESeq2 on PopSeq datasets and
   counting "local regression fit was automatically substituted". If often, port `locfit`'s method (a later design).
2. Exact R `lbfgsb`: port it later if the fallback rows turn out to change calls.

Settled: the m − p ≤ 3 prior is deterministic (`set.seed(2)`); counts are rounded once per gene half away from zero in
both the CLI and the page (`qrnamap-core` counts.rs), so both engines see the same integers.

## Review

An independent review (Opus 5.5, 2026-11-27) checked the first draft against the DESeq2 1.50.2 source. It found the
Cook's-distance step, the independent-filtering rule, the handling of non-converged GLM rows and the m − p ≤ 3 prior wrong
or incomplete, the IRLS and dispersion details imprecise (ridge scale, deviance convergence, μ floor, line-search rules),
the mean-fit fallback unsafe, the tolerances too loose to catch those errors, and many cases missing from the fixtures.
This revision takes all of those findings.

## Risks

- Fidelity of the line search, IRLS and the R smoothers; the per-stage fixtures find which step differs.
- Our special functions and `exp`/`ln` differ from R's in the last bits; the tight tolerances are relative, and the wasm32
  run of the fixtures catches platform differences.
- The filter's discrete jump: a difference of one rejection can change j; the exact-j tolerance makes that visible.
- Scope: shrinkage, transforms, the local fit and an exact `lbfgsb` are later designs.
