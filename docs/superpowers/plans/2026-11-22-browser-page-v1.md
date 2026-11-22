# qRNAmap browser page, version 1: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A static page at `https://qrnamap.fungalgenomes.org/` where a visitor picks a genome, adds local FASTQ files or loads a test set, runs kallisto in the browser, and gets a mapping summary, a gene count table, a correlation heatmap, a PCA, and downloads.

**Architecture:** `qrnamap-core` gets the shared logic (gene matrix, QC numbers). A new `qrnamap-wasm` crate wraps `qrnamap-core`, `qseq-io` and `genome-catalog` with `wasm-bindgen`. A web worker runs the vendored kallisto WASM and `qrnamap-wasm`. `site/` holds the landing page (`/`) and the tool (`/analyze/`). The Pages workflow builds the WASM package and assembles `dist-site/`, as `../qmutfind` does.

**Tech Stack:** Rust (wasm32-unknown-unknown, wasm-bindgen), vendored kallisto 0.52.0 WASM, plain JavaScript modules, Playwright (Chromium, Firefox, WebKit) in GitHub Actions, pixi.

**Spec:** `docs/superpowers/specs/2026-11-22-browser-page-v1-design.md` (read it first). Evidence: `docs/browser/spike0.md`. Parent: `docs/superpowers/specs/2026-11-19-qrnamap-design.md`.

## Global Constraints

- Quantification mode only (kallisto). No DE, no alignment mode, no free ENA accession box. See the spec non-goals.
- `qseq-io`, `genome-catalog`, `qrnamap-core` stay pure (no I/O, no JSON dependency) and build for `wasm32-unknown-unknown`. JSON and `wasm-bindgen` live only in `qrnamap-wasm`.
- Acceptance rule: for each test set, with the same reads and preview size, the page's gene count table equals the CLI's `counts.tsv` (same bytes after zstd is removed). Gene summing, rounding and matrix building use one shared code path in `qrnamap-core`.
- Every table, plot and download made from a preview says "preview of the first N pairs". File names carry `preview`.
- Limits (N, hard limit, sizes) come from Task 3 measurements, not from guesses. State "UNVERIFIED" for anything not run. Browser claims are scoped: headless Chromium, Firefox, WebKit on a GitHub runner, stated datasets, stated dates.
- Vendored kallisto files: provenance file plus a CI SHA-256 check. Upstream notices ship beside them. Source pin: genomics-wasm commit `0b015f63e8048eb57b72ffd2706b79b803862b7c`, kallisto 0.52.0 revision `4e9f29cf3b021260415430c057a22469ca081391`, Emscripten 3.1.64. Our build: fork run 44939577353, artifact `quant-compiled-kallisto`.
- Copied `qmutfind` modules (`gz.js`, `reads.js`, `reference-panel.js`, `persist.js`, `theme.js`, plus `vendor/fflate.js` and its licence): record the `qmutfind` commit (`4404e1c93f47ea527287b93d6b496070563b193a` at the time of writing) in each file header. Local `qmutfind` is at `/bigdata/stajichlab/jstajich/projects/qmutfind` (read only).
- Reads never leave the browser. No secrets, keys or tokens in any file or workflow.
- Tools on the cluster: put `<worktree>/.pixi/envs/default/bin` on PATH; never `pixi run`; `CARGO_TARGET_DIR=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-page`; the new worktree has no pixi env, so Task 1 runs `pixi install` in the background with a log under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/`. `pixi.toml` uses `[workspace]`. No `BASH_SOURCE` in scripts that may run under SLURM. Jobs on the cluster use node-local `$SCRATCH` and default memory.
- Allowed outward actions: push branch `page-v1` to `origin`, and run workflows from it by dispatch or push (`gh workflow run`). Pushing or merging `main` is the user's decision. Do not delete branches or change repository settings. Browser tests may download public data from NCBI and ENA inside CI.
- Run `cargo fmt --all -- --check` and `cargo test --workspace` (with `QRNAMAP_REQUIRE_TOOLS=1`) before each Rust commit. The wasm check must pass: `cargo build -p qseq-io -p genome-catalog -p qrnamap-core --target wasm32-unknown-unknown`.
- Docs in Simplified Technical English. Measured numbers only. Commit trailer on every commit: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.

## Review Focus

- A preview shown as a full run: any table, plot or file without the "preview of the first N pairs" label or `preview` in its name.
- Bad input must give a clear message and no hang: FASTA names that do not match the GFF3; a GFF3 with no usable transcripts; a truncated or corrupt gzip; mates with different record counts; a FASTQ with a broken record; duplicate or unsafe sample names; an empty file.
- Small designs: one sample, two samples (no PCA, no correlation heatmap with a clear message), a sample with zero total counts.
- Interrupted runs: cancel during index or quant; reload during a run; saved results from another genome or another preview size must not be shown as current.
- Network failure for NCBI or ENA (offline, 404, 5xx): a message, not a blank page. Hard-limit and out-of-memory (worker crash) give a message that names the command-line route.
- A changed vendored kallisto file must fail CI (hash check).

## File Structure

- Create `crates/qrnamap-wasm/` (`Cargo.toml`, `src/lib.rs`), `crates/qrnamap-core/src/qc.rs`, `scripts/build-wasm-pkg.sh`, `scripts/build-site.sh`, `scripts/check-vendor.sh`.
- Create `site/` files: `index.html` (landing, replaces the placeholder), `analyze/index.html`, `analyze/app.js`, `analyze/worker.js`, `style.css`, copied modules, `testsets.json`, `vendor/kallisto/` (with `PROVENANCE.md`).
- Modify `Cargo.toml` (workspace member), `pixi.toml` (tasks, `wasm-bindgen-cli` pin equal to `Cargo.lock`), `.github/workflows/pages.yml`, `.github/workflows/ci.yml` (hash check), new `.github/workflows/page.yml` (browser tests), `README.md`.
- Create `docs/browser/page-v1-measurements.md`, `validation/browser-page/` (test-set expected tables, scripts).

---

### Task 1: Shared logic in core, and the `qrnamap-wasm` crate

**Files:** `crates/qrnamap-core/src/qc.rs`, `counts.rs` or a new module for the gene matrix, `lib.rs`, `crates/qrnamap-wasm/`, `Cargo.toml`, `pixi.toml`, `scripts/build-wasm-pkg.sh`; tests in `crates/qrnamap-core/tests/qc.rs` and `crates/qrnamap-wasm/tests/`.

**Interfaces — Produces (core):**
`pub enum QcError { TooFewSamples { need: usize, have: usize }, ZeroLibrary(String) }` (Display text names the sample).
`pub fn log_cpm(m: &CountMatrix) -> Result<Vec<Vec<f64>>, QcError>`: `log2(count / library_size * 1e6 + 1)`, genes by samples.
`pub fn sample_correlation(m: &CountMatrix) -> Result<Vec<Vec<f64>>, QcError>`: Pearson on `log_cpm`, over genes with total count above 0; needs 2 samples.
`pub struct Pca { pub scores: Vec<[f64; 2]>, pub variance_explained: [f64; 2] }` and `pub fn pca(m: &CountMatrix, top_genes: usize) -> Result<Pca, QcError>`: top genes by variance of `log_cpm` (page uses 500; fewer if the matrix has fewer), centred, needs 3 samples. (My call: this is a plain log-CPM PCA, not DESeq2 variance-stabilised; the page labels it so.)
**Existing to reuse:** `gene_counts_from_abundance`, `assemble`, `CountMatrix::to_tsv`, `build_transcripts`, `write_fasta`, `tx2gene_tsv`, `gene_table_tsv`, `genome_catalog::builtin`, `qseq_io::{read_fasta, parse_gff3}`, kallisto argv from `toolcall::kallisto_index/kallisto_quant`. Read how `crates/qrnamap-cli/src/pipeline.rs` builds the gene list and the matrix. If that step lives in the CLI crate, move it into `qrnamap-core` (`count_matrix(per_sample, model) -> CountMatrix`) and make the CLI call it, with the CLI output byte-identical (existing tests stay green).
**Interfaces — Produces (wasm crate, `wasm-bindgen`, names fixed, JSON strings or typed arrays at the builder's choice):** `catalog()` (genomes for the menu: accession, label, organism, size, fasta/gff URLs); `build_transcriptome(fasta_bytes, gff3_text) -> {fasta, tx2gene, genes_tsv}` or an error string; `gene_counts(abundance_tsv, tx2gene_tsv) -> json`; `assemble_matrix(samples_json) -> {tsv, ...}`; `qc(matrix_tsv) -> {correlation, pca} | error`; `kallisto_args(...)` that returns the same argv as `toolcall` for index and quant; `de_command(sheet_name, counts_name, design, contrast) -> string`.

- [ ] **Step 1: Write failing tests.** Core `qc.rs` tests with a hand-computed 3-gene, 3-sample matrix: `log_cpm` values; `sample_correlation` is 1 on the diagonal and symmetric; two identical samples give correlation 1; `pca` on a matrix with two clear groups separates them on PC1; errors: 1 sample (`TooFewSamples`), 2 samples for PCA, a zero-library sample (`ZeroLibrary` names it). Wasm crate test (Node, `wasm-bindgen-test` or the pkg-node script pattern of `../qmutfind/scripts/test-wasm.sh`): `build_transcriptome` on the yeast fixture FASTA and GFF3 from `validation/` or the test data gives the same `transcripts.fa` and `tx2gene` bytes as `qrnamap transcripts` on the same input; `gene_counts` plus `assemble_matrix` on the kallisto `abundance.tsv` files in `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures/native/` give the same `counts.tsv` bytes as the CLI path on the same files (use the `yeast.genes`/native fixture the Task 2 spike made; state which files you used).
- [ ] **Step 2: Run them. Expected: FAIL (items undefined).**
- [ ] **Step 3: Implement** the interfaces above. `qrnamap-wasm` depends on the pure crates and `wasm-bindgen`/`serde_json`. `scripts/build-wasm-pkg.sh <web|nodejs> <dir>` after `../qmutfind/scripts/build-wasm-pkg.sh` (arguments only, no `BASH_SOURCE`). Add pixi tasks `build-wasm-web` and `test-wasm`. Pin `wasm-bindgen-cli` to the `wasm-bindgen` version in `Cargo.lock`.
- [ ] **Step 4: Run** `cargo test --workspace`, fmt check, the wasm32 check, the pkg build and the Node wasm test. Expected: all pass. Commit (`feat: qrnamap-wasm crate and shared QC logic in core`).

---

### Task 2: Vendor kallisto and its runtime

**Files:** `site/vendor/kallisto/` (`kallisto.mjs`, `kallisto.wasm`, runtime loader files, notices, `PROVENANCE.md`), `scripts/check-vendor.sh`, `.github/workflows/ci.yml`.

**Produces:** a vendored directory the worker can load, and a hash check that fails on any change.

- [ ] **Step 1: Get the files** from the fork artifact `quant-compiled-kallisto` of run 44939577353 (`gh run download` into a new temp dir outside the worktree; treat as untrusted data, do not run it). Also read upstream `runtime/` at the pin to find the loader (`ToolRunner`) the page needs; vendor only what the worker needs, with its licence (MIT, per the repository `LICENSE`).
- [ ] **Step 2: Write `PROVENANCE.md`:** upstream repo and pin, kallisto revision, Emscripten version, our build run URL, artifact name, and a SHA-256 for every vendored file; the upstream notices (`THIRD_PARTY_LICENSES.txt`, kallisto and Bifrost licences) copied next to the files.
- [ ] **Step 3: Write `scripts/check-vendor.sh <dir>`** (absolute path argument): recompute SHA-256 values and compare with `PROVENANCE.md`; exit 1 on a mismatch, a missing file or an extra file. Test: it passes on the real directory and fails on a temp copy with one byte changed. Add it to `ci.yml`.
- [ ] **Step 4: Commit** (`feat(site): vendor kallisto 0.52.0 WASM with provenance`). Size must be under 3 MiB in total; report the sizes.

---

### Task 3: Measurements that fix the limits

**Files:** `docs/browser/page-v1-measurements.md`, `validation/browser-page/measure/` (scripts), `.github/workflows/page.yml` (measure job), small result TSVs.

**Consumes:** Task 2 files and Task 1 wasm package (build them in the workflow). **Produces:** the numbers the later tasks read: preview size `N`, the hard limit, the test-set read-fetch method.

- [ ] **Step 1: Build a measurement harness** in the style of `validation/browser-spike0/` (Node + Playwright, headless Chromium, Firefox, WebKit, one fresh browser per run, whole-browser PSS or the probe the spike used; reuse its scripts where they fit). Workflow `page.yml` with a `measure` job, dispatchable. Push branch `page-v1` (allowed) and dispatch.
- [ ] **Step 2: Measure kallisto quantification** (index of the yeast transcripts, then quant) at several read counts of the SRR3396381 subset (for example 100k, 200k, 500k, 1M, 2M pairs): time and peak memory, 3 runs each, per browser. Find the largest N that finishes in a stated time budget (propose 5 minutes) and under 1.5 GiB (propose; state the reason), and the count where it fails or passes 2 GiB. A failure is a result; record it.
- [ ] **Step 3: Truncated gzip.** From ENA, fetch the first part of `SRR3396381_1.fastq.gz` and `_2` by byte range, in the browser, and run them through kallisto. Record pass or fail and the exact error. Also record whether paired mates cut at different byte counts break pairing (they must be cut at matching record counts: describe the method that makes N pairs from two truncated files, for example read both streams to N records and stop). If truncated ENA reads cannot work, record that and choose the fallback in the spec (host small subset files as release assets of this repository; public data).
- [ ] **Step 4: The other open items:** read counts and file sizes of the yeast set's six samples (ENA portal API) and of the *C. neoformans* PRJNA976751 samples used for the test set; whether the *C. neoformans* paper's gene identifiers match the RefSeq H99 annotation (GCF_000149245.1) (the paper's Supplementary Table S3 may need manual download: if it cannot be fetched, say so and mark it UNVERIFIED); download size and transcript-building time for the three largest catalog genomes and for *S. cerevisiae*, *A. fumigatus*, *C. neoformans* H99.
- [ ] **Step 5: Write `docs/browser/page-v1-measurements.md`:** methods, versions, dates, all values, the chosen `N` (default) and hard limit, the test-set method, and every failed or UNVERIFIED item. Update the "To measure" section of the spec to point to it. Commit (`docs: measurements for the page limits`).

---

### Task 4: Site shell, landing page and build

**Files:** `site/index.html`, `site/analyze/index.html`, `site/style.css`, copied modules, `scripts/build-site.sh`, `.github/workflows/pages.yml`.

**Produces:** `/` with a link to the tool and the other items below; `/analyze/` with the three-step shell (Genome, Samples, Results; steps enabled in order); the copied modules with source headers; a build script and the Pages workflow.

- [ ] **Step 1: Landing page** (replaces the placeholder): what qRNAmap is; a clear "Analyze reads" button to `/analyze/`; the two test sets (names, organism, contrast); links to the CLI install steps, the repository, `docs/browser/spike0.md` and `docs/browser/page-v1-measurements.md`; a plain statement of what version 1 does not do (no DE, no alignment, preview for large inputs). Keep light and dark themes. A test (Node, HTML parse or Playwright) checks that the button link resolves in `dist-site`.
- [ ] **Step 2: Copy the `qmutfind` modules** named in Global Constraints. Add a header line to each: source path and commit. Keep their licence files. Change only what the page needs; list every change in `site/COPIED.md`.
- [ ] **Step 3: `scripts/build-site.sh`** after `../qmutfind/scripts/build-site.sh`: assemble `dist-site/`, stamp file names with a content hash (the `__BUILD__` scheme), write `CNAME` (`qrnamap.fungalgenomes.org`), copy the wasm package and `vendor/kallisto`, fail if a `__BUILD__` token is left. Update `pages.yml`: add the pixi and wasm build steps (as `../qmutfind` `pages.yml` does, same pinned action versions), keep the trigger paths, add `crates/**` and `scripts/**`.
- [ ] **Step 4: Run the build locally** (the wasm package from Task 1) and serve `dist-site` with a static server in a headless-browser smoke test (CI job in `page.yml`): both pages load in the three engines with no console errors. Commit (`feat(site): landing page, tool shell and site build`).

---

### Task 5: Genome step (menu, download, upload, transcriptome, index)

**Files:** `site/analyze/app.js`, `site/analyze/worker.js`, copied `reference-panel.js`, tests under `validation/browser-page/`.

**Consumes:** Task 1 wasm exports (`catalog`, `build_transcriptome`), Task 2 kallisto, Task 3 limits. **Produces:** a transcript set and a kallisto index in the worker, with progress, cancel and clear errors.

- [ ] **Step 1: Failing browser tests** (Playwright, run in `page.yml` on a local static server, small fixtures in `validation/browser-page/fixtures/` made from a tiny public genome or the yeast mitochondrial genome so they stay small): choose a catalog genome (use a small one; the test may stub the NCBI response with the real file bytes from fixtures) and reach "index ready"; upload FASTA + GFF3 and reach "index ready"; error cases from Review Focus: FASTA names not matching GFF3, GFF3 with no transcripts, truncated gzip, empty file, a genome above the size limit from Task 3 (clear message that names the command-line route); cancel during download and during index (state returns to "choose genome"); offline fetch error.
- [ ] **Step 2: Implement.** Menu from `catalog()` (searchable, NCBI RefSeq first); fetch with progress and byte counts; gunzip through `gz.js`; `build_transcriptome`; `kallisto index` in the worker with `-t 1` (the build requires it). Show genome label, accession, number of genes and transcripts. Remember the last genome in `persist.js` storage.
- [ ] **Step 3: Run the tests in the three engines in CI. Expected: pass.** Commit (`feat(site): genome step`).

---

### Task 6: Samples step and quantification

**Files:** `site/analyze/app.js`, `worker.js`, `copied reads.js`, tests.

**Consumes:** Task 5 index; Task 1 `gene_counts`, `assemble_matrix`, `kallisto_args`, sample-sheet validation. **Produces:** per-sample abundance, gene counts and mapping summary, for N pairs.

- [ ] **Step 1: Failing browser tests:** a sample table (name, condition column, mate 1 and mate 2 file pickers; single-end allowed with fragment length mean and SD fields, as the CLI requires); unsafe and duplicate sample names rejected with the CLI's reserved-name rules (read `RESERVED_OUTPUT_NAMES` in `pipeline.rs`); mates with different record counts give a clear error; the preview reads exactly N pairs from each mate (stream, stop at N, never read the whole file); a file above the hard limit is refused with the CLI route; a sample finishes and shows reads processed, pseudoaligned and the fraction (from `run_info.json`).
- [ ] **Step 2: Implement:** stream local files in chunks, decompress with `gz.js`, cut at N pairs (N from Task 3, user can lower it, never raise it above the hard limit), pass FASTQ to the vendored kallisto in the worker one sample at a time, free memory after each sample, per-sample progress and cancel. The page never holds a whole sample.
- [ ] **Step 3: Run the tests in the three engines. Expected: pass.** Commit (`feat(site): samples step and kallisto quantification`).

---

### Task 7: Results step

**Files:** `site/analyze/app.js`, `results.js` (new), `persist.js` use, tests.

**Consumes:** Task 6 per-sample results; Task 1 `qc`. **Produces:** the results view and downloads.

- [ ] **Step 1: Failing tests (browser, with the fixture samples):** mapping summary table; count table (sortable, search box, shows the preview label); correlation heatmap and PCA drawn as SVG; one or two samples show the clear message instead of the PCA/heatmap (Review Focus); downloads: `counts.preview-N.tsv` equals the bytes `assemble_matrix` returned, `samples.tsv` (the sheet format the CLI reads), and `de-command.txt` containing the exact `qrnamap de` command with the user's file names; reload restores saved results with the same genome and N, and does not restore results made with another genome or N.
- [ ] **Step 2: Implement** with plain SVG (no chart library). Every table and plot carries "preview of the first N pairs" (or "all pairs" when the input was below N). Label the PCA "PCA of log2(CPM+1), top 500 variable genes (not DESeq2-normalised)".
- [ ] **Step 3: Run the tests in the three engines. Expected: pass.** Commit (`feat(site): results view and downloads`).

---

### Task 8: Test sets

**Files:** `site/testsets.json`, `site/analyze/testsets.js`, `validation/browser-page/testsets/`, tests.

**Consumes:** Task 3 method for fetching reads; Task 5 and 6. **Produces:** two "Load test set" entries: yeast AmpB (GSE80357) and *C. neoformans* H99 Lkh1 (PRJNA976751, *lkh1* deletion against wild type, basal, 3 against 3).

- [ ] **Step 1: Write `testsets.json`** (name, organism, catalog accession, sample table with ENA run accessions, expected read counts from the ENA portal API, N, fetch method from Task 3, the contrast, the paper reference, the GEO series). Pick the 3+3 samples from the study metadata and record why.
- [ ] **Step 2: Failing tests:** the button fills the genome and the sample table; each sample is fetched (range or fallback files), cut to N pairs, and its record count equals N; ENA failure (offline, 404) gives a clear message; the test set run ends in the same results view.
- [ ] **Step 3: Implement** the fetch (browser `fetch` with `Range`, record count check). If Task 3 chose the fallback, create the subset FASTQ files from the public runs on the cluster, with ENA read-count checks, compress them, and publish them as assets of a GitHub release of this repository (outward action, allowed: release `page-testsets-v1` on this repository only; list the digests in `validation/browser-page/testsets/manifest.tsv`).
- [ ] **Step 4: Run the tests. Expected: pass.** Commit (`feat(site): test sets`).

---

### Task 9: Acceptance against the CLI, docs, and the final checks

**Files:** `validation/browser-page/` (CLI references, compare script), `docs/browser/page-v1.md`, `README.md`, `.github/workflows/page.yml`.

- [ ] **Step 1: CLI references.** For each test set, run `qrnamap quant`/`run` with the native kallisto for the same reads cut to the same N pairs (use the same preview reads: write the N-pair FASTQ from the browser's cut rule with a script so both sides get identical bytes; record SHA-256), and keep `counts.tsv` hashes in `validation/browser-page/expected/`. Cluster work under SLURM with default memory and `$SCRATCH`.
- [ ] **Step 2: Compare.** A Playwright job runs both test sets in Chromium, Firefox and WebKit and compares the page's count table bytes with the CLI reference. Pass means identical. A difference is reported with the first differing gene and blocks the task until explained. Record per set and engine: pass, fail, skipped.
- [ ] **Step 3: Check the *C. neoformans* result against the paper** in a stated way, without claiming a match the data do not show: count genes called by the CLI `de` command at adjusted p < 0.05 and |log2 fold change| of at least 2 on the CLI table, and report the overlap with the paper's list if Task 3 got Supplementary Table S3; otherwise say UNVERIFIED. This is a report, not a gate.
- [ ] **Step 4: Docs.** `docs/browser/page-v1.md`: what the page does, limits (N and the hard limit from Task 3), browsers tested, test-set results, known gaps. README: link to the page and the doc. Update the memory notes only if the user asks.
- [ ] **Step 5: Final checks:** `cargo fmt --all -- --check`, `cargo test --workspace` with `QRNAMAP_REQUIRE_TOOLS=1`, the wasm32 check, `scripts/check-vendor.sh`, the full `page.yml` run green on branch `page-v1`. Commit (`docs: browser page version 1 results`). Do not merge or push `main`; report the branch, the run IDs and the results table.
