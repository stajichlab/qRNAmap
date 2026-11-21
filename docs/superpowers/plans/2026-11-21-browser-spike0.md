# Browser spike 0: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure whether kallisto, HISAT2, featureCounts and DESeq2 run correctly and fit in a browser, using GitHub Actions in a fork of `genomics-wasm`, and write a go or no-go per tool.

**Architecture:** A fork `stajichlab/genomics-wasm` pinned at commit 0b015f6 builds the tools with upstream workflows. Our spike scripts live in this repository under `validation/browser-spike0/`. A sync script copies them to branch `qrnamap-spike0` of the fork, where a spike workflow runs parity, size, isolation and CORS tests. Results come back with `gh` into `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/` and a summary goes in `docs/browser/spike0.md`.

**Tech Stack:** GitHub Actions (ubuntu-24.04), Emscripten 3.1.64 (via upstream scripts), Node 22, Playwright (Chromium, Firefox, WebKit), pixi native tools on the cluster, `gh`.

**Spec:** `docs/superpowers/specs/2026-11-21-browser-spike0-design.md` (read it first). Parent: `docs/superpowers/specs/2026-11-19-qrnamap-design.md` sections 9, 10, 15.

## Global Constraints

- Only these outward actions are allowed: create the fork `stajichlab/genomics-wasm`, push branch `qrnamap-spike0` to it, and run workflows in it. No push to any other remote, no change to `comparativechrono/genomics-wasm`. If the fork cannot be created (permissions), stop and report BLOCKED.
- Pin upstream at commit `0b015f63e8048eb57b72ffd2706b79b803862b7c`. Do not rewrite upstream build steps. No patches to upstream tools.
- Measured numbers only. Write UNVERIFIED for anything not run. A failed tool is reported as failed; it is not tuned until it passes.
- Same inputs and same tool versions on both sides of every comparison: kallisto 0.52.0, HISAT2 2.2.2 (note any native difference), featureCounts 2.1.1. DESeq2: upstream native uses R 4.6.0 (CI) and our CLI uses 1.50.2; state the versions on both sides of each DE comparison, with numeric tolerances (not equality).
- Pass criteria (parent design section 15): kallisto and featureCounts counts exactly equal; HISAT2 SAM records equal apart from `@PG`; the run completes in a real browser without an out-of-memory error; the 12 Mb yeast index builds in 10 minutes or less and under 2 GiB.
- No secrets in any file or workflow. Large data stays out of git (manifest with SHA-256 only). Compress kept outputs with zstd.
- Native reference jobs on the cluster: SLURM, default `--mem-gb 8` unless measured otherwise, node-local `$SCRATCH` for temp I/O, never `BASH_SOURCE` in SLURM scripts, run tools with the pixi env bin on PATH (never `pixi run`), `pixi.toml` uses `[workspace]`.
- Docs: Simplified Technical English. Commit trailers on every commit: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01CoSBcPmuXVDpBrrPNsnYtc`.

## Review Focus

- The fork stays at the pinned commit (record `git rev-parse` of the fork base and of every artifact's source commit).
- A WASM tool and its native reference read byte-identical inputs (SHA-256 recorded for both).
- Browser memory is measured by a stated method per browser, not copied from Node numbers.
- Isolation and CORS results are per URL and per browser, with the exact response headers seen.
- A failed or skipped test appears in the write-up as failed or skipped, with its log reference.

## File Structure

- Create `validation/browser-spike0/` (scripts, manifest, sync script, README). Create `docs/browser/spike0.md`.
- In the fork, branch `qrnamap-spike0`: `.github/workflows/qrnamap-spike0.yml` and `spike0/` copied by the sync script.
- Outputs under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/` (not in git).

---

### Task 1: Fork, pin and baseline the upstream builds

**Files:** `validation/browser-spike0/README.md`, `validation/browser-spike0/baseline.tsv`.

**Produces:** a fork with branch `qrnamap-spike0` at the pinned commit; upstream workflow run IDs and artifacts for kallisto, HISAT2, featureCounts, DESeq2.

- [ ] **Step 1: Check `gh auth status` has rights to create a repository under `stajichlab`.** If not, stop with BLOCKED and the exact error.
- [ ] **Step 2: Create the fork** (`gh repo fork comparativechrono/genomics-wasm --org stajichlab --clone=false`), create branch `qrnamap-spike0` at the pinned commit, enable Actions on the fork.
- [ ] **Step 3: Run the upstream workflows** `quant.yml` (kallisto only: restrict the matrix by workflow input or by dispatching and cancelling the salmon job), `ci.yml` (HISAT2), `rna-analysis.yml` (featureCounts) and `deseq2.yml` on the branch. Wait with bounded checks, not tight polling.
- [ ] **Step 4: Record in `baseline.tsv`:** tool, workflow, run ID, conclusion, wall time, artifact names, source commit. For HISAT2 and featureCounts compare the SHA-256 of the built package with the release asset and report equal or different. Report upstream's own validation results (Node and browsers) as "upstream data, not ours".
- [ ] **Step 5: If a build fails,** record the failing step and log excerpt and continue with the other tools. DESeq2 fallback if its build fails or exceeds the runner limits: use the release `deseq2-package.zip` and write "not rebuilt".
- [ ] **Step 6: Commit** (`docs: spike 0 baseline of the upstream builds`).

---

### Task 2: Fixtures and native references

**Files:** `validation/browser-spike0/fixtures.sh`, `native_refs.sh`, `manifest.tsv`; data under `qRNAmap_data/spike0/fixtures/`.

**Produces:** a small yeast read subset, the two genomes with annotation and transcript sets, and native outputs for every tool.

- [ ] **Step 1: Read subset.** From the AMB reads already in `qRNAmap_data/amb`, take a deterministic subset (fixed seed, stated number of pairs, about 200,000) for one sample; gzip it. Record SHA-256 and the source accessions with their ENA read counts.
- [ ] **Step 2: Genomes.** Download S. cerevisiae R64 (about 12 Mb) and A. fumigatus Af293 (about 29 Mb) with GFF3 from NCBI using the URLs from `genome-catalog`; record sizes and SHA-256. Build the transcript FASTA with `qrnamap transcripts` (the CLI on `main`).
- [ ] **Step 3: Native references** (SLURM job, pixi env on PATH): `kallisto index` and `kallisto quant` (0.52.0, `-t 1`, `--plaintext`); `hisat2-build` and `hisat2` for both genomes (time and `/usr/bin/time -v` peak RSS); `featureCounts` on the HISAT2 output; DESeq2 1.50.2 on a small count matrix from the AMB run. Keep outputs zstd-compressed.
- [ ] **Step 4: Write `manifest.tsv`** (file, bytes, SHA-256, how made). Test: a `check_manifest.sh` that recomputes the hashes and exits non-zero on a mismatch. Run it.
- [ ] **Step 5: Commit** scripts and manifest (`feat(validation): browser spike 0 fixtures and native references`).

---

### Task 3: Parity in Node and in browsers (spike 0a)

**Files:** `validation/browser-spike0/workflow/qrnamap-spike0.yml`, `parity.mjs`, `compare_parity.py`, `sync.sh`.

**Consumes:** Task 1 artifacts and fork; Task 2 fixtures (read subset uploaded as a release asset of the fork; genomes downloaded in the workflow).

**Produces:** parity tables per tool and engine.

- [ ] **Step 1: Write `sync.sh`** (absolute paths from arguments, no `BASH_SOURCE`) that copies `validation/browser-spike0/{workflow,*.mjs,*.py}` to the fork branch and pushes only branch `qrnamap-spike0`.
- [ ] **Step 2: Write `parity.mjs`** that runs the WASM kallisto, HISAT2 and featureCounts on the fixture (use the upstream `validate*.mjs` scripts as templates for loading the packages) in Node and in Playwright engines `chromium`, `firefox`, `webkit`; and DESeq2 in Node. Output a TSV per tool and engine.
- [ ] **Step 3: Write `compare_parity.py`** with checks: kallisto `abundance.tsv` `est_counts` and `tpm` columns identical; featureCounts count table identical; HISAT2 SAM records identical after removing `@PG`; DESeq2 results within stated tolerances (report max absolute and relative differences). A pytest-style test with a tiny hand-made pair (equal and unequal) shows the comparer fails on a difference.
- [ ] **Step 4: Run the workflow** in the fork, download the artifacts, run `compare_parity.py` against the Task 2 native references. Record per tool and engine: pass, fail, skipped, with the first difference.
- [ ] **Step 5: Commit** (`feat(validation): browser spike 0 parity runs`) with small result TSVs.

---

### Task 4: Index size and time (spike 0b)

**Files:** `validation/browser-spike0/index_size.mjs`, results TSVs.

- [ ] **Step 1: Check first** which commands each WASM package supports (for example whether `hisat2-build` exists in the package). If HISAT2 index building is not in the package, record that and test index loading instead; do not build one by other means.
- [ ] **Step 2: Measure** HISAT2 index (12 Mb and 29 Mb) and kallisto index (Af293 transcript set): wall time and peak memory in Node (process RSS and WASM heap) and in Chromium (CDP `Performance.getMetrics` and process memory) and Firefox (state the method and its limit). Run each three times; report all three values.
- [ ] **Step 3: Apply the criteria** (12 Mb index: 10 minutes or less and under 2 GiB). A tool that fails is "CLI only" in the write-up. Record the 2 GiB WASM heap limit from the kallisto README as the stated limit, not a measurement.
- [ ] **Step 4: Commit** (`feat(validation): browser spike 0 index size and time`).

---

### Task 5: Isolation and cross-origin downloads (spike 0c, 0d)

**Files:** `validation/browser-spike0/isolation.mjs`, results TSVs.

- [ ] **Step 1: Test** with a local server in two modes (no headers; `Cross-Origin-Opener-Policy: same-origin` plus `Cross-Origin-Embedder-Policy: require-corp`) and the service-worker approach of upstream, in Chromium, Firefox and WebKit: `crossOriginIsolated` value; featureCounts with threads; kallisto and HISAT2 without isolation.
- [ ] **Step 2: Fetch from the page** in each mode and browser: an NCBI genome file, an NCBI GFF3 file, an ENA FASTQ (full and a byte range), and an Ensembl file. Record per URL: status, `Access-Control-Allow-Origin`, `Cross-Origin-Resource-Policy`, and whether the body could be read. Record the exact URLs and date.
- [ ] **Step 3: Commit** (`feat(validation): browser spike 0 isolation and CORS`).

---

### Task 6: Licence check (spike 0e)

**Files:** section in `docs/browser/spike0.md`.

- [ ] **Step 1: Read** `LICENSE`, `THIRD_PARTY.md` and each tool's `THIRD_PARTY_LICENSES.txt` and `README.md` in the fork at the pinned commit. List per tool: tool licence, bundled library licences, and the repository licence as written.
- [ ] **Step 2: Write facts only.** State what the texts say about redistribution of built bundles (for example GPL tools in a static page). Give no legal conclusion; list it as the user's decision.

---

### Task 7: Write-up and decision table

**Files:** `docs/browser/spike0.md`, `README.md` (status line), memory note.

- [ ] **Step 1: Write `docs/browser/spike0.md`:** methods, versions, run IDs, a table with one row per tool and columns parity (Node, Chromium, Firefox, WebKit), index time and memory, isolation need, and go or no-go with the reason. Include measured numbers only; list every failed or skipped test.
- [ ] **Step 2: State plainly** which tools can enter the browser build, what must stay CLI-only, and what remains UNVERIFIED. No recommendation beyond the table; the next step (the first page) is the user's decision.
- [ ] **Step 3: Run** `cargo fmt --all -- --check` and the workspace tests (no code changes expected; confirm). Commit (`docs: browser spike 0 results`).
