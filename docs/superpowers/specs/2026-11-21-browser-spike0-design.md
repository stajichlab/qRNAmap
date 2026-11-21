# Browser spike 0: design

Date: 2026-11-21. Status: approved by the user on 2026-11-21 with option F (fork). Parent design: `2026-11-19-qrnamap-design.md` (sections 9, 10, 15).

## Intent

The parent design says the browser build follows "spike 0" and only for tools that pass it. Spike 0 measures whether kallisto, HISAT2, featureCounts and DESeq2 run correctly and fit in a browser. This round does spike 0 only. It builds no page and no `qrnamap-wasm` crate. The output is a measured go or no-go per tool, written down.

Why: the page design (limits, genome mode, hosting) depends on these numbers. No number exists yet. The design records that `genomics-wasm` was never built or run by us.

## Goals

- For each of kallisto 0.52.0, HISAT2 2.2.2, featureCounts 2.1.1 and DESeq2: run the WASM build on yeast data and compare with the native tool from pixi.
- Measure index time and peak memory for a 12 Mb (yeast) and a 29 Mb (Af293) genome: HISAT2 index, kallisto index of the transcript set.
- Test cross-origin isolation, and browser download from NCBI and ENA, and Ensembl (CORS).
- Read the licences of `genomics-wasm` and of the bundled tools.
- Write all results in `docs/browser/spike0.md` with commands, versions, run IDs and artifact hashes.

## Non-goals

- No web page, no `qrnamap-wasm`, no worker code for the product.
- No Salmon in the browser (the CLI Salmon design excludes it).
- No change to defaults or to the CLI.
- No fix of upstream `genomics-wasm` code. A failure removes a tool from the browser build. It does not start a patch effort.

## Anti-goals

- A "pass" without the native comparison, or from a tool built from different sources than the pinned ones.
- Numbers from a different data set than the stated one, or estimated numbers written as measured.

## Facts checked on 2026-11-21

- Release `v0.1.0-rc.1` of `comparativechrono/genomics-wasm` (source commit 0b015f6) has HISAT2 2.2.2 (`core-packages.zip`), featureCounts 2.1.1 (`rna-analysis-packages.zip`) and a DESeq2 package (437 MB). It has no kallisto and no Salmon. The repository has build scripts for both (`tools/kallisto`, `tools/salmon`, workflow `quant.yml`).
- GitHub shows the repository licence as "Other". Its `LICENSE` file is not yet read.
- This cluster node has no Emscripten. The kallisto README requires Emscripten 3.1.64.
- `qRNAmap` is a private repository. Actions minutes there are billed to the account plan. UNVERIFIED: the plan's quota.

## Approaches

1. Release zips first, build only kallisto locally. Fastest. Not chosen.
2. **Build the tools in GitHub Actions (chosen by the user).** Clean Emscripten environment, same as upstream CI. Costs a CI loop per question. Browser tests also run in Actions, where Chromium and Firefox are available.
3. Build everything locally. Slowest. Not chosen.

## Design

- Pin the source: `comparativechrono/genomics-wasm` at commit 0b015f6 (the release commit). Use its own workflows for the builds (`quant.yml` for kallisto; the hisat2, featurecounts and deseq2 workflows). Do not rewrite the build steps.
- A new spike workflow takes the build artifacts, runs the native tools (pixi environment of this repository) and the WASM tools on the same inputs, and writes comparison tables and timings. Inputs: the yeast AMB subset used earlier (small, with checked accessions) and the two genomes.
- Node runs the parity tests. Playwright runs Chromium and Firefox for the browser tests: peak memory (Chromium CDP heap and process memory), isolation headers, CORS fetches. WebKit is best effort.
- Pass criteria (from the parent design, section 15): counts exactly equal for kallisto and featureCounts; HISAT2 SAM records equal apart from `@PG`; the run completes in a real browser without an out-of-memory error; the 12 Mb yeast index builds in 10 minutes or less and under 2 GiB. DESeq2 compares with stated numeric tolerances, not equality (parent design, section 15).
- Artifacts are checked by SHA-256 against the release manifest where the release has the same tool.

## Where CI runs (decided: option F, 2026-11-21)

- **Option F (recommended): fork `comparativechrono/genomics-wasm` under `stajichlab` and run the spike there.** The upstream workflows run unchanged. A public fork has free Actions minutes. The spike workflow is added on a branch of the fork. This creates a public repository under the lab organization.
- Option P: put the spike workflow in the private `qRNAmap` repository and check out upstream at the pinned commit. Keeps everything private. Uses billed minutes. UNVERIFIED how much the full DESeq2 (webR) build needs.

The user chose option F. The fork is created under `stajichlab`. Option P is not used.

## Left to the builder

Exact spike workflow layout, the choice of yeast subset, Playwright memory probes, and the table layout of `docs/browser/spike0.md`.

## Risks

- The DESeq2 webR build may be too large or slow for a free runner. Fallback: use the released `deseq2-package.zip` and record that it was not rebuilt.
- kallisto uses `MEMORY64=2` (experimental in Emscripten). The README says the output is ordinary 32-bit WASM. Spike 0a verifies this in practice.
- The service-worker headers approach for GitHub Pages is untested. Spike 0c tests it.
