# Browser spike 0

This directory holds the evidence for spike 0 of a future browser (WebAssembly) build of qRNAmap.
Spike 0 tests whether kallisto, HISAT2, featureCounts and DESeq2 builds from `comparativechrono/genomics-wasm` work.

## Fork and pin

- Fork: https://github.com/stajichlab/genomics-wasm (public fork of `comparativechrono/genomics-wasm`).
- Upstream pin: commit `0b015f63e8048eb57b72ffd2706b79b803862b7c` (release `v0.1.0-rc.1`).
- Work branch `qrnamap-spike0` in the fork was created at the pin. At creation, the fork `main` was also at the pin.
- All baseline runs used the unmodified upstream workflow files from the pin. No workflow edit, no patch.
- Allowed outward actions: create the fork, push `qrnamap-spike0` to the fork, run workflows in the fork. Nothing else.

## Files

- `baseline.tsv`: one row per tool with the run ID, conclusion, wall time, artifact names, source commit and notes.
- Release copies used for the hash comparison: `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0-recon/` (`core-packages.zip`, `rna-analysis-packages.zip`).
- DESeq2 fallback: `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/baseline/release/deseq2-package.zip`.
- Artifacts are not in git. They are in `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/baseline/`.

## Baseline result (all measured, 2026-11-21)

- kallisto (quant.yml): build job passed in 248 s. Salmon is out of scope. It ran in the same matrix because the workflow was not edited.
- HISAT2 (ci.yml): passed. The package zip hash equals the release asset.
- featureCounts (rna-analysis.yml): passed. The package zip hash differs from the release asset. Only release.json, SHA256SUMS and the validation reports differ. The WASM file is identical.
- DESeq2 (deseq2.yml): failed. An upstream source URL now returns HTTP 404. The run uploaded only `deseq2-validation-node`. DESeq2 is not rebuilt. The release `deseq2-package.zip` is the fallback.
- Browser validation results (Node, Chromium, Firefox, WebKit) are upstream data, not ours.

## Fixtures and native references (Task 2, measured 2026-11-21)

Data are in `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/fixtures/` (not in git). `manifest.tsv` here is a copy of
`fixtures/manifest.tsv` (file, bytes, SHA-256, how made). Check a fixture directory with
`check_manifest.sh /absolute/path/to/fixtures`. It exits 1 on a missing file or a mismatch.

- Reads: 200,000 pairs from SRR3396381 (AmpB_1, ENA read_count 8,405,534, both mates checked). Seed 20261121 (`subset_reads.py`).
- Genomes: S. cerevisiae R64 (GCF_000146045.2) and A. fumigatus Af293 (GCF_000002655.1) from NCBI through `qrnamap fetch`.
- Native tools (pixi env): kallisto 0.52.0, HISAT2 2.2.2, featureCounts 2.1.1, DESeq2 1.50.2 with R 4.5.3.
  The WASM comparison for DESeq2 uses R 4.6.0 upstream (CI), so the R versions differ.
- The Af293 HISAT2 run uses the yeast reads. Only 0.10% of reads align. It measures index build cost, not alignment quality.
- featureCounts uses a SAF file made by `qrnamap transcripts` (6,477 yeast genes), `-p --countReadPairs -s 0`.
- DESeq2 input: the six kallisto gene-count tables of the AMB run (6,477 genes), design `~ condition`, AmpB vs Control, no shrinkage.
- Times and peak RSS: `native/native_times.tsv` in the fixture directory. Run-to-run repeatability of HISAT2 SAM output was not tested.
- Scripts: `fixtures.sh` (stages `fetch`, `build`), `native_refs.sh`, `subset_reads.py`, `merge_counts.py`, `deseq2_ref.R`,
  `make_manifest.sh`, `how_made.tsv`, `check_manifest.sh`.
- Job requests (see `submit_fixtures_and_native.sh`): partition short, 2 CPUs, 8 GB, 2 h, for both jobs. Measured peak MaxRSS:
  111040K (fixtures job 32356703, 1m12s) and 1231516K (native job 32356735, 8m24s). The 8 GB request was not a measured need.
  The peak was 1.2 GB, so the default memory would likely have been enough. This was not tested.
- Exclude kallisto `run_info.json` from any comparison. It holds a start time and absolute paths.
- `fixtures.sh build` falls back to the hard-coded ENA read count (8,405,534) when ENA is unreachable. The record-count
  check in `subset_reads.py` is the check that enforces it.

## Parity runs (Task 3, measured 2026-11-21)

- Workflow: `workflow/qrnamap-spike0.yml`. `sync.sh SRC_DIR FORK_DIR "message" [--push]` copies it to
  `.github/workflows/` of the fork and copies `*.mjs`, `*.py`, `deseq2_ref.R` and `manifest.tsv` to `spike0/`.
  It pushes only branch `qrnamap-spike0`. It stops when the branch differs from the pin outside these two paths.
- Fork head: `955963bd1b67587e401b3c5c8acae3eb804bd01e` (two commits on the pin). Run 44944424373.
- Tools: the compiled artifacts of the Task 1 runs (kallisto 44939577353, HISAT2 44939585268,
  featureCounts 44939587731). Their `.wasm` and `.mjs` SHA-256 equal the local baseline copies.
  DESeq2: release `deseq2-package.zip` (DESeq2 1.52.0, R 4.6.0, webR), not rebuilt.
- Inputs: release `qrnamap-spike0-fixtures-v1` of the fork (reads, transcripts, SAF, native yeast SAM, count matrix)
  and the yeast genome from NCBI. Every file is checked against `manifest.tsv` in the workflow.
  Each WASM tool reads the same bytes as its native reference (`results/run-*/<engine>/staged.sha256`).
- featureCounts reads the native HISAT2 SAM. This tests featureCounts alone.
- `parity.mjs` writes `<tool>.steps.tsv` per tool and engine. `compare_parity.py` writes `results/parity_results.tsv`.
  Test: `python3 -m pytest validation/browser-spike0/test_compare_parity.py`.
- Result: kallisto, HISAT2 and featureCounts pass in Node, Chromium, Firefox and WebKit.
  DESeq2 (Node) fails the stated tolerance (absolute 1e-12 or relative 1e-6) on 1 of 38,203 numeric values.
  The table has 6,477 genes x 6 columns = 38,862 cells. 659 cells are NA on both sides (same NA pattern),
  so 38,203 numeric values were compared. DESeq2 in browsers was not run.
- DESeq2 versions: native DESeq2 1.50.2 with R 4.5.3, WASM DESeq2 1.52.0 with R 4.6.0. No same-version native
  reference was run. The failure cannot be attributed to WASM or to the version difference.
- kallisto index bytes differ between all four engines and from the native index. The quant output
  (`abundance.tsv`) is identical in all cases. The cause is UNVERIFIED.
- Any push to the fork branch without `[skip ci]` starts 11 upstream workflows (`on: push` with no filter).
  They were cancelled by hand after the trigger push of run 44944424373. The spike workflow is now registered
  (ID 426260115). Dispatch by workflow ID works: Task 4 run 44946369959 was started with
  `gh workflow run 426260115 --ref qrnamap-spike0 -f suite=index` (no push, no upstream runs).
- Full outputs: `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/parity/run-44944424373/` (SAM files zstd).

## Index size and time (Task 4, measured 2026-11-21)

### Commands in the packages

- HISAT2 package: `build` runs `hisat2-build-s` (README of `tools/hisat2` at the pin). So index building was measured,
  not index loading. kallisto package: `index` exists (README of `tools/kallisto`).
- Both built modules have `getHeapMax=()=>2147483648` in their `.mjs` (`results/index-run-*/<engine>/heap-max.txt`).
  The kallisto README states "The WebAssembly heap can grow to 2 GiB". This 2 GiB is a stated limit. It is not a measurement.

### Method

- Workflow job `index-size` in `workflow/qrnamap-spike0.yml` (matrix node, chromium, firefox, webkit).
  Input `suite` (`parity`, `index`, `all`) selects jobs on dispatch. After this change a push runs only the index job.
  The Task 3 parity jobs (`parity`, `deseq2-node`) run only on dispatch with `suite=parity` or `suite=all`.
- Run 44946369959 (dispatch by workflow ID 426260115, ref `qrnamap-spike0`), fork head
  `4333f4f56d5a57cd69817aaefe87c21f2e9c15bb`. All four jobs passed. One run in total.
- Runner: ubuntu-24.04 GitHub runner, 16,367,524 kB RAM, 3,145,724 kB swap (`free-before.txt`). Node v22.23.5.
  Chromium 156.0.8078.4, Firefox 157.0, WebKit 27.2. `crossOriginIsolated` was true in all three browsers.
- Tool files: same `.wasm` and `.mjs` SHA-256 as Task 3. Upstream runtime `ToolRunner` (one worker per run).
- Inputs (staged SHA-256 equal to local decompression of the fixtures):
  - `yeast.fa` 12,310,392 bytes, 12,157,105 bases, 17 sequences (gunzip of the NCBI R64 `.fna.gz`).
  - `afum.fa` 29,752,955 bytes, 29,384,958 bases, 8 sequences (gunzip of the NCBI Af293 `.fna.gz`).
  - `yeast.transcripts.fa` 6,478 transcripts. `afum.transcripts.fa` 10,086 transcripts (zstd -d of the fixtures).
  - `afum.transcripts.fa.zst` was added to the fixture release `qrnamap-spike0-fixtures-v1` (GitHub digest equals `manifest.tsv`).
- Arguments equal the native arguments: `build -p 1 X.fa X.hs`; `index -i X.idx X.transcripts.fa`.
- `index_size.mjs`: 4 cases x 3 runs per engine. Each run uses a fresh process: a child Node process, or a new
  browser launch. Per-run timeout 45 min (not reached; the longest run was 84.41 s).
- `index_size_wrap.mjs` is given to `ToolRunner` as the module URL. It calls the real tool factory with the same
  options and wraps `callMain`. It reads `Module.HEAPU8.length` before and after `callMain`. The tool code is not changed.
- `summarize_index_size.py` writes `results/index_size_runs.tsv` (48 runs) and `results/index_size_summary.tsv`.

### Memory probes and their limits

| probe | engines | method | limit |
|---|---|---|---|
| WASM linear memory | all | `Module.HEAPU8.length` after `callMain`, in the tool worker (`wasm_heap_final_bytes`) | This is the final size. It equals the peak only if linear memory never shrinks during the run; WebAssembly memory can only grow, but this was not measured here. It excludes MEMFS files (JS arrays), input copies and engine overhead. |
| `node_maxrss_kb` | node | `process.resourceUsage().maxRSS` of the child process | Kernel peak RSS of the whole process (main thread, worker, WASM, input and output copies). |
| `node_rss_peak_kb` | node | `process.memoryUsage().rss` every 250 ms | Can miss a peak shorter than 250 ms. |
| `tree_pss_*_kb` | browsers | sum of `Pss` in `/proc/PID/smaps_rollup` over all descendants of the harness process, every 250 ms | Whole browser (all processes), not only the tool. Can miss a peak shorter than 250 ms. Baseline = after page load and input fetch, before the run. Pss was readable for all processes in all runs. |
| `cdp_js_heap_*` | chromium | CDP `Performance.getMetrics` `JSHeapUsedSize`, every 1 s | Page main thread only. It does not include the worker or WASM memory. |

- All browsers ran headless: `index_size.mjs` and `parity.mjs` call `launch({ headless: true })`.
- Firefox and WebKit have no CDP in Playwright. For them, the process probe (`/proc` Pss) and the WASM probe are the only methods.
- `performance.measureUserAgentSpecificMemory()` was not used.

### Results (3 runs each; run_s = `ToolRunner.run`, includes worker start and file copies)

| case | engine | run_s | WASM heap MiB | peak memory MiB (probe) | outcome |
|---|---|---|---|---|---|
| HISAT2 yeast 12 Mb | node | 5.36, 5.39, 5.33 | 113.3 x3 | 298.9, 300.4, 305.5 (maxRSS) | finished x3 |
| HISAT2 yeast 12 Mb | chromium | 6.95, 6.98, 7.00 | 113.3 x3 | 449.7, 450.3, 447.0 (tree Pss) | finished x3 |
| HISAT2 yeast 12 Mb | firefox | 7.38, 7.38, 7.36 | 113.3 x3 | 670.4, 669.8, 671.7 (tree Pss) | finished x3 |
| HISAT2 yeast 12 Mb | webkit | 7.46, 7.49, 7.43 | 113.3 x3 | 675.6, 667.1, 667.0 (tree Pss) | finished x3 |
| HISAT2 Af293 29 Mb | node | 11.71, 11.99, 12.14 | 134.7 x3 | 412.3, 424.7, 413.5 | finished x3 |
| HISAT2 Af293 29 Mb | chromium | 14.60, 14.65, 14.66 | 134.7 x3 | 531.4, 531.1, 533.9 | finished x3 |
| HISAT2 Af293 29 Mb | firefox | 15.29, 15.40, 16.21 | 134.7 x3 | 754.2, 717.3, 718.8 | finished x3 |
| HISAT2 Af293 29 Mb | webkit | 15.97, 15.88, 15.79 | 134.7 x3 | 834.9, 835.3, 834.2 | finished x3 |
| kallisto yeast | node | 19.44, 19.34, 19.58 | 117.1 x3 | 289.8, 280.3, 274.7 | finished x3 |
| kallisto yeast | chromium | 20.95, 20.78, 21.14 | 117.1 x3 | 369.9, 371.2, 371.7 | finished x3 |
| kallisto yeast | firefox | 36.80, 37.03, 36.80 | 117.1 x3 | 650.0, 653.5, 649.7 | finished x3 |
| kallisto yeast | webkit | 21.77, 21.38, 22.28 | 117.1 x3 | 824.4, 823.9, 828.8 | finished x3 |
| kallisto Af293 | node | 43.96, 44.18, 44.12 | 204.4 x3 | 411.0, 412.0, 415.4 | finished x3 |
| kallisto Af293 | chromium | 46.99, 47.08, 47.15 | 204.4 x3 | 483.0, 491.9, 481.1 | finished x3 |
| kallisto Af293 | firefox | 84.41, 84.37, 83.89 | 204.4 x3 | 781.1, 780.2, 770.2 | finished x3 |
| kallisto Af293 | webkit | 49.52, 49.79, 49.42 | 204.4 x3 | 1170.4, 1176.5, 1170.1 | finished x3 |

- Browser tree Pss baseline (before the run, browser idle with the input loaded): chromium 187 to 198 MiB,
  firefox 442 to 521 MiB, webkit 364 to 412 MiB. All values: `results/index_size_summary.tsv`.
- Chromium CDP JS heap used (page main thread): 0.7 to 0.8 MiB peak in all runs.
- Native reference (Task 2, cluster, `/usr/bin/time -v`, one run): HISAT2 yeast 8.83 s, 107.0 MiB peak RSS;
  HISAT2 Af293 24.63 s, 128.2 MiB; kallisto yeast 21.7 s, 136.5 MiB. kallisto Af293: not measured natively.
  Native and WASM ran on different machines, so the times are not a direct speed comparison.
- HISAT2 index files: identical SHA-256 in all 12 runs per genome (4 engines x 3). Total bytes equal the native
  index (yeast 22,818,973; Af293 48,273,818). Native index bytes were not kept, so no byte comparison with native.
- kallisto index files: all 12 runs per transcript set give different bytes and slightly different sizes, also
  between runs in the same engine (yeast 5,831,519 to 5,832,013 bytes; Af293 10,611,076 to 10,611,476 bytes).
  The cause is not known. Task 3 showed identical quant output from such indexes.

### Criteria

- 12 Mb yeast index in 600 s or less and under 2 GiB: pass in node, chromium, firefox and webkit (all 3 runs).
  Browser tree Pss peaks include the idle browser baseline (all browser processes), not only the tool.
  The largest yeast value on any probe is 675.6 MiB (webkit tree Pss); the longest run is 7.49 s.
- The other cases are within 600 s and 2 GiB on every probe in every engine (information, not a criterion).
  The largest value is 1176.5 MiB (kallisto Af293, webkit tree Pss, whole browser).
- No run failed, ran out of memory, crashed or timed out. No tool is "CLI only" on the basis of this test.
- Failures: none. Skipped: none. Not tested: genomes larger than 29 Mb, multi-threaded builds, devices with less
  memory than the 16 GB runner, browsers other than the Playwright builds.
