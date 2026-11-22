# Browser spike 0 results

Date of all measurements: 2026-11-21. Status: spike 0 is complete. The next step (the first page) is the user's decision.
Evidence: `validation/browser-spike0/` (README, `baseline.tsv`, `results/*.tsv`). Raw artifacts are not in git. They are in
`/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/`.

## Summary

- Spike 0 measured the WebAssembly (WASM) builds of kallisto 0.52.0, HISAT2 2.2.2, featureCounts 2.1.1 and DESeq2 1.52.0
  from `genomics-wasm`. It compared their output with the native tools on the same input bytes.
- Where: GitHub Actions in the public fork https://github.com/stajichlab/genomics-wasm, branch `qrnamap-spike0`.
  Runner: GitHub-hosted `ubuntu-24.04`, 16,367,524 kB RAM. Node v22.23.5. Playwright browsers, all headless:
  Chromium 156.0.8078.4, Firefox 157.0, WebKit 27.2.
- Upstream pin: commit `0b015f63e8048eb57b72ffd2706b79b803862b7c` (release `v0.1.0-rc.1`). Our additions on the fork branch
  end at head `7cba1e97b8d25175789a8ed83c042a3485284648`. They touch only `.github/workflows/qrnamap-spike0.yml` and `spike0/`.
- Fixtures: 200,000 read pairs of SRR3396381 (AmpB_1, seed 20261121), S. cerevisiae R64 (GCF_000146045.2, 12,157,105 bases)
  and A. fumigatus Af293 (GCF_000002655.1, 29,384,958 bases) from NCBI.
- Result in one line per tool (details in the decision table):
  - kallisto: all measured criteria pass in Node, Chromium, Firefox and WebKit. Index bytes are not reproducible (cause unknown).
  - HISAT2: all measured criteria pass in Node, Chromium, Firefox and WebKit.
  - featureCounts: all measured criteria pass in Node, Chromium, Firefox and WebKit. It needs cross-origin isolation.
  - DESeq2: 1 of 38,203 values is outside the stated tolerance in Node, with a version difference. Browsers were not run.
- Scope of every "pass" in this file: headless Playwright browsers, one 16 GB GitHub runner, one read subset, one date.
  It is not a test of user devices, headed browsers, mobile browsers or full-size samples.

## Methods and versions

### Tool builds

| tool | WASM build | build run | native tool (pixi, cluster) |
|---|---|---|---|
| kallisto | 0.52.0, built by upstream `quant.yml` | 44939577353 (success) | kallisto 0.52.0 |
| HISAT2 | 2.2.2, built by upstream `ci.yml`; zip SHA-256 equals the release asset | 44939585268 (success) | HISAT2 2.2.2 |
| featureCounts | 2.1.1 (Subread), built by upstream `rna-analysis.yml`; `.wasm` equals the release, zip metadata differs | 44939587731 (success) | featureCounts 2.1.1 |
| DESeq2 | 1.52.0 with R 4.6.0 (webR). Release `deseq2-package.zip`, not rebuilt, not validated by us | 44939590596 (failure) | DESeq2 1.50.2 with R 4.5.3 |

- The upstream DESeq2 build failed at step "Assemble pinned R and DESeq2 browser files". The DESeq2 1.52.0 source URL on
  `bioc-release.r-universe.dev` redirects to a URL that returns HTTP 404 (checked 2026-11-21). So the spike used the release
  zip (SHA-256 `ed28cb91a15a18629f627708f76a18b90b4ed69089e7f298575be634ae55ad28`, equal to the release digest).
- DESeq2 versions differ between the two sides: native 1.50.2 with R 4.5.3, WASM 1.52.0 with R 4.6.0.
  No same-version native reference was run.
- All builds used the unmodified upstream workflow files at the pin. No patch to any tool.

### Spike runs (fork workflow `qrnamap-spike0.yml`)

| part | run | fork head | what it measured |
|---|---|---|---|
| 0a parity | 44944424373 | `955963bd1b67587e401b3c5c8acae3eb804bd01e` | WASM output against native output, Node and 3 browsers; DESeq2 in Node only |
| 0b index size | 44946369959 | `4333f4f56d5a57cd69817aaefe87c21f2e9c15bb` | index time and memory, 4 cases x 3 runs x 4 engines (48 runs) |
| 0c, 0d isolation and CORS | 44963830025 | `7cba1e97b8d25175789a8ed83c042a3485284648` | isolation modes, tools with and without isolation, 13 cross-origin requests |

- Each part ran once. No reruns.
- Each WASM tool read the same bytes as its native reference (SHA-256 recorded in `results/run-*/<engine>/staged.sha256`).
- Parity comparisons (`compare_parity.py`): kallisto `abundance.tsv` `est_counts` and `tpm` identical; HISAT2 SAM records and
  header identical apart from `@PG`; featureCounts count table and summary identical apart from the first comment line
  (program path); DESeq2 per value absolute difference 1e-12 or less, or relative difference 1e-6 or less. The DESeq2 tolerance
  was written before the run and was not changed after it.
- featureCounts parity read the native HISAT2 SAM. It tests featureCounts alone.
- The featureCounts output SHA-256 in Node (`5186a514...`) differs from the browser value (`8e375993...`) because the first comment line holds the program path (worker path in Node, `./this.program` in browsers). Parity ignores that line.
- The parity browser runs used `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`
  (`crossOriginIsolated` true). Part 0c tested the tools without these headers.
- Native references ran on the UCR cluster (SLURM jobs 32356703 and 32356735), with one thread. Native and WASM times come
  from different machines. They are not a speed comparison.

### Memory probes

| engine | probe | what it includes |
|---|---|---|
| Node | `process.resourceUsage().maxRSS` of a fresh child process | whole process peak RSS |
| Chromium, Firefox, WebKit | sum of `Pss` from `/proc/PID/smaps_rollup` over the browser process tree, every 250 ms | whole browser, idle baseline included (Chromium 187 to 198 MiB, Firefox 442 to 521 MiB, WebKit 364 to 412 MiB); can miss a peak shorter than 250 ms |
| all | `Module.HEAPU8.length` after `callMain` | final WASM linear memory only |
| Chromium | CDP `JSHeapUsedSize`, every 1 s | page main thread only (0.7 to 0.8 MiB); not the worker or WASM memory |

- The 2 GiB WASM heap limit (`getHeapMax` 2147483648 in both `.mjs` files; kallisto README) is a stated limit, not a measurement.

## Decision table

Criteria (parent design `2026-11-19-qrnamap-design.md`, section 11b; the spike plan calls it section 15; the parent design
has no section 15): (1) counts exactly equal for kallisto and featureCounts; HISAT2 SAM records equal apart from `@PG`;
(2) the run completes in a real browser without an out-of-memory error; (3) the 12 Mb yeast index builds in 10 minutes or
less and under 2 GiB. DESeq2 uses stated numeric tolerances, not equality. A failure removes a tool from the browser build,
not from the CLI.

Index columns: range over 3 runs in each of Node, Chromium, Firefox and WebKit. Memory: Node maxRSS, browser tree Pss (whole browser).
Per-engine values are in the next table.

| tool | parity Node / Chromium / Firefox / WebKit | index time, s | index peak memory, MiB | needs isolation headers (measured) | result against criteria | verdict |
|---|---|---|---|---|---|---|
| kallisto 0.52.0 | pass / pass / pass / pass (`abundance.tsv` byte-identical, 6,478 targets) | yeast transcripts 19.34 to 37.03; Af293 transcripts 43.96 to 84.41 | yeast 274.7 to 828.8; Af293 411.0 to 1,176.5; WASM heap 117.1 and 204.4 | no (ran in all 3 browsers without headers) | (1) pass, (2) pass, (3) not the criterion case (genome index); transcript indexes within 600 s and 2 GiB | `go with a stated caveat` (criteria 1-3, 200,000-pair subset, headless runner): index bytes differ between every run, cause unknown; quant output identical |
| HISAT2 2.2.2 | pass / pass / pass / pass (430,922 records, header without `@PG` identical) | yeast 12 Mb 5.33 to 7.49; Af293 29 Mb 11.71 to 16.21 | yeast 298.9 to 675.6; Af293 412.3 to 835.3; WASM heap 113.3 and 134.7 | no (ran in all 3 browsers without headers) | (1) pass, (2) pass, (3) pass: max 7.49 s and 675.6 MiB | `go` (criteria 1-3, 200,000-pair subset, headless runner): all criteria pass |
| featureCounts 2.1.1 | pass / pass / pass / pass (6,477 genes and summary identical) | no index; not measured | not measured | yes: fails without isolation in all 3 browsers (`SharedArrayBuffer is not defined`), also with `-T 1` | (1) pass, (2) pass with isolation, (3) not applicable | `go with a stated caveat` (criteria 1-2, 200,000-pair subset, headless runner): needs COOP/COEP headers or the service worker with one reload |
| DESeq2 1.52.0 (WASM) vs 1.50.2 (native) | fail / not run / not run / not run | no index; not measured | not measured | not tested | tolerance: fail on 1 of 38,203 values (YPR148C pvalue, rel 5.36e-06 against 1e-6); (2) not tested | `not decided` (Node only, one count matrix, one run): browser run not done, and the failure cannot be split from the version difference |

Why DESeq2 is `not decided` and not `no-go` or `go with a stated caveat`:
- Criterion (2) needs a real browser run. DESeq2 did not run in any browser. So the criteria cannot give a `go` of any kind.
- The tolerance test failed on one value. The two sides use different DESeq2 and R versions. This run cannot show whether
  WASM or the version difference causes the failure. So the evidence does not support `no-go` either.
- Other DESeq2 facts from the same run: genes, columns and NA pattern equal (659 NA cells on both sides); max relative
  difference per column: baseMean 0, log2FoldChange 6.18e-09, lfcSE 4.6e-07, stat 4.58e-07, pvalue 5.36e-06, padj 3.57e-07;
  genes with padj < 0.05: 2,042 native and 2,042 WASM.

### Index time and memory per engine (3 runs each, all finished)

| case | engine | run time, s | peak memory, MiB (probe) | WASM heap, MiB |
|---|---|---|---|---|
| HISAT2 yeast 12 Mb | Node | 5.33 to 5.39 | 298.9 to 305.5 (maxRSS) | 113.3 |
| HISAT2 yeast 12 Mb | Chromium | 6.95 to 7.00 | 447.0 to 450.3 (tree Pss) | 113.3 |
| HISAT2 yeast 12 Mb | Firefox | 7.36 to 7.38 | 669.8 to 671.7 (tree Pss) | 113.3 |
| HISAT2 yeast 12 Mb | WebKit | 7.43 to 7.49 | 667.0 to 675.6 (tree Pss) | 113.3 |
| HISAT2 Af293 29 Mb | Node | 11.71 to 12.14 | 412.3 to 424.7 (maxRSS) | 134.7 |
| HISAT2 Af293 29 Mb | Chromium | 14.60 to 14.66 | 531.1 to 533.9 (tree Pss) | 134.7 |
| HISAT2 Af293 29 Mb | Firefox | 15.29 to 16.21 | 717.3 to 754.2 (tree Pss) | 134.7 |
| HISAT2 Af293 29 Mb | WebKit | 15.79 to 15.97 | 834.2 to 835.3 (tree Pss) | 134.7 |
| kallisto yeast transcripts | Node | 19.34 to 19.58 | 274.7 to 289.8 (maxRSS) | 117.1 |
| kallisto yeast transcripts | Chromium | 20.78 to 21.14 | 369.9 to 371.7 (tree Pss) | 117.1 |
| kallisto yeast transcripts | Firefox | 36.80 to 37.03 | 649.7 to 653.5 (tree Pss) | 117.1 |
| kallisto yeast transcripts | WebKit | 21.38 to 22.28 | 823.9 to 828.8 (tree Pss) | 117.1 |
| kallisto Af293 transcripts | Node | 43.96 to 44.18 | 411.0 to 415.4 (maxRSS) | 204.4 |
| kallisto Af293 transcripts | Chromium | 46.99 to 47.15 | 481.1 to 491.9 (tree Pss) | 204.4 |
| kallisto Af293 transcripts | Firefox | 83.89 to 84.41 | 770.2 to 781.1 (tree Pss) | 204.4 |
| kallisto Af293 transcripts | WebKit | 49.42 to 49.79 | 1,170.1 to 1,176.5 (tree Pss) | 204.4 |

- Run time is `ToolRunner.run`: it includes worker start and file copies. Source: `validation/browser-spike0/results/index_size_summary.tsv`.
- Native reference (cluster, one run, `/usr/bin/time -v`): HISAT2 yeast 8.83 s, 107.0 MiB; HISAT2 Af293 24.63 s, 128.2 MiB;
  kallisto yeast 21.7 s, 136.5 MiB; kallisto Af293 not measured.
- HISAT2 index files: identical SHA-256 in all 12 runs per genome. Total bytes equal the native index (yeast 22,818,973;
  Af293 48,273,818). Native index bytes were not kept, so no byte comparison with native.
- No run failed, ran out of memory, crashed or timed out (48 of 48 finished).

### What can enter the browser build, what stays CLI-only

- Criteria met (scope as in the summary): kallisto, HISAT2, featureCounts (featureCounts only with isolation).
- Criteria not decided: DESeq2.
- CLI-only because of a measured failure: none.
- CLI-only by design, not tested: Salmon (excluded by `2026-11-20-salmon-support-design.md`).
- The parent design rule "if 0b fails for HISAT2, genome mode becomes CLI-only" does not apply: 0b passed for HISAT2.
- The per-sample read size cap (parent design 11d and 11e, spike 0f) was not measured. The fixture is 200,000 pairs.

## Isolation and downloads

### Isolation (part 0c, run 44963830025)

Modes: `none` (no header); `coop-coep` (COOP `same-origin`, COEP `require-corp`, CORP `same-origin`); `service-worker`
(no server header; the page registers a worker that adds the three headers, then reloads). Local server on `http://127.0.0.1`.

| mode | `crossOriginIsolated` (all 3 browsers) | featureCounts `-T 1`, `-T 2`, `-T 4` | kallisto index + quant | HISAT2 build + align |
|---|---|---|---|---|
| none | false | error in all 3 browsers (`SharedArrayBuffer is not defined`) | ok in all 3 browsers | ok in all 3 browsers |
| coop-coep | true | ok in all 3 browsers | ok | ok |
| service-worker | false on first load, true after one reload | ok after the reload | ok | ok |

- Output SHA-256 of every ok step equals the part 0a output, in every mode and browser, for `-T 1`, `-T 2` and `-T 4`.
- The test does not show that more than one featureCounts thread ran in parallel. No speed-up is claimed.
- Service worker: the worker file is a copy written for this test (whole-origin scope). Its fetch logic follows the upstream
  `examples/rna-analysis/isolation-worker.js`. The upstream example page was not run. Isolation needed one page reload after
  registration in all 3 browsers.
- `new WebAssembly.Memory({shared: true})` succeeded in every mode, also without isolation. The failure comes from the
  missing `SharedArrayBuffer` global.

### Cross-origin downloads (part 0d, run 44963830025)

Every request: `fetch(url, {mode: 'cors', credentials: 'omit'})` from a page on `http://127.0.0.1:PORT`. 13 requests x 3 modes x 3 browsers = 117 rows
(`results/cors_results.tsv`). 108 rows readable, 9 rows failed (all Ensembl Genomes).

| source | requests | result in all 3 modes and 3 browsers | headers seen |
|---|---|---|---|
| NCBI FTP over https (`ftp.ncbi.nlm.nih.gov`) | genome `.fna.gz`, Range 0-1023, GFF3 `.gff.gz` | readable; status 200 and 206; full-file SHA-256 equals `manifest.tsv` | `Access-Control-Allow-Origin: *`; no CORP seen |
| ENA (`www.ebi.ac.uk` portal API, `ftp.sra.ebi.ac.uk`) | filereport; small FASTQ HEAD and GET (1,280,744 bytes); large FASTQ HEAD and 2 Range GETs | readable; status 200 and 206 with the exact byte counts | `Access-Control-Allow-Origin: *`; no CORP seen; `Content-Range` not readable from script |
| Ensembl (`ftp.ensembl.org`) | yeast Mito GFF3, Mito FASTA, `CHECKSUMS` | readable; status 200 | `Access-Control-Allow-Origin: *`; no CORP seen |
| Ensembl Genomes (`ftp.ensemblgenomes.org`) | `pub/fungi/current/gff3/` | failed: TLS certificate name mismatch in all 3 browsers and in `node:https`. This is not a CORS result | none |
| FungiDB | none | NOT TESTED (needs an API key) | not tested |

- Under `require-corp` (modes `coop-coep` and `service-worker` after reload), every readable URL stayed readable in `cors` mode.
  No-cors requests to these hosts were not tested. None of these hosts sent a CORP header.
- The large ENA file (669,296,728 bytes) was not downloaded in full. Only its HEAD and two byte ranges were read.
- An empty header column means Playwright saw no such header. It is not a verified absence for every host.

## Hosting implications (facts only)

- featureCounts in this build needs cross-origin isolation. On a static host, that means COOP/COEP response headers, or a
  service worker that adds them after one reload (measured on `http://127.0.0.1` only).
- kallisto and HISAT2 ran without isolation headers in all 3 browsers.
- Whether GitHub Pages can send COOP/COEP headers was not measured.
- The service worker was not tested on a real https host. No real host was tested.
- The placeholder site https://qrnamap.fungalgenomes.org/ is live (HTTP 200, recorded after its deploy on 2026-11-21). It is
  deployed from `site/` on `main` by `.github/workflows/pages.yml` (commit `f834efc`). The site files set no isolation header
  and register no service worker. The spike runs did not test this site. Its live headers are listed under "Not tested" below.
- The DESeq2 package zip is about 438 MB (inner zip 437,694,988 bytes). Host size limits were not checked.

## Failed, skipped, not tested and UNVERIFIED items

Failed:
- DESeq2 parity in Node: 1 of 38,203 values outside the stated tolerance (log: `results/run-44944424373/deseq2-node/`, job
  `deseq2-node` of run 44944424373). Version confound: 1.52.0/R 4.6.0 (WASM) against 1.50.2/R 4.5.3 (native).
- Upstream DESeq2 build: run 44939590596 failed (source URL returns HTTP 404). The package was not rebuilt.
- Ensembl Genomes request: TLS certificate name mismatch (not CORS), 9 of 9 combinations.
- featureCounts without isolation: error in all 3 browsers. This is the isolation finding. It matches the parent design statement "featureCounts needs the headers".

Skipped:
- DESeq2 in Chromium, Firefox and WebKit (Node only in part 0a). So DESeq2 browser memory and OOM behaviour are UNVERIFIED.
- Salmon: out of scope. It ran in the upstream `quant.yml` matrix only because that workflow was not edited. Its artifacts were not used.

Not tested (UNVERIFIED):
- A real https host, GitHub Pages headers, and a service worker on a non-localhost origin. Whether GitHub Pages can send COOP/COEP headers is not measured.
- The live placeholder site in a browser. One header check exists, made by the reviewer, not by a spike run: on 2026-11-22 00:10 UTC,
  `curl -sI https://qrnamap.fungalgenomes.org/` gave HTTP/2 200, server GitHub.com via Fastly, `access-control-allow-origin: *`,
  `cache-control: max-age=600`, and no `cross-origin-opener-policy`, `cross-origin-embedder-policy` or `cross-origin-resource-policy` header.
- No-cors requests (for example plain `<img>` or `<script>`) to NCBI, ENA or Ensembl under COEP `require-corp`.
- DESeq2 with and without cross-origin isolation (DESeq2 ran in Node only).
- The upstream `examples/` pages, COEP `credentialless`, headed browsers, mobile browsers, user devices with less than 16 GB.
- Full-size read downloads in the browser (the 669 MB ENA file), and the browser read count check against ENA.
- The per-sample read size cap and chunked processing (spike 0f), and browser gzip decompression speed.
- Genomes larger than 29 Mb; multi-threaded index builds.
- featureCounts on the WASM HISAT2 SAM (it read the native SAM; the two SAMs are equal apart from `@PG`).
- FungiDB downloads; other Ensembl Genomes hosts or paths; large Ensembl files; other `Origin` values.
- A same-version native DESeq2 reference (DESeq2 1.52.0 with R 4.6.0).
- The licence table of the other 41 R packages in the DESeq2 lock file.
- Repeatability of native HISAT2 and kallisto output (one native run each).

Unexplained observations (cause UNVERIFIED):
- kallisto index bytes differ between all WASM runs (also within one engine) and from the native index. Sizes vary by up to
  about 500 bytes. Quant output from these indexes is identical. Whether two native index runs give the same bytes was not tested.
- kallisto Af293 has no native reference (index time and memory not measured natively).
- WebKit: after service worker registration and before the reload, a dedicated worker reported `crossOriginIsolated` true
  while `SharedArrayBuffer` was not defined in it. The page was not isolated.
- Chromium logged `net::ERR_ABORTED` for the 2 ENA HEAD requests while the page fetch got status 200.
- HISAT2 yeast index build in Node took 5.3 s in part 0b and 8.72 s in part 0a (different runner instances).
- The failed upstream DESeq2 build log has 14 other `curl: (22) ... 404` lines (15 in total). Why only one was fatal is UNVERIFIED.

## Decisions for the user

1. Licences: the 7 questions at the end of the licence section below (page licence, hosting of the bundles, source archives,
   notices, upstream licence clarification, DESeq2 package licence table, pin by release asset or by commit).
2. Hosting: which host for the page, and how it gets cross-origin isolation for featureCounts (host headers, or the service
   worker with one reload). GitHub Pages header ability and a real https host are not measured.
3. DESeq2: which path to take to reach a decision.
4. Genome mode limits: which genome size and read size limits the page shows.
5. Whether a kallisto index made in the browser may differ in bytes from run to run, given identical quant output in this test.
6. The first page: whether and when to start it. This file gives no recommendation beyond the decision table.

## Licences (spike 0e)

This section lists facts from licence texts. It gives no legal conclusion. The decisions are in the last list.

### Scope and sources read

- Repository: `comparativechrono/genomics-wasm`, pinned commit `0b015f63e8048eb57b72ffd2706b79b803862b7c`.
- Files read through `gh api` (read-only) at the pin: `LICENSE`, `THIRD_PARTY.md`, `README.md`, and for kallisto, HISAT2, featureCounts and DESeq2 the files in `tools/<tool>/`.
- Bundles that the spike used, listed with `unzip -l` or `find`, under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/baseline/`:
  - `quant/quant-compiled-kallisto/dist/kallisto/` (unpacked build output).
  - `ci/tool-packages/hisat2-2.2.2-wasm.zip`.
  - `rna-analysis/rna-analysis-packages/featurecounts-2.1.1-wasm.zip`.
  - `release/deseq2-package.zip` (contains one inner file, `deseq2-webr-1.52.0.zip`, about 438 MB).
- I did not run any file from the repository or the bundles.
- I did not check that the spike's bundles were built from the pinned commit. The kallisto `sources.lock.json` in the bundle has the same `revision` as the pin.

### Repository licence (as written)

- `LICENSE` (repo root), first line: "MIT License". Copyright line: "Copyright (c) 2026 genomics-wasm contributors".
- `LICENSE`, last paragraph: "This licence covers original build infrastructure, JavaScript runtime and validation code only. Upstream programs, libraries, and code derived from them retain their own licences."
- `README.md`, heading "Licence": "compiled packages are not all MIT-licensed."
- `THIRD_PARTY.md`, last paragraph: "The repository-level MIT licence applies only to original infrastructure, runtime and validation code."
- `THIRD_PARTY.md` heading "Upstream software", first sentence: the WebAssembly packages "contain upstream software, not reimplementations of the algorithms."
- `THIRD_PARTY.md`, last paragraph: "Each binary package includes the relevant notices, build recipes, portability changes and the corresponding upstream sources."
- GitHub API field `license` for the repository: key `other`, name "Other", `spdx_id` `NOASSERTION`.
- The texts say the repository MIT licence covers the repository's own code. They say it does not cover the tool programs. They do not say that the compiled bundles are MIT.
- The earlier note in the qRNAmap design that "genomics-wasm's own LICENSE is MIT" matches the root `LICENSE` text. The texts do not say why GitHub shows "Other".

### Per tool

#### kallisto 0.52.0

- Tool licence as written:
  - `tools/kallisto/README.md`, last paragraph: "kallisto and Bifrost use the BSD 2-Clause licence."
  - `tools/kallisto/LICENSE`, first line: "BSD 2-Clause License". Copyright line: "Copyright (c) 2017, Nicolas Bray, Harold Pimentel, Páll Melsted and Lior Pachter".
  - `THIRD_PARTY.md`, entry "kallisto": "BSD-2-Clause".
- Bundled library licences (from `tools/kallisto/THIRD_PARTY_LICENSES.txt`; headings are file names in that file):
  - `LICENSE.bifrost`: "BSD 2-Clause License", copyright "[2018], [Pall Melsted, Guillaume Holley]".
  - `LICENSE.CRoaring`: "Apache License Version 2.0, January 2004".
  - `LICENSE.BooPHF`: "The MIT License (MIT)", copyright 2015 Guillaume Rizk.
  - `LICENSE.libdivide`: "libdivide is made available under two licenses" (zlib and Boost).
  - `src/kseq.h`: "The MIT License", copyright Attractive Chaos.
  - `src/unordered_dense.h`: header says "ankerl::unordered_dense" version 3.1.1 and `SPDX-License-Identifier: MIT`.
  - `ext/bifrost/src/libpopcnt.h`: notice text, copyright 2016 - 2017 Kim Walisch.
  - `ext/bifrost/src/wyhash.h`: "released under the Unlicense".
  - File text, line 3: "HDF5, HTSlib and zlib-ng are not linked into this build."
  - zlib 1.3.1: pinned in `sources.lock.json`. The text "Zlib license is supplied separately as LICENSE.zlib" is in `THIRD_PARTY_LICENSES.txt`. `LICENSE.zlib` is not in `tools/kallisto/` at the pin. It is in the built `dist/kallisto/`.
  - Emscripten: `sources.lock.json` pins `"emscripten": "3.1.64"`. `tools/kallisto/` has no `LICENSE.emscripten` at the pin. The built `dist/kallisto/` has no `LICENSE.emscripten` either.
- Source versions (`tools/kallisto/sources.lock.json`): kallisto `0.52.0`, revision `4e9f29cf3b021260415430c057a22469ca081391`; zlib `1.3.1`; Emscripten `3.1.64`.
- Redistribution text: `THIRD_PARTY_LICENSES.txt` line 3: "The exact upstream source archive accompanies this package." `README.md`: "the exact source archives accompany each package."
- Files in the spike's built `dist/kallisto/`: `LICENSE`, `LICENSE.bifrost`, `LICENSE.zlib`, `THIRD_PARTY_LICENSES.txt`, `sources.lock.json`, `browser.patch`, `build-manifest.json`, `kallisto-4e9f29cf...tar.gz`, `zlib-1.3.1.tar.gz`, `kallisto.mjs`, `kallisto.wasm`, `serial-thread.h`. Separate `LICENSE.BooPHF`, `LICENSE.CRoaring` and `LICENSE.libdivide` files are not in that folder. Their texts are inside `THIRD_PARTY_LICENSES.txt`.
- No kallisto release zip was in the baseline folder. I only listed the unpacked build output.

#### HISAT2 2.2.2

- Tool licence as written:
  - `tools/hisat2/README.md`, last paragraph: "HISAT2 and the port are distributed under GPL-3.0-or-later."
  - `tools/hisat2/LICENSE`: "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007" (full text, 674 lines).
  - `THIRD_PARTY.md`, entry "HISAT2": "including code from Bowtie 2; GPL-3.0-or-later."
- Bundled component notices (`tools/hisat2/THIRD_PARTY_LICENSES.txt`):
  - Line 3: "Full upstream sources and GPL-3.0 license are included in this package."
  - `tinythread.h`: copyright 2010-2012 Marcus Geelnard, with a zlib-style permission text; modified by Nigel Dyer 2016.
  - `third_party/cpuid.h`: Free Software Foundation, with the "GCC Runtime Library Exception, version 3.1" (full exception text in the file).
  - Bowtie 2 code: "Copyright 2011, Ben Langmead", under "GNU General Public License ... version 3".
  - `qsufsort.c`: "Copyright 1999, N. Jesper Larsson"; "This software may be used freely for any purpose. However, when distributed, the original source must be clearly stated".
  - The file has no zlib, htslib or Emscripten section. `sources.lock.json` lists no zlib dependency.
- Source versions (`tools/hisat2/sources.lock.json`): HISAT2 `2.2.2`, revision `99583d7536b9ee017ac07de8834017a3bf99a2fe`; Emscripten `3.1.64`.
- Redistribution text: `README.md`: "The package includes that archive, the patch, dispatcher, license notices and build script."
- In `hisat2-2.2.2-wasm.zip`: the root `LICENSE` starts with "MIT License". `dist/hisat2/LICENSE` and `tools/hisat2/LICENSE` start with "GNU GENERAL PUBLIC LICENSE". The zip also holds `THIRD_PARTY.md`, `source-manifest.json`, `dist/hisat2/THIRD_PARTY_LICENSES.txt`, and `build/downloads/hisat2-99583d75...tar.gz`.

#### featureCounts 2.1.1 (Subread)

- Tool licence as written:
  - `tools/featurecounts/README.md`, heading "Sources and licences": "Subread is distributed under GPL version 3 or later."
  - Same section: "The build scripts and small port helpers use the repository's MIT licence; this does not change the licence of the compiled program."
  - `tools/featurecounts/LICENSE`: "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007".
  - `tools/featurecounts/THIRD_PARTY_LICENSES.txt`: "Licence: GNU General Public License, version 3 or later."
  - `THIRD_PARTY.md`, entry "featureCounts": "GPL-3.0".
- Bundled libraries:
  - zlib 1.3.1: "Licence: zlib licence. LICENSE.zlib accompanies the compiled distribution". `LICENSE.zlib` is not in `tools/featurecounts/` at the pin. It is in `dist/featurecounts/` in the bundle.
  - Emscripten 3.1.64: `LICENSE.emscripten` first line: "Emscripten is available under 2 licenses, the MIT license and the University of Illinois/NCSA Open Source License." `THIRD_PARTY_LICENSES.txt`: the `mailbox.js` override "derives from src/library_pthread.js" in Emscripten; "their notices remain in that SDK's source distribution."
- Source versions (`tools/featurecounts/sources.lock.json`): Subread `2.1.1` (`subread-2.1.1-source.tar.gz`); zlib `1.3.1`; Emscripten `3.1.64`.
- Redistribution text: `README.md`: "The distribution includes its source lock, build manifest, patches, complete Subread and zlib source archives, and licence notices."
- In `featurecounts-2.1.1-wasm.zip`: root `LICENSE` starts with "MIT License". `dist/featurecounts/LICENSE` starts with "GNU GENERAL PUBLIC LICENSE". The zip also has `dist/featurecounts/LICENSE.emscripten`, `LICENSE.zlib`, `THIRD_PARTY_LICENSES.txt`, `sources.lock.json`, and three source tarballs: `build/downloads/subread-2.1.1-source.tar.gz`, `build/downloads/zlib-1.3.1.tar.gz` and `build/rna-analysis-sources/hisat2-source.tar.gz`. It also has test-data licence files `datasets/rna-analysis/provenance/HISAT2_LICENSE` and `datasets/rna-analysis/provenance/yeast-UPSTREAM_LICENSE.txt`.

#### DESeq2 1.52.0 (webR, R 4.6.0)

- Files at the pin: `tools/deseq2/` has `README.md`, `sources.lock.json`, `licenses/Noto-OFL.txt`, `licenses/Noto-font-notices.txt`. It has no `LICENSE`, no `LICENSE*` file and no `THIRD_PARTY_LICENSES.txt` at the pin (not found at the pin).
- Tool licence as written:
  - `tools/deseq2/README.md`, last paragraph: "DESeq2 is LGPL-3.0-or-later. The webR binary distribution is GPL-3.0 and dependencies retain their individual licences."
  - Same paragraph: "The archive includes upstream licence notices, corresponding sources and build information. Original JavaScript and validation infrastructure in this repository are MIT-licensed."
  - `sources.lock.json`, package `DESeq2` 1.52.0: `"license": "LGPL (>= 3)"`.
  - `THIRD_PARTY.md`, entry "DESeq2": "LGPL-3.0-or-later."
- GPL parts, as written:
  - `THIRD_PARTY.md`, entry "webR": "distribution binaries are GPL-3.0. Build scripts and additional application software are MIT unless stated otherwise upstream. R and bundled libraries retain their respective licences."
  - `THIRD_PARTY.md`, entry "R": "GPL-2.0-or-later."
  - `sources.lock.json`, `r` section: `"license": "GPL (>= 2)"`; R source `R-4.6.0.tar.gz`.
  - `sources.lock.json`, package `locfit` 1.5-9.12: `"license": "GPL (>= 2)"`. `THIRD_PARTY.md`: "`locfit` is GPL-2.0-or-later; the build changes its compiler configuration without changing the statistical implementation."
  - `sources.lock.json` lists 43 packages, each with a `license` field. I read only DESeq2, locfit and R. I did not tabulate the other 41 licences.
- Other bundled licences: `tools/deseq2/licenses/Noto-font-notices.txt`: Noto fonts "licensed under the SIL Open Font License, Version 1.1". The inner zip has `dist/deseq2/webr/LICENSE.md`, a webR "Software Licensing" list (R GPL v2 and v3, PCRE2 BSD, XZ public domain, jsTree MIT, libgfortran GPL v2 and v3, Dragonegg GPL v2, LLVM Flang Apache 2.0, xterm.js MIT, and more; I read the first 8 entries).
- Source versions (`tools/deseq2/sources.lock.json`): DESeq2 `1.52.0`; R `4.6.0`; webR `0.6.0` (npm `webr-0.6.0.tgz` and `webr-0.6.0-source.tar.gz`); Emscripten source `5.0.7` appears under `build/deseq2/sources/` in the bundle.
- Redistribution text: see the README paragraph above.
- In `deseq2-package.zip`: one file, `deseq2-webr-1.52.0.zip`. That inner zip has 313 entries. Its root `LICENSE` starts with "MIT License" and matches the repository `LICENSE`. It has `THIRD_PARTY.md` (119 lines; it differs from the repository file at the pin, I did not diff the content), `dist/deseq2/webr/LICENSE.md`, `tools/deseq2/licenses/Noto-*.txt`, and 75 files under `build/deseq2/sources/` (including `R-4.6.0.tar.gz` and `webr-0.6.0-source.tar.gz`). It has no `LICENSE` file for DESeq2 itself and no `THIRD_PARTY_LICENSES.txt`.

### Summary table (texts as written)

| Tool | Version | Tool licence as written | Source archive in bundle |
|---|---|---|---|
| kallisto | 0.52.0 | BSD 2-Clause | yes (kallisto, zlib) |
| HISAT2 | 2.2.2 | GPL-3.0-or-later | yes (HISAT2) |
| featureCounts | 2.1.1 | GPL-3.0 or later | yes (Subread, zlib) |
| DESeq2 | 1.52.0 | LGPL-3.0-or-later; webR binaries GPL-3.0; R GPL-2.0-or-later | yes (75 source files, including R and webR) |

### What the texts say and do not say

- All four READMEs and `THIRD_PARTY.md` say each package ships notices and the upstream sources. No text states a rule about hosting the bundles on a web page or about fetching them from another site.
- No text in the files I read says how a user or a third party may redistribute the built bundles. The repository `LICENSE` says it does not cover upstream code.
- The texts give MIT for the repository code and the tool licences above for the tools. They do not say which licence applies to a zip that holds both. What I found in the bundles: HISAT2 zip and featureCounts zip each have an MIT root `LICENSE` and a GPL `LICENSE` in `dist/`. The DESeq2 inner zip has an MIT root `LICENSE` and `dist/deseq2/webr/LICENSE.md` (a webR licence list); it has no GPL `LICENSE` file. For kallisto, no zip exists in the baseline folder; the unpacked build output has `LICENSE` (BSD 2-Clause) and no MIT root `LICENSE`.
- `THIRD_PARTY.md` does not list the Emscripten runtime as its own entry. The featureCounts files name it (MIT and University of Illinois/NCSA). The kallisto and HISAT2 files at the pin do not name it.
- GitHub shows "Other" (`NOASSERTION`). The root `LICENSE` text says MIT.
- The DESeq2 directory at the pin has no tool `LICENSE` file. The README gives the licence statement.

### qRNAmap side (facts)

- qRNAmap has a `LICENSE` file at the repository root. Its text is "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007".
- `Cargo.toml` (workspace) has `license = "GPL-3.0-only"`. `README.md` line 4 says "Licence: GPL-3.0-only."
- `docs/superpowers/specs/2026-11-19-qrnamap-design.md` (lines 174 and 178) says: qRNAmap is GPL-3.0-only (user, 2026-11-19); qmutfind is GPL-3.0-only, verified in its `Cargo.toml`.
- The same design (line 210) says "A page that serves GPL binaries also offers the matching source archives, which the bundles include." It marks this as proposed.
- Size of the DESeq2 package zip: about 438 MB (inner zip 437,694,988 bytes). GitHub Pages size limits are not checked in this section.

### Questions for the user (decision)

1. Which licence does qRNAmap use for the page code that loads the tool bundles? The repo file says GPL-3.0-only now. Do you keep that?
2. May the qRNAmap page host (bundle) the tool bundles, or must the page fetch them from the genomics-wasm releases at run time?
3. If the page hosts them: do you want to ship the source archives from each bundle with the page, and where (same site, or a link)?
4. Which attribution text and notice files must the page show for kallisto, HISAT2, featureCounts, DESeq2, webR, R and the Emscripten runtime? Where in the page does it go?
5. Do you want to ask the genomics-wasm author to clarify the repository licence (GitHub shows "Other") and the status of the compiled bundles?
6. For DESeq2: the bundle has no tool `LICENSE` file at the pin, and 41 package licences are not tabulated here. Do you want a full licence table for the 43 R packages before release?
7. Should the project pin the bundles by release asset or by commit, given that the licence texts differ between files at the pin and in the zips (for example `THIRD_PARTY.md` in the DESeq2 zip)?
