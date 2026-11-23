# Browser page version 1: measurements that fix the limits

Date of the runs: 2026-11-22 (UTC). Task 3 of `docs/superpowers/plans/2026-11-22-browser-page-v1.md`.
Spec: `docs/superpowers/specs/2026-11-22-browser-page-v1-design.md` ("To measure before fixing numbers").
Scripts: `validation/browser-page/measure/`. Workflow: `.github/workflows/page.yml`, job `measure`.
Result tables (small TSVs): `validation/browser-page/measure/results/`.

## Scope

- All browser values come from headless Playwright browsers on GitHub-hosted `ubuntu-24.04` runners
  (4 CPUs, 16,367,520 kB RAM, 3,145,724 kB swap). Playwright 1.64.0: Chromium 156.0.8078.4, Firefox 157.0, WebKit 27.2.
  Node v22.23.5. They are not the release builds that users run, and the runner has more memory than many user devices.
- Reads: SRR3396381 (yeast, GSE80357 AmpB_1), the first N pairs. Transcripts: the yeast set of spike 0 (6,478 transcripts).
- Genomes: six genomes of the built-in catalog (list below), from NCBI on 2026-11-22.
- One network (GitHub runners to NCBI, ENA and GitHub) on one day. Download speeds and failures can differ on other days and networks.
- A value marked UNVERIFIED was not measured. A failure is a result and is reported.

## Decisions for the later tasks

- **Changed 2026-11-22 (user decision): N depends on the read length** (section "Page quantification path: peak memory"). The values below hold for
  101-base reads (SRR3396381) only.
- **Preview size N (default): 200,000 read pairs per sample.** Reason: it is the largest tested count that met the N rule
  (all 3 runs finished, 300 s or less per sample, sample peak 1.5 GiB or less) in Chromium, Firefox and WebKit, in mode `plain` and in mode `wplain`.
  At 500,000 pairs, Firefox (`plain`: 1,731 to 1,938 MiB; `wplain`: 2,683 to 2,788 MiB) and WebKit (`plain`: 4,329 to 5,781 MiB; `wplain`: 2,040 to 2,505 MiB) were above 1.5 GiB.
  Time is not the limit at this count: a sample of 200,000 pairs took 4.7 to 10.0 s (preparation plus quant) in the three browsers in these two modes.
- **Hard limit: 200,000 read pairs per sample.** Reason: it is also the largest tested count where all runs finished with a sample peak of 2 GiB or less
  in all three browsers, in every mode. At the next tested count (500,000), WebKit had at least one run above 2 GiB or not finished in every mode, and so had Firefox in modes `wgz` and `wplain`.
  So the page cannot offer a count above N with these data. Larger previews passed in single browsers (tables below), for example Chromium in mode `wgz`
  up to 3,000,000 pairs under 2 GiB; a per-browser limit was not chosen.
- **How the page gives reads to kallisto:** plain FASTQ cut at N records (modes `plain` or `wplain`). The gzip repack of qmutfind `reads.js` (mode `gz`)
  met the N rule only at 100,000 pairs (WebKit preparation peak 1,327 to 1,535 MiB), and in a worker (mode `wgz`) it met the N rule at no count
  (Firefox preparation peak above 5 GiB from 100,000 pairs; WebKit 1,348 to 1,693 MiB at 100,000 pairs).
- **Test-set reads: fetch from ENA in the browser and cut both mates at N records** (read the stream, count N x 4 newlines, keep, cancel).
  The cut worked in all three browsers on the bytes ENA sent. The stream method without a Range header succeeded in 17 of 30 sessions
  (17 of the 19 sessions where ENA answered the Range requests; it failed twice with a network error during the stream). A Range header is not needed for this method.
  ENA did not answer in 11 of 30 browser sessions on the runners, so the page needs a clear error and the test-set tests in CI must not depend on ENA.
  The fallback in the spec (the page fetches GitHub release assets) does not work: all three browsers blocked the request (CORS) in all 9 sessions.
- **Genome menu limit:** no limit is needed for the catalog. The three largest catalog genomes (63.4 to 67.2 Mb) finished download, transcript building
  and index in all engines, with a peak of 1,796.8 MiB or less.
- **kallisto output check:** Firefox wrote a damaged `abundance.tsv` in 4 of 130 kallisto runs (about 3%), with no error (cause UNVERIFIED).
  So a byte-identical comparison of Firefox output with the CLI (the Task 9 acceptance rule) can fail without a fault in the page code.
  Task 9 must re-run or check the file before it reports a mismatch. Needed page behaviour: check each `abundance.tsv`, and re-run the sample once
  when the check fails. This behaviour was not tested (see "Run-to-run output of kallisto").

## Runs

| run | trigger | head | what | result |
|---|---|---|---|---|
| [44984006776](https://github.com/stajichlab/qRNAmap/actions/runs/44984006776) | push | `b682fac` | first check (2 counts, 1 run) | failed: the full ENA FASTQ download stopped (curl error 18, `transfer closed with 647661716 bytes remaining to read`, about 21.6 MB in 480 s). No measurement ran. |
| [44984897294](https://github.com/stajichlab/qRNAmap/actions/runs/44984897294) | push | `c4abdce` | first check (100,000 and 200,000 pairs, 1 run per case) | all jobs passed. Used for the first sight of the Firefox output fault (below). |
| [44985823017](https://github.com/stajichlab/qRNAmap/actions/runs/44985823017) | dispatch | `c4abdce` | full: quant (modes `gz`, `plain`), truncated, genomes; 3 runs per case | see the tables |
| [44986429632](https://github.com/stajichlab/qRNAmap/actions/runs/44986429632) | dispatch | `3a4a968` | determinism (100,000 pairs) and truncated; 3 runs | see the tables |
| [44990656616](https://github.com/stajichlab/qRNAmap/actions/runs/44990656616) | dispatch | `5d307df` | quant (modes `wgz`, `wplain`) in the three browsers, truncated with the GitHub release check; 3 runs | see the tables |

The merged result tables in `results/` come from runs 44985823017, 44986429632 and 44990656616; `results/truncated.tsv` also holds the
truncated-test rows of run 44984897294 (1 session per browser). Run 44984897294 quant rows are not in `results/quant_runs.tsv`.

Every push to `page-v1` also starts `ci.yml` (it has `on: push` with no branch filter). Commits with `[skip ci]` started no workflow.

## Method

The harness follows `validation/browser-spike0/index_size.mjs`. The spike files were not changed.

- `measure.mjs` runs one suite in one engine. Each run uses a new process: a child Node process, or a new browser launch.
  A local static server (`http://127.0.0.1:PORT`) serves the vendored kallisto (`site/vendor/kallisto/`), the harness files and the data.
  It sends no isolation headers. kallisto runs through the vendored `ToolRunner` (`site/vendor/kallisto/runtime/index.mjs`), one worker per command, `-t 1`.
- kallisto arguments are those of `qrnamap-core` `toolcall`: `index -i yeast.idx yeast.transcripts.fa`, then
  `quant -i yeast.idx -o kq --plaintext -t 1 r_1.fastq[.gz] r_2.fastq[.gz]`.
- One quant run: index build, then read preparation (the first N pairs of each mate from the full SRR3396381 files on the local server), then quant.
  A run is "finished" only when kallisto exits with 0, both outputs exist, both mates have N records, and `run_info.json` `n_processed` equals N.
- Read preparation modes (all cut at N records per mate, counted as N x 4 newlines, as in qmutfind `reads.js`):

| mode | where | how | what kallisto gets |
|---|---|---|---|
| `gz` | page main thread (Node: the process) | qmutfind `reads.js` `keepFirstReads` (fflate gunzip, then fflate gzip level 3 of the kept records) | `.fastq.gz` |
| `plain` | page main thread (Node: the process) | qmutfind `gz.js` `VerifyingGunzip` (fflate), kept records stay plain | `.fastq` |
| `wgz` | a dedicated worker; result arrays transferred to the page | as `gz` | `.fastq.gz` |
| `wplain` | a dedicated worker; result arrays transferred to the page | `DecompressionStream('gzip')`, kept records stay plain | `.fastq` |

- Counts: 100,000, 200,000, 500,000, 1,000,000, 2,000,000, 3,000,000, 4,000,000, 6,000,000 and 8,405,534 (all pairs), 3 runs each.
  Stop rules (harness rules, not results): larger counts are skipped when no run of a count finished, or when a run of a count passed 6 GiB.
  A run is stopped when the whole-browser Pss (or the Node child RSS) passes 10 GiB (outcome `mem_guard`). Per-run timeout 30 min. Job budget 300 min.
- Memory probes (as in spike 0; limits in `validation/browser-spike0/README.md`):
  `tree Pss` = sum of `Pss` over all processes below the harness (the whole browser), sampled every 250 ms, per phase
  (index, preparation, quant). Node: `maxRSS` of the child process (whole run). WASM heap: `HEAPU8.length` after `callMain` (`heap_wrap.mjs`).
- The peak of one sample is the higher of the preparation peak and the quant peak (tree Pss). The index is built once per genome,
  so the time of one sample is preparation plus quant. Both include the idle browser (Chromium about 165 to 170 MiB after page load).

### Rules for N and the hard limit (thresholds from the brief; the peak definition changed once, see the quant section)

- N: the largest count where all 3 runs finished in Chromium, Firefox and WebKit, each sample took 300 s or less (preparation plus quant),
  and the sample peak was 1.5 GiB (1,572,864 kB) or less.
  Reason for 1.5 GiB: it keeps 0.5 GiB below 2 GiB, the stated WebAssembly heap limit of this build, and the peaks are whole-browser values
  on a 16 GB runner; user devices can have less memory. This margin is a choice, not a measurement.
- Hard limit: the largest count where all 3 runs finished in all three browsers with a sample peak of 2 GiB (2,097,152 kB) or less.
- `summarize.py` applies both rules per mode (`results/limits.tsv`).

## How kallisto takes the reads, and the memory this needs

From the vendored runtime code (read, not changed):

- `ToolRunner.run` posts the input files to a new worker with `postMessage` without a transfer list (comment in `index.mjs`: "Structured cloning leaves caller-owned input buffers intact"). This makes a copy in the worker.
- The worker (`execute.mjs`) writes each file into the Emscripten file system with `FS.writeFile` before `callMain`. MEMFS keeps its own copy.
- The caller keeps its array until `run` returns, because `ToolRunner.run` holds it.
- So during quant the reads exist up to three times in JavaScript memory: page, worker message, MEMFS. They are outside the WASM heap.
  The runtime also has a `stdin` option, but kallisto paired-end quant reads two files; `stdin` was not tested.
- kallisto reads `.fastq.gz` itself (zlib). The WASM heap after quant was 115,998,720 bytes (110.6 MiB) in every finished run of runs
  44985823017 and 44990656616, at every count and in every engine (the index heap was 122,748,928 bytes, 117.1 MiB). The heap did not grow with the read count. The 2 GiB WASM heap limit was not reached.
- Measured cost per pair: in Chromium the quant peak grew by about 547 bytes per pair from 1,000,000 to 4,000,000 pairs in mode `wgz` (input about 178 bytes per pair as gzip, both mates),
  and by about 1,975 bytes per pair from 1,000,000 to 3,000,000 pairs in mode `plain` (input 555 bytes per pair). That is about 3.1 and 3.6 times the input bytes,
  which fits the three copies above (see the quant tables).
- Failure mode of large plain inputs in Chromium: at 4,000,000 pairs (1,123,555,584 bytes per mate) quant stopped in all 3 runs with
  `DataCloneError: Failed to execute 'postMessage' on 'Worker': Data cannot be cloned, out of memory.` (modes `plain` and `wplain`).
- No streaming of reads into kallisto was tested. A way to avoid the three copies (for example a transfer list, or a worker that reads the
  file itself) is not part of the vendored runtime and was not tested.

## Quantification at several read counts (Step 2)

Runs 44985823017 (modes `gz`, `plain`; Node and the three browsers) and 44990656616 (modes `wgz`, `wplain`; browsers). Each cell: runs finished of 3;
sample time in seconds (preparation plus quant, finished runs); sample peak in MiB (browsers: the higher of the preparation peak and the quant peak,
whole-browser tree Pss; Node: maxRSS of the process). `skipped`: not run because of a stop rule. Per-run values: `results/quant_runs.tsv`;
per-case values with the index time, the quant time and the baseline: `results/quant_summary.tsv`; the rules applied: `results/limits.tsv`.

#### Mode `gz` (main thread; qmutfind reads.js; gzip level 3 to kallisto), run 44985823017

| pairs | node | chromium | firefox | webkit |
|---|---|---|---|---|
| 100,000 | 3/3; 5.3 to 5.6 s; 311 to 318 MiB | 3/3; 5.6 to 5.7 s; 420 to 448 MiB | 3/3; 7.9 to 8.3 s; 667 to 762 MiB | 3/3; 7.4 to 7.8 s; 1,327 to 1,535 MiB |
| 200,000 | 3/3; 9.4 to 9.5 s; 372 to 387 MiB | 3/3; 9.9 s; 479 to 532 MiB | 3/3; 14.2 to 14.3 s; 723 to 777 MiB | 3/3; 10.8 to 12.6 s; 1,431 to 1,981 MiB |
| 500,000 | 3/3; 21.7 to 21.8 s; 585 to 619 MiB | 3/3; 22.5 to 22.6 s; 630 to 666 MiB | 3/3; 32.9 to 33.6 s; 1,086 to 1,135 MiB | 3/3; 144.4 to 145.7 s; 7,032 to 8,399 MiB |
| 1,000,000 | 3/3; 42.5 to 42.7 s; 921 to 927 MiB | 3/3; 43.8 to 44.0 s; 936 to 974 MiB | 3/3; 64.5 to 65.5 s; 1,500 to 1,546 MiB | 3/3; 105.9 to 165.9 s; 7,231 to 8,568 MiB |
| 2,000,000 | 3/3; 81.7 to 83.7 s; 1,529 to 1,539 MiB | 3/3; 86.7 to 87.3 s; 1,545 to 1,551 MiB | 3/3; 126.8 to 127.8 s; 2,198 to 2,496 MiB | 3/3; 203.6 to 205.5 s; 7,849 to 8,809 MiB |
| 3,000,000 | 3/3; 123.8 to 127.1 s; 2,145 to 2,148 MiB | 3/3; 129.9 to 131.3 s; 2,154 to 2,159 MiB | 3/3; 192.1 to 193.8 s; 2,974 to 3,259 MiB | 2/3; 221.5 to 225.6 s; 8,166 to 10,282 MiB; mem_guard |
| 4,000,000 | 3/3; 164.4 to 166.7 s; 2,721 to 2,747 MiB | 3/3; 174.0 to 174.6 s; 2,761 to 2,763 MiB | 3/3; 260.1 to 261.1 s; 3,433 to 4,227 MiB | skipped |
| 6,000,000 | 3/3; 250.5 to 254.4 s; 3,958 to 3,962 MiB | 3/3; 262.6 to 263.8 s; 5,040 to 5,043 MiB | 3/3; 399.9 to 406.8 s; 5,481 to 6,033 MiB | skipped |
| 8,405,534 | 3/3; 359.4 to 361.0 s; 5,448 to 5,457 MiB | 3/3; 369.6 to 372.7 s; 6,945 to 6,948 MiB | 3/3; 570.1 to 577.9 s; 7,778 to 8,548 MiB | skipped |

#### Mode `plain` (main thread; fflate; plain FASTQ to kallisto), run 44985823017

| pairs | node | chromium | firefox | webkit |
|---|---|---|---|---|
| 100,000 | 3/3; 4.6 s; 432 to 454 MiB | 3/3; 3.7 to 3.8 s; 597 to 634 MiB | 3/3; 5.9 to 6.0 s; 788 to 824 MiB | 3/3; 5.6 to 5.9 s; 1,034 to 1,129 MiB |
| 200,000 | 3/3; 7.7 to 7.8 s; 687 to 699 MiB | 3/3; 6.2 to 6.3 s; 722 to 753 MiB | 3/3; 9.7 to 10.0 s; 1,105 to 1,109 MiB | 3/3; 7.8 to 9.2 s; 1,234 to 1,327 MiB |
| 500,000 | 3/3; 16.8 to 17.2 s; 1,266 to 1,276 MiB | 3/3; 13.8 s; 1,248 to 1,260 MiB | 3/3; 21.6 to 21.9 s; 1,731 to 1,938 MiB | 3/3; 53.8 to 64.5 s; 4,329 to 5,781 MiB |
| 1,000,000 | 3/3; 31.6 to 31.9 s; 2,216 to 2,222 MiB | 3/3; 26.5 to 26.6 s; 2,175 to 2,177 MiB | 3/3; 41.5 to 42.4 s; 3,106 to 3,242 MiB | 3/3; 68.3 to 72.8 s; 4,884 to 6,299 MiB |
| 2,000,000 | 3/3; 62.8 to 63.3 s; 4,166 to 4,173 MiB | 3/3; 51.7 to 52.2 s; 4,058 to 4,063 MiB | 3/3; 85.3 to 88.2 s; 5,289 to 5,879 MiB | 3/3; 77.2 to 90.3 s; 6,057 to 6,459 MiB |
| 3,000,000 | 3/3; 92.8 to 93.9 s; 6,091 to 6,100 MiB | 3/3; 77.2 to 77.6 s; 5,942 to 5,943 MiB | 3/3; 127.6 to 132.7 s; 7,227 to 8,051 MiB | 0/3; 8,772 to 8,777 MiB; timeout |
| 4,000,000 | 3/3; 121.5 to 125.3 s; 8,031 to 8,047 MiB | 0/3; 4,262 to 4,481 MiB; oom | skipped | skipped |
| 6,000,000 | skipped | skipped | skipped | skipped |
| 8,405,534 | skipped | skipped | skipped | skipped |

#### Mode `wgz` (worker; qmutfind reads.js; gzip level 3 to kallisto), run 44990656616

| pairs | chromium | firefox | webkit |
|---|---|---|---|
| 100,000 | 3/3; 5.6 to 5.7 s; 376 to 380 MiB | 3/3; 9.6 to 9.9 s; 5,161 to 5,354 MiB | 3/3; 9.0 to 9.1 s; 1,348 to 1,693 MiB |
| 200,000 | 3/3; 9.7 to 10.1 s; 442 to 446 MiB | 3/3; 15.4 to 15.8 s; 5,171 to 5,221 MiB | 3/3; 12.3 to 13.3 s; 1,535 to 1,799 MiB |
| 500,000 | 3/3; 22.5 to 22.6 s; 574 to 583 MiB | 3/3; 35.0 to 35.1 s; 5,116 to 5,372 MiB | 0/3; 8,521 to 8,576 MiB; timeout |
| 1,000,000 | 3/3; 43.9 to 44.1 s; 836 to 840 MiB | 3/3; 65.2 to 67.5 s; 5,229 to 5,540 MiB | skipped |
| 2,000,000 | 3/3; 86.3 to 87.1 s; 1,362 to 1,366 MiB | 3/3; 132.7 to 133.5 s; 4,980 to 5,785 MiB | skipped |
| 3,000,000 | 3/3; 128.8 to 129.3 s; 1,886 to 1,888 MiB | 3/3; 200.8 to 203.8 s; 4,867 to 4,881 MiB | skipped |
| 4,000,000 | 3/3; 170.6 to 173.6 s; 2,402 to 2,403 MiB | 3/3; 266.8 to 269.4 s; 4,085 to 4,293 MiB | skipped |
| 6,000,000 | 3/3; 260.0 to 264.3 s; 4,508 to 4,511 MiB | 3/3; 409.3 to 412.5 s; 5,146 to 5,579 MiB | skipped |
| 8,405,534 | 3/3; 371.1 to 372.2 s; 6,202 to 6,204 MiB | 3/3; 586.6 to 593.0 s; 6,354 to 6,390 MiB | skipped |

#### Mode `wplain` (worker; DecompressionStream; plain FASTQ to kallisto), run 44990656616

| pairs | chromium | firefox | webkit |
|---|---|---|---|
| 100,000 | 3/3; 2.9 to 3.1 s; 525 to 527 MiB | 3/3; 5.0 to 5.2 s; 878 to 929 MiB | 3/3; 3.8 to 3.9 s; 937 to 980 MiB |
| 200,000 | 3/3; 4.7 to 5.0 s; 652 to 657 MiB | 3/3; 8.6 to 8.7 s; 1,366 to 1,375 MiB | 3/3; 5.9 to 6.0 s; 1,206 to 1,233 MiB |
| 500,000 | 3/3; 10.1 to 10.2 s; 1,117 to 1,118 MiB | 3/3; 18.7 to 18.9 s; 2,683 to 2,788 MiB | 3/3; 13.2 to 13.7 s; 2,040 to 2,505 MiB |
| 1,000,000 | 3/3; 19.9 to 20.4 s; 1,913 to 1,914 MiB | 3/3; 37.9 to 38.9 s; 5,561 to 5,747 MiB | 3/3; 25.3 to 27.5 s; 3,307 to 3,312 MiB |
| 2,000,000 | 3/3; 39.0 to 39.8 s; 3,529 to 3,530 MiB | 3/3; 79.1 to 80.7 s; 7,113 to 7,216 MiB | 3/3; 53.2 to 60.8 s; 6,067 to 6,070 MiB |
| 3,000,000 | 3/3; 57.6 to 59.7 s; 5,144 to 5,145 MiB | skipped | 0/3; 8,692 to 8,788 MiB; timeout |
| 4,000,000 | 0/3; 3,317 to 3,565 MiB; oom | skipped | skipped |
| 6,000,000 | skipped | skipped | skipped |
| 8,405,534 | skipped | skipped | skipped |

Index build in the same runs (yeast transcripts, before each quant): Chromium 20.2 to 21.8 s, Firefox 35.2 to 39.2 s, WebKit 20.8 to 25.3 s, Node 20.6 to 30.6 s
(`index_s` in `results/quant_runs.tsv`).

Result of the rules (`results/limits.tsv`):

| mode | N rule: largest count in all three browsers | hard rule: largest count | first count failing the N rule (Chromium, Firefox, WebKit) |
|---|---|---|---|
| `gz` | 100,000 | 200,000 | 2,000,000, 1,000,000, 200,000 |
| `plain` | 200,000 | 200,000 | 1,000,000, 500,000, 500,000 |
| `wgz` | none | none | 3,000,000, 100,000, 100,000 |
| `wplain` | 200,000 | 200,000 | 1,000,000, 500,000, 500,000 |

- Failures: Chromium `plain` and `wplain` at 4,000,000 pairs: `DataCloneError ... out of memory` (3 of 3 each). WebKit `plain` and `wplain` at 3,000,000 and
  `wgz` at 500,000: timeout after 30 min (3 of 3 each), with 8,521 to 8,788 MiB Pss. WebKit `gz` at 3,000,000: one run stopped by the 10 GiB guard.
- The sample peak rule was first written with the quant peak only. After run 44985823017 showed WebKit preparation peaks of 7,032 to 9,063 MiB (mode `gz`, 500,000 pairs and more), the rule was
  changed to the higher of the preparation and quant peaks, before run 44990656616. With the quant peak only, N would also be 200,000 in modes `gz`, `plain`
  and `wplain` (`quant_peak_mib` in `results/quant_summary.tsv`).
- The tree Pss includes memory that the browser has not freed yet. The harness did not force garbage collection. In WebKit and Firefox the baseline before quant
  was sometimes much higher than after page load (for example WebKit `gz` 1,000,000 pairs: 488.5 to 3,572.0 MiB). Behaviour under memory pressure was not tested.
- Read counts: every finished run had `n_processed` equal to N and N records in each mate.
- Same reads, same result: for each count, all finished runs of all engines and modes gave one `abundance.tsv` SHA-256 (`results/quant_abundance_agreement.tsv`,
  runs 44985823017 and 44990656616). The Firefox output fault below did not occur in these runs.
- Reads in CI came from the NCBI SRA file of SRR3396381 (`prepare_reads.sh`, fasterq-dump 3.2.0). Both mates have 8,405,534 records (equal to the ENA `read_count`).
  The sequence and quality lines of the records in the first 1 MiB of each ENA file equal the first records of the SRA-derived files
  (13,696 and 13,590 records compared; `results/reads_check.run-44985823017.tsv`). The header lines were not compared.

## Truncated gzip from ENA, and the test-set read fetch

Suite `truncated` (browsers). Runs 44984897294 (1 session per browser), 44985823017, 44986429632 and 44990656616 (3 sessions per browser each): 30 sessions.
URLs: `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR339/001/SRR3396381/SRR3396381_{1,2}.fastq.gz`. Each session asks for bytes 0 to 8,388,607 of each mate
(`Range`, HTTP 206), then runs the tests below on those bytes. One test streams the whole URL without `Range` and stops early.
Values below are identical in every session where the fetch worked, unless counts are given. In the passing stream tests, 3,839,689 to 4,303,689 bytes
were read before the cancel, and the plain bytes equalled those of the Range bytes.

| test | Chromium | Firefox | WebKit |
|---|---|---|---|
| `DecompressionStream('gzip')` on the cut file | error at the end: `TypeError: Compressed input was truncated.` | `TypeError: The input is ended without reaching the stream end` | `TypeError: Incomplete compressed input.` |
| data before that error | mate 1: 25,328,351 bytes, 109,190 complete records; mate 2: 25,195,052 bytes, 108,620 records | same | same |
| fflate (`gz.js` `VerifyingGunzip`), push without end | no error; same bytes and records as above | same | same |
| fflate, end of stream (`push(..., true)`) | `Error: unexpected EOF` | same | same |
| kallisto on the two cut `.gz` files (same byte count) | exit 0, no error or warning in stderr; `n_processed` 108,620 (the smaller record count of the two mates) | same | same |
| kallisto, mate 2 cut at 6,291,456 bytes (81,481 records), mate 1 at 8,388,608 (109,190 records) | exit 0, no warning; `n_processed` 81,481 | same | same |
| N-pair method: `keepFirstReads` to 50,000 records on each mate, stop | pass: 50,000 and 50,000 records; read names equal record by record for all 50,000 pairs | pass | pass |
| kallisto on those 50,000 pairs | pass: `n_processed` 50,000, `n_pseudoaligned` 47,100 | pass (see the Firefox output fault) | pass |
| same cut with `DecompressionStream` instead of fflate | same plain bytes (SHA-256 equal) | same | same |
| stream the whole URL without `Range`, stop at 50,000 records, cancel (run only in the 19 sessions where both ranges came) | pass 7 of 7 | pass 8 of 9; 1 fail: `TypeError: Error in input stream` (run 44986429632, session 1) | pass 2 of 3; 1 fail: `TypeError: Load failed` (run 44985823017, session 2) |
| ask for more records than the cut file holds (218,380) | `Error: unexpected EOF` | same | same |
| ENA filereport from the page | HTTP 200, `read_count` 8405534 | same | same |
| GitHub release asset `yeast_sub_1.fastq.gz` of `stajichlab/genomics-wasm` (run 44990656616 only) | fail x3: `TypeError: Failed to fetch` | fail x3: `TypeError: NetworkError when attempting to fetch resource.` | fail x3: `TypeError: Load failed` |

- Pairing: kallisto pairs records by position. It gives no error when one mate has fewer records; it stops at the shorter mate.
  So two cut files with different record counts give a silent short run. The page must cut both mates at the same record count itself.
- Method that gives exactly N pairs from two cut streams (tested): read each mate stream through a gzip decoder, count newlines to N x 4,
  keep the bytes up to that newline, cancel the stream. Do not call the decoder end on a cut stream (it throws). Check that both mates have N records.
  The method works with fflate and with `DecompressionStream`, from a Range request or from a plain request that is cancelled.
  Overall, the stream method (no Range) succeeded in 17 of 30 sessions: 17 of the 19 sessions where ENA answered the Range requests.
- If the stream ends before N records, fflate throws `unexpected EOF` at the end. The page must show this as "the file has fewer than N pairs" or as a damaged download.
- GitHub release assets: the request redirects to `release-assets.githubusercontent.com`. WebKit logged
  `Cross-origin redirection to https://release-assets.githubusercontent.com/... due to access control checks.` and, for the Range request,
  `Origin http://127.0.0.1:PORT is not allowed by Access-Control-Allow-Origin. Status code: 403`. So a page could not read that release asset in these
  three browsers (tested from an `http://127.0.0.1` origin only). The tested asset belongs to `stajichlab/genomics-wasm`
  (`measure.mjs`, `ghAsset`), not to this repository. That release assets of this repository behave the same is an inference (same GitHub download
  path and redirect host), not a measurement. On that basis the fallback in the spec (the test-set button fetches release assets) does not work as written.

ENA availability from the runners (each session needs both mates):

| run | Chromium sessions with both mates fetched | Firefox | WebKit |
|---|---|---|---|
| 44984897294 | 1 of 1 | 1 of 1 | 1 of 1 |
| 44985823017 | 2 of 3 | 3 of 3 | 2 of 3 |
| 44986429632 | 1 of 3 | 2 of 3 | 0 of 3 |
| 44990656616 | 3 of 3 | 3 of 3 | 0 of 3 |
| total | 7 of 10 | 9 of 10 | 3 of 10 |

- In all, 19 of 30 sessions fetched both 8 MiB ranges. Errors: Chromium `TypeError: Failed to fetch` (log: `net::ERR_FAILED`, `net::ERR_CONTENT_LENGTH_MISMATCH`);
  Firefox `TypeError: Content-Length header of network response exceeds response Body.` and `TypeError: Error in input stream`;
  WebKit `TypeError: Load failed`, with `Status code: 403` and no CORS header in the log of run 44990656616.
- Speed of a successful 8 MiB range: 29.2 to 1,813 kB/s. The stream test read about 4 MB at 152 to 6,452 kB/s.
- In CI, a full download of the ENA FASTQ with `curl` failed (run 44984006776). The NCBI SRA file (1,082,915,326 bytes, from the AWS open data URL of the NCBI locate API) downloaded in
  16.9 s (run 44984897294), 43.0 s (run 44985823017), 15.5 s (run 44986429632) and 19.6 s (run 44990656616).
- Not tested: ENA from a user's network or from an `https` page origin; retries; other ENA hosts.

### Test-set fetch method (decision)

- Method: fetch each mate from ENA, decompress as a stream, cut at N records on both mates, cancel the stream. Check both mates have N records
  and that the read names match record by record. Then give kallisto plain FASTQ (see "Decisions").
- Bytes to read from ENA for N = 200,000 pairs: about 15.4 MB per mate (the 8,388,608-byte cut held 109,190 records of mate 1, 76.8 bytes per record; computed, not measured).
- The page must show a clear error when ENA fails (the errors above), and offer a retry. A retry was not tested.
- Fallback: GitHub release assets cannot be read from the page (CORS, all 9 sessions). A fallback that serves the subset files from the same origin as the page
  (for example, files added to the GitHub Pages site at build time) was not tested: UNVERIFIED. The spec "To measure" item 2 now points here.
- Browser tests of the test sets in CI should serve the subset files from the local test server, not from ENA, because ENA failed in 11 of 30 sessions.

## Run-to-run output of kallisto (determinism suite, and a Firefox output fault)

Suite `determinism`, run 44986429632: per engine, 3 sessions; each session builds one index, prepares the first 100,000 pairs as `gz` and as `plain`,
and runs quant twice per input with that same index. The index and every `abundance.tsv` were kept (artifact only).

- The index bytes differed in every session (12 of 12 distinct SHA-256 values), as in spike 0.
- Chromium, WebKit and Node: all 36 `abundance.tsv` files were identical (SHA-256 `eabe3c98...`), for both inputs and with every index.
- Firefox: 10 of 12 files equal that value. Two files (session 1, `plain`, quant 1 and quant 2; same index and same input bytes) differ by one line each,
  and the change is a repeated number token (`check_abundance.py`):
  - line 6119: `rna-NM_001183914.1  1491  1327.5  9  31.1145` became `... 999  31.1145` (`est_counts` 9 written as 999; the sum of `est_counts` became 95,268.997 instead of 94,278.997);
  - line 3167: `... 40  204.998` became `... 40  204.998204.998204.998` (the line no longer parses as numbers).
- The same kind of difference was seen earlier: run 44984897294, Firefox `gz` at 100,000 pairs (SHA-256 `46b1209c...` against `eabe3c98...`), and
  the truncated test in that run (Firefox, 50,000 pairs, `d04e452f...` against `8e9924b4...`). Those files were not kept.
- Count over all kallisto outputs that were compared (runs 44984897294 to 44990656616): Firefox 4 damaged of 130 (quant 1 of 4 in run 44984897294,
  0 of 87 in runs 44985823017 and 44990656616; truncated tests 1 of 27; determinism 2 of 12). Chromium, WebKit and Node: 0 of 269.
- Cause: UNVERIFIED. The fault is in the written text, not in the counts (other lines and `n_pseudoaligned` are equal).
  It was not seen in Chromium, WebKit or Node.
- Effect on the page: an `abundance.tsv` from Firefox can be wrong without an error. A check that would have caught both cases here:
  every data line parses as 4 numbers, the row count equals the transcript count, and the sum of `est_counts` equals `n_pseudoaligned` within rounding
  (94,278.997 against 94,279). This check was computed on these files only; a page that checks and repeats the run was not tested.
- Effect on the acceptance test (Task 9): a byte-identical comparison of Firefox output with the CLI can fail without a code fault
  (4 of 130 Firefox kallisto runs here, about 3%; cause UNVERIFIED). Task 9 must re-run the sample or check the file before it reports a mismatch.
  Needed page behaviour: check each `abundance.tsv` as above, and re-run the sample once when the check fails (untested).

## Genomes: download size, gunzip and transcript building

Suite `genomes`, run 44985823017, 3 runs per genome and engine. Files from the catalog URLs (NCBI FTP over https).
Each run: gunzip of the `.fna.gz` and `.gff.gz` with fflate (`VerifyingGunzip`) and with `DecompressionStream` (browsers) or the same in Node,
`build_transcriptome` of `qrnamap-wasm` (web package in browsers, Node package in Node), then `kallisto index` of the transcripts.

| genome | accession | why | genome size (catalog) | `.fna.gz` bytes | `.gff.gz` bytes | FASTA bytes after gunzip | genes | transcripts | transcript FASTA bytes | CI download s (fna, gff) |
|---|---|---|---|---|---|---|---|---|---|---|
| Saccharomyces cerevisiae S288C | GCF_000146045.2 | named | 12,071,326 | 3,843,460 | 2,203,839 | 12,310,392 | 6,477 | 6,478 | 9,226,153 | 0.88, 0.59 |
| Aspergillus fumigatus Af293 | GCF_000002655.1 | named | 29,384,958 | 9,003,172 | 1,984,609 | 29,752,955 | 10,086 | 10,086 | 19,061,020 | 0.72, 0.59 |
| Cryptococcus neoformans H99 | GCF_000149245.1 | named | 18,891,193 | 5,942,468 | 1,996,629 | 19,153,849 | 8,338 | 9,189 | 18,154,107 | 0.71, 0.58 |
| Fulvia fulva Race5_Kim | GCF_020509005.1 | 3 largest | 67,169,167 | 20,677,971 | 2,527,762 | 68,009,607 | 14,993 | 14,993 | 19,661,494 | 1.23, 0.62 |
| Laccaria bicolor S238N-H82 | GCF_000143565.1 | 3 largest | 64,877,444 | 18,136,286 | 3,795,369 | 65,761,131 | 18,261 | 18,261 | 21,626,169 | 0.82, 0.66 |
| Colletotrichum lupini IMI 504893 | GCF_023278565.1 | 3 largest | 63,407,421 | 19,663,029 | 3,250,524 | 64,237,885 | 18,681 | 18,681 | 26,301,685 | 0.79, 0.66 |

Times in seconds, 3 runs. Peaks: whole-browser tree Pss (browsers) or sampled RSS (Node), MiB.

| genome | engine | gunzip fflate s | gunzip DecompressionStream s | build_transcriptome s | kallisto index s | index WASM heap MiB | peak during gunzip and build | peak during index | page fetch of both files from NCBI s (run 1) |
|---|---|---|---|---|---|---|---|---|---|
| S. cerevisiae | node | 0.25, 0.25, 0.25 | 0.11, 0.14, 0.12 | 0.20, 0.20, 0.20 | 21.16, 21.20, 21.46 | 117.1 | 239.7, 240.6, 252.8 | 421.2, 437.9, 450.8 | - |
| S. cerevisiae | chromium | 0.26, 0.25, 0.27 | 0.08, 0.08, 0.08 | 0.19, 0.20, 0.20 | 24.64, 24.27, 24.30 | 117.1 | 332.0, 325.8, 326.3 | 578.5, 511.0, 511.6 | 0.54 |
| S. cerevisiae | firefox | 0.26, 0.27, 0.26 | 0.10, 0.10, 0.09 | 0.20, 0.20, 0.20 | 38.16, 37.76, 38.61 | 117.1 | 586.5, 520.7, 508.4 | 805.9, 804.5, 802.9 | 0.38 |
| S. cerevisiae | webkit | 0.37, 0.37, 0.37 | 0.12, 0.10, 0.10 | 0.17, 0.17, 0.18 | 22.22, 22.17, 21.54 | 117.1 | 499.9, 431.6, 521.9 | 858.4, 837.0, 835.6 | 1.43 |
| A. fumigatus | node | 0.44, 0.44, 0.44 | 0.26, 0.26, 0.28 | 0.60, 0.60, 0.60 | 48.68, 48.19, 48.31 | 204.4 | 493.1, 494.4, 496.8 | 755.5, 756.0, 751.0 | - |
| A. fumigatus | chromium | 0.55, 0.51, 0.52 | 0.17, 0.17, 0.17 | 0.47, 0.46, 0.47 | 55.41, 53.82, 53.94 | 204.4 | 532.6, 523.2, 505.6 | 830.4, 770.7, 772.4 | 0.70 |
| A. fumigatus | firefox | 0.52, 0.52, 0.50 | 0.21, 0.20, 0.20 | 0.50, 0.49, 0.49 | 85.31, 84.20, 86.06 | 204.4 | 832.3, 757.0, 754.1 | 1006.7, 992.6, 989.0 | 0.35 |
| A. fumigatus | webkit | 0.78, 0.76, 0.78 | 0.20, 0.19, 0.20 | 0.38, 0.36, 0.38 | 50.71, 49.58, 50.03 | 204.4 | 724.2, 726.3, 729.9 | 1418.0, 1392.8, 1392.1 | 1.08 |
| C. neoformans | node | 0.43, 0.44, 0.43 | 0.27, 0.26, 0.26 | 0.68, 0.67, 0.67 | 38.69, 38.83, 39.10 | 212.8 | 425.4, 437.7, 446.8 | 746.3, 764.6, 765.0 | - |
| C. neoformans | chromium | 0.50, 0.51, 0.49 | 0.16, 0.15, 0.15 | 0.53, 0.53, 0.53 | 43.59, 43.33, 43.61 | 212.8 | 615.6, 563.9, 582.4 | 804.7, 782.1, 781.8 | 0.70 |
| C. neoformans | firefox | 0.52, 0.53, 0.47 | 0.21, 0.19, 0.19 | 0.59, 0.57, 0.56 | 69.97, 68.91, 67.89 | 212.8 | 855.1, 797.6, 822.2 | 1026.5, 1018.6, 1018.6 | 0.43 |
| C. neoformans | webkit | 0.72, 0.71, 0.74 | 0.20, 0.19, 0.19 | 0.42, 0.41, 0.43 | 38.71, 39.79, 39.76 | 212.8 | 780.4, 718.4, 733.7 | 1444.0, 1424.3, 1422.8 | 1.06 |
| F. fulva | node | 0.81, 0.80, 0.81 | 0.50, 0.47, 0.49 | 0.81, 0.81, 0.81 | 48.49, 48.03, 48.26 | 248.6 | 723.2, 725.4, 663.5 | 986.8, 986.5, 1001.2 | - |
| F. fulva | chromium | 1.00, 1.00, 0.97 | 0.31, 0.32, 0.32 | 0.67, 0.64, 0.72 | 54.10, 54.33, 54.10 | 248.6 | 766.8, 691.7, 691.7 | 980.4, 974.8, 974.5 | 3.17 |
| F. fulva | firefox | 1.08, 0.99, 0.97 | 0.43, 0.39, 0.41 | 0.65, 0.66, 0.63 | 83.00, 84.00, 83.44 | 248.6 | 1001.4, 983.1, 991.5 | 1166.8, 1132.3, 1132.9 | 0.87 |
| F. fulva | webkit | 1.46, 1.44, 1.44 | 0.45, 0.42, 0.42 | 0.47, 0.48, 0.49 | 48.56, 47.51, 48.74 | 248.6 | 976.9, 996.0, 1034.2 | 1639.7, 1728.9, 1635.4 | 1.81 |
| L. bicolor | node | 0.91, 0.91, 0.91 | 0.60, 0.60, 0.59 | 1.33, 1.32, 1.33 | 49.12, 49.70, 49.80 | 302.9 | 835.5, 812.0, 835.5 | 1099.3, 953.7, 972.1 | - |
| L. bicolor | chromium | 1.15, 1.13, 1.13 | 0.36, 0.37, 0.37 | 1.02, 1.00, 1.01 | 57.62, 58.18, 57.80 | 302.9 | 971.2, 876.9, 970.0 | 1114.2, 1171.6, 1100.5 | 0.79 |
| L. bicolor | firefox | 1.13, 1.15, 1.11 | 0.43, 0.42, 0.44 | 1.11, 1.08, 1.08 | 87.81, 87.22, 89.38 | 302.9 | 1227.7, 1210.5, 1213.4 | 1462.4, 1385.6, 1430.6 | 0.51 |
| L. bicolor | webkit | 1.70, 1.70, 1.67 | 0.42, 0.46, 0.43 | 0.79, 0.79, 0.80 | 50.40, 49.73, 51.47 | 302.9 | 1016.8, 976.1, 1008.6 | 1796.8, 1785.0, 1765.6 | 1.24 |
| C. lupini | node | 0.84, 0.85, 0.84 | 0.54, 0.53, 0.48 | 1.08, 1.09, 1.09 | 67.23, 65.83, 65.22 | 299.0 | 729.9, 755.2, 755.5 | 986.2, 982.1, 986.3 | - |
| C. lupini | chromium | 1.02, 1.02, 1.02 | 0.36, 0.33, 0.33 | 0.85, 0.86, 0.85 | 74.12, 73.59, 74.22 | 299.0, 358.8 | 822.7, 827.5, 871.8 | 1148.8, 1055.4, 1056.7 | 0.25 |
| C. lupini | firefox | 1.02, 0.99, 1.02 | 0.39, 0.40, 0.39 | 0.93, 0.88, 0.87 | 113.57, 112.64, 119.46 | 299.0, 358.8 | 1063.4, 1083.4, 1115.0 | 1310.6, 1337.1, 1347.9 | 0.60 |
| C. lupini | webkit | 1.53, 1.56, 1.57 | 0.42, 0.43, 0.46 | 0.66, 0.66, 0.69 | 64.29, 68.33, 70.75 | 299.0 | 937.7, 973.9, 934.7 | 1713.2, 1684.9, 1702.8 | 1.38 |

- All 72 runs finished (6 genomes x 4 engines x 3 runs). The transcript FASTA had one SHA-256 per genome over all engines and runs.
- The fflate and `DecompressionStream` outputs had the same SHA-256 in every run. Every NCBI file here is one gzip member (`VerifyingGunzip` count 1 for each file).
  `DecompressionStream` was 2.37 to 4.05 times faster than fflate in the browsers and 1.52 to 2.27 times in Node
  (in Node the `DecompressionStream` column is the Node built-in class).
- The page fetch from NCBI (run 1 of each browser) gave HTTP 200 and the same SHA-256 as the CI copy for all 6 genomes in all 3 browsers.
- Largest values: index 119.46 s (*C. lupini*, Firefox); peak 1,796.8 MiB during the index (*L. bicolor*, WebKit); index heap 358.8 MiB
  (*C. lupini*, one run in Chromium and one in Firefox; the other runs 299.0 MiB). All are below 2 GiB.
- The catalog sets `max_genome_size` 70,000,000 and a size warning above 50,000,000 bases. The three largest catalog genomes (63.4 to 67.2 Mb) finished
  in all engines. No catalog genome was too large in this test. A genome size limit from these data: none needed below 67.2 Mb on this runner.
- Not tested: genomes uploaded by users above 67.2 Mb; devices with less memory.

## Page genome path: peak memory

- The genome step of the page (Task 5) does more than the harness above: the page downloads both files, gives them to a page worker,
  and the page worker starts the kallisto worker (a nested worker) through `ToolRunner`. The harness ran `build_transcriptome` in the page
  and started the kallisto worker from the page.
- Peak memory and time of this page path: UNVERIFIED. No run measured them.
- Job that measures them: `.github/workflows/page-genome-memory.yml` (dispatch only). It loads the real page in Chromium, Firefox and WebKit,
  chooses Fulvia fulva Race5_Kim (GCF_020509005.1, the page size limit) in the catalog menu, downloads it from NCBI, and records the time to
  "index ready" and the whole-browser tree Pss peak (sampler copied from `measure.mjs`), 3 runs per engine (`validation/browser-page/tests/genome-memory.mjs`).
- Values to compare with (harness, section "Genomes"): index peak 1,728.9 MiB for F. fulva (WebKit, the highest of its 3 runs);
  the highest index peak of all genomes was 1,796.8 MiB (L. bicolor, WebKit).

## Page quantification path: peak memory

- The samples step of the page (Task 6) differs from the harness above. The page worker reads each local file twice (pass 1 counts the bytes of the first N records,
  pass 2 copies them into one array of that exact size), with `gz.js` (fflate) for gzip files. The page worker then starts the kallisto worker (a nested worker)
  through `ToolRunner`, with the index and the plain reads. The harness prepared the reads in the page or in a dedicated worker and started kallisto from the page.
- First run: workflow run 45049907640, N = 200,000 pairs for every case. Quant-phase whole-browser Pss, as reported to this task (the artifact TSVs
  were not read in this task, so these numbers are not checked here): H99 150-base plain: Chromium 1,204 MiB, Firefox 1,663 to 1,759 MiB, WebKit 1,993 to 2,044 MiB;
  yeast 101-base: WebKit up to 1,508 MiB. The 150-base values are above the N rule (1.5 GiB in all three browsers). So N must depend on the read length.
- User decision (2026-11-22): the preview size N depends on the read length. The page takes N from `previewLimit(readLength)`
  (`site/analyze/preview-limit.js`; read length = the longer first record of the two mates). Table now: reads of up to 101 bases: 200,000 pairs
  (measured, Task 3 and this section); longer reads: 75,000 pairs (measured with the prep worker, see "Results" below). The page field "Read pairs per sample" can only lower N.
- N rule: all 3 runs finish, at most 300 s per sample, and a whole-browser Pss of at most 1.5 GiB (1,572,864 kB) in all three browsers.

### Results (CI runs; highest value over 3 runs; tree Pss in kB; times all under 15 s per sample)

Yeast 101-base reads, N = 200,000 (run 45063764607): Chromium 861,919, Firefox 1,366,284, WebKit 1,483,056 (quant phase). The rule is met.
(Run 45099586699 used the prep worker, see below; yeast has not been re-measured with it yet.)

H99 150-base reads, default path (transcripts built in the page worker), WebKit:

| N | quant phase, run 45063764607 | quant phase, run 45094297581 | highest of all phases, run 45094297581 |
|---|---|---|---|
| 50,000 | not run | 1,552,121 | 1,573,977 |
| 75,000 | not run | 1,571,705 | 1,593,941 |
| 100,000 | 1,611,032 | 1,655,965 | 1,655,965 |
| 125,000 | 1,752,310 | not run | not run |
| 150,000 | 1,854,361 | not run | not run |

The same N gave 1,611,032 and 1,655,965 kB in two runs, so the run-to-run difference is about 50 MB. Chromium and Firefox met the rule at all these N
(Firefox 1,492,485 at N = 150,000; Chromium 1,024,984). In WebKit the default path met the rule in no case when all phases count.

Prep worker (run 45099586699, commit 6d638e6; the page reads the genome files and builds the transcripts in a worker that ends before the kallisto index starts;
page option `?prep_worker=0` turns it off). Highest of all phases (genome, index, prep, quant), 3 runs each:

| H99 150-base | Chromium | Firefox | WebKit |
|---|---|---|---|
| default path, N = 50,000 | 797,884 | 1,148,810 | 1,586,569 |
| prep worker, N = 50,000 | 645,805 | 989,096 | 1,491,448 |
| prep worker, N = 75,000 | 666,758 | 1,094,731 | 1,467,306 |

Genome-step time was the same with and without the prep worker (Chromium about 41 s, Firefox 118 to 124 s, WebKit 59 to 63 s).
The prep worker is the page default.

Final run 45111397838 (commit 2fa8411, prep worker on, 3 runs per case, all finished, at most 16 s per sample). Highest of all phases (kB):

| Case | Chromium | Firefox | WebKit |
|---|---|---|---|
| yeast, 101 bases, N = 200,000, gzip | 781,325 | 1,250,067 | 1,450,575 |
| yeast, 101 bases, N = 200,000, plain | 779,092 | 1,255,986 | 1,424,803 |
| H99, 150 bases, N = 75,000 | 670,358 | 1,080,402 | 1,423,276 |
| H99, 150 bases, N = 100,000 | 769,371 | 1,167,731 | 1,498,957 |
| H99, 150 bases, N = 75,000, `?prep_worker=0` | 812,133 | 1,285,185 | 1,575,591 |

The rule (1,572,864 kB) is met in all three browsers by every case with the prep worker. Without it, WebKit at N = 75,000 is over the rule by 2,727 kB.
Decision: N = 75,000 for reads longer than 101 bases. Run 45111397838 met the rule at N = 100,000 (peak 1,498,957 kB), but the margin was small
(73,907 kB, 4.7%). N = 75,000 was chosen for safety (peak 1,423,276 kB, margin 149,588 kB, 9.5%). N = 125,000 with the prep worker is UNVERIFIED.
Reads longer than 150 bases and other genomes are not measured.

- Job that measures them: `.github/workflows/page-quant-memory.yml` (script `validation/browser-page/tests/quant-memory.mjs`, reads from
  `validation/browser-page/tests/quant-memory-reads.sh`). It starts by `workflow_dispatch` (only when the file is on the default branch), or by a push to branch
  `page-v1` that changes `validation/browser-page/memory-run.trigger` or the workflow file. It loads the real page in Chromium, Firefox and WebKit, 3 runs per engine and case.
  Per run: the catalog genome from NCBI through the menu (index built in the page), then one sample through the file chooser, counted at N pairs.
  N for yeast: the page value for 101-base reads (200,000). N for H99: 100,000, 125,000 and 150,000, set with the page option `?preview_pairs=N`
  (it replaces the read-length limit; the reads go through the same page path). The TSV has the columns `n_option`, `read_length` and `n` (the N the page used).
  It records the whole-browser tree Pss peak per phase (genome, prep, quant; sampler copied from `measure.mjs`), the time of each phase,
  the pairs, n_processed and n_pseudoaligned, and the re-runs.
- Cases: yeast SRR3396381, the first 1,000,000 pairs (release asset, SHA-256 checked against `measure/reads_prefix.tsv`), given as gzip and as plain FASTQ,
  with S. cerevisiae S288C (GCF_000146045.2); C. neoformans H99 SRR23037513 (BioProject PRJNA922535), the first 200,000 pairs from NCBI SRA (`fastq-dump -X 200000`),
  plain FASTQ, with H99 (GCF_000149245.1), at the three candidate N. The H99 cases are left out of a run when the SRA download fails on the runner.
- Read lengths: 101 bases for SRR3396381 (both mates, all 1,000,000 records of the prefix files). For SRR23037513 the first 200,000 records of each mate are
  150 bases long (counted on the cluster on 2026-11-22; the job writes the counts again in `reads_check.tsv`). Before run 45049907640, reads longer than
  101 bases were not measured: every value in the sections above uses SRR3396381.
- Values to compare with (mode `plain`, 200,000 pairs, run 44985823017, `results/quant_runs.tsv`, sample peak = the higher of the preparation and quant
  tree Pss peaks): Chromium 722 to 753 MiB, Firefox 1,105 to 1,109 MiB, WebKit 1,234 to 1,327 MiB.
- A script check (not a measurement for this document): 1 Chromium run per case on the cluster node (not a runner), 2026-11-22. All finished with 200,000 pairs
  processed; n_pseudoaligned was 188,616 for yeast gzip and yeast plain (equal to the harness value at 200,000 pairs) and 192,174 for H99 plain.

## Read counts and file sizes (ENA portal API, 2026-11-22)

Query: `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=<run or PRJNA976751>&result=read_run&fields=...&format=tsv`.
Full rows: `results/yeast_gse80357_ena.tsv`, `results/cneo_prjna976751_ena.tsv`. All runs are paired-end.

Yeast amphotericin B set (GSE80357; run accessions from `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb/samples.tsv`), Illumina HiSeq 2000:

| sample | run | ENA sample title | read pairs | `_1.fastq.gz` bytes | `_2.fastq.gz` bytes |
|---|---|---|---|---|---|
| Control_1 | SRR3396385 | Matched control for AMB only rep 1 | 9,625,471 | 765,392,684 | 776,514,673 |
| Control_2 | SRR3396386 | Matched control for AMB only rep 2 | 10,302,142 | 816,554,593 | 828,166,555 |
| Control_3 | SRR3396387 | Matched control for AMB only rep 3 | 9,997,148 | 796,139,461 | 808,975,406 |
| AmpB_1 | SRR3396381 | AMB only rep 1 | 8,405,534 | 669,296,728 | 680,167,572 |
| AmpB_2 | SRR3396382 | AMB only rep 2 | 9,755,935 | 781,154,812 | 791,042,140 |
| AmpB_3 | SRR3396384 | AMB only rep 3 | 10,175,785 | 815,031,018 | 822,459,056 |

*C. neoformans* H99 Lkh1 set (PRJNA976751), basal condition, 3 against 3 (the ENA sample titles name the condition):

| group | run | ENA sample title | read pairs | `_1.fastq.gz` bytes | `_2.fastq.gz` bytes |
|---|---|---|---|---|---|
| wild type | SRR24750781 | H99 basal sample_1 | 18,134,154 | 597,146,669 | 597,682,343 |
| wild type | SRR24750780 | H99 basal sample_2 | 15,089,131 | 506,584,759 | 508,140,072 |
| wild type | SRR24750777 | H99 basal sample_3 | 18,176,463 | 601,705,504 | 603,814,811 |
| *lkh1* deletion | SRR24750779 | lkh1 mutant basal sample_1 | 16,437,474 | 547,724,942 | 550,950,870 |
| *lkh1* deletion | SRR24750774 | lkh1 mutant basal sample_2 | 18,221,839 | 593,288,172 | 594,487,222 |
| *lkh1* deletion | SRR24750773 | lkh1 mutant basal sample_3 | 17,635,731 | 585,386,266 | 586,984,165 |

The project has 6 more runs (rapamycin 60 min, wild type and *lkh1*); they are not in the test set.

## *C. neoformans* gene identifiers against RefSeq H99 (GCF_000149245.1)

- Supplementary Table S3 of the paper (Front Cell Infect Microbiol 2024, PMID 38774630, PMC11106425) was fetched on 2026-11-22 from the Europe PMC
  supplementary files service (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11106425/supplementaryFiles`, zip SHA-256 `58bac5b4...c2814`,
  `Table_3.xlsx` SHA-256 `7b723ed2...b44b8d7`). The PMC download link returned an HTML page ("Preparing to download"), not the file.
  The files are on `/bigdata` (`qRNAmap_data/page-v1/dl_s3/`), not in git.
- Sheet `lkhl_basal_vs_H99_basal`: 409 genes (276 with a positive fold change, 133 negative). This equals the 276 up and 133 down in the spec.
  Columns include `Gene_ID` (NCBI GeneID), `locus_tag` (`CNAG_xxxxx`) and `Protein_ID` (`XP_...`).
- `s3_vs_gff.py` against the RefSeq H99 GFF3 (catalog URL, SHA-256 `69ee1730...19cf3`): all 409 locus tags are genes in the GFF3, and all 409
  (`locus_tag`, `GeneID`) pairs equal the GFF3 pairs. 293 of 409 `Protein_ID` values are CDS `protein_id` values in the GFF3
  (S3 `gene_biotype`: 326 protein_coding, 81 ncRNA, 1 rRNA, 1 tRNA; 33 `Protein_ID` cells list several proteins; `results/cneo_s3_ids.tsv`).
- The gene names that `qrnamap transcripts` makes from this GFF3 (CLI, debug build of this branch, run on the cluster) are the locus tags:
  8,338 genes, 9,189 transcripts. All 409 S3 locus tags are among them. So the paper's identifiers match the RefSeq H99 annotation and the qRNAmap gene names.
- The S3 count columns are named `H99_basal_3`, `_6`, `_7` and `lkhl_basal_3`, `_6`, `_7`. The ENA titles say `sample_1` to `sample_3`.
  Which ENA run is which S3 column is UNVERIFIED.

## Failed, UNVERIFIED and not tested

Failed (measured):
- Run 44984006776: full ENA FASTQ download from the runner (curl error 18). Fixed by the NCBI SRA source in CI.
- ENA from the browsers: 11 of 30 sessions did not get both 8 MiB ranges (errors in the truncated section).
- GitHub release asset from the page: blocked in Chromium, Firefox and WebKit (9 of 9 sessions).
- Chromium, 4,000,000 plain pairs (`plain`, `wplain`): `DataCloneError ... out of memory`.
- WebKit: 3,000,000 pairs (`plain`, `wplain`) and 500,000 pairs (`wgz`): timeout after 30 min; 3,000,000 pairs (`gz`): 1 run stopped at 10 GiB.
- Firefox wrote a damaged `abundance.tsv` in 4 of 130 kallisto runs.

UNVERIFIED (not measured, or cause not known):
- The cause of the Firefox output fault.
- The cause of the different kallisto index bytes in every run (as in spike 0); the quant output was equal.
- Which ENA run belongs to which S3 count column (`_3`, `_6`, `_7`) of the *C. neoformans* table.
- A same-origin fallback for the test-set files; retries of ENA requests; ENA from user networks.
- Memory and time on user devices, under memory pressure, or with less than 16 GB; release builds of the browsers; an `https` page origin;
  the real GitHub Pages host.
- Phone-sized screens (spec "To measure" item 5): not part of this task, not measured.

Not tested:
- Streaming reads into kallisto (the vendored runtime copies whole files into MEMFS); the runtime `stdin` option; a way to avoid the three copies.
- Counts between the tested counts (for example 300,000 pairs).
- Single-end reads; other datasets than SRR3396381; reads longer than 101 bases.
- User-uploaded genomes above 67.2 Mb.
- Peaks shorter than 250 ms can be missed by the Pss and RSS sampling.

## Files

- `validation/browser-page/measure/`: `measure.mjs` (harness), `steps.mjs` (steps shared by page and Node), `prep_worker.mjs`, `heap_wrap.mjs`,
  `qm/` (qmutfind `gz.js`, `reads.js`, `vendor/fflate.js` and its licence, copied without change from qmutfind commit
  `4404e1c93f47ea527287b93d6b496070563b193a`, one header line added), `prepare_reads.sh`, `prepare_genomes.py`, `summarize.py`,
  `check_abundance.py`, `s3_vs_gff.py`, `package.json` and `package-lock.json` (Playwright 1.64.0).
- `validation/browser-page/measure/results/`: the merged TSVs of the runs above.
- Full artifacts (logs, indexes, abundance files): the workflow runs (30-day retention) and `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/page-v1/runs/`.
