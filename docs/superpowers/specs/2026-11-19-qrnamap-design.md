# qRNAmap: design

Date: 2026-11-19. Status: draft for approval. Target repository: `stajichlab/qRNAmap` (not created yet).

## 1. Intent

qRNAmap maps bulk RNA-seq reads and tests for differential expression. It runs as a CLI and as a web page (WebAssembly, WASM) from one code base. It follows the design of `../qmutfind`.

Two mapping modes:
- **Fast mode:** pseudoalignment of reads to a transcript set (kallisto first, Salmon second; the user picks one per run).
- **Genome mode:** spliced alignment to the genome (HISAT2), then read counting per feature (featureCounts).

Both modes end in a gene-by-sample count table. DESeq2 then tests a model the user chooses.

Why: the lab works on fungal genomes. A user should go from a genome choice, a sample sheet and FASTQ files to counts and a DE table, in a browser or on the cluster, with the same steps and the same results.

## 2. Goals, non-goals, anti-goals

Goals
- One sample sheet drives the whole pipeline in the CLI and in the browser.
- The genome is loaded as in qmutfind (menu, NCBI download, user file, `.gz` and `.zst`).
- The transcript FASTA is built from genome FASTA plus GFF3, on the fly. A downloaded transcript file is also accepted.
- Outputs: per-sample abundance, a gene count matrix, mapping summaries, DESeq2 results, normalized counts, QC plots.
- Native CLI and WASM give the same counts for the same small input. A test checks this.

Non-goals (version 1)
- No new aligner or quantifier algorithm. We run the upstream tools.
- No single-cell, no long-read RNA, no novel transcript assembly (StringTie is later work).
- No human-scale genomes in the browser. The CLI handles larger genomes.
- No server. All browser computation is local.

Anti-goals
- A browser result that differs from the CLI result without a test that shows and explains it.
- An API key in browser code.
- Reporting speed or accuracy numbers that we did not measure.

## 3. Constraints and facts (verified 2026-11-19 unless marked)

- `comparativechrono/genomics-wasm` builds upstream tools for WASM: kallisto 0.52.0 (newest upstream release is v0.52.0), HISAT2 2.2.2 (upstream is 2.2.3), featureCounts 2.1.1, Salmon 2.8.0, minimap2 2.31, StringTie 3.0.3, DESeq2 1.52.0 (webR). Repository facts: 0 stars, GitHub shows no licence (its own `LICENSE` is said to be MIT, see qmutfind PLAN section 4), last push 2026-10-24. I have not built or run any of these tools.
- Each tool call runs in a fresh worker with a fresh virtual filesystem (`ToolRunner`). Index files are passed in again for every call. Heap limit 2 GiB. All are single-threaded except featureCounts (1 to 4 threads) and Salmon.
- featureCounts needs shared memory, so the page needs the headers `Cross-Origin-Opener-Policy: same-origin` and `Cross-Origin-Embedder-Policy: require-corp`. kallisto and HISAT2 do not need them. UNVERIFIED: whether these headers break the genome downloads from NCBI and reads from ENA (cross-origin isolation blocks resources that lack CORP or CORS headers).
- kallisto in WASM uses Emscripten `MEMORY64=2`, which Emscripten marks experimental. Bootstrap tables differ between native and WASM. The main abundances matched in the package's checks.
- HISAT2 package says: use small references. A genome-scale index in a browser is not established. UNMEASURED for a 29 Mb fungal genome.
- The webR DESeq2 bundle does not include `tximport`, `apeglm`, `ashr` or `glmGamPoi`.
- Genome sources (from qmutfind PLAN section 12c, measured 2026-11-12): NCBI FTP files carry `Access-Control-Allow-Origin`, so the browser can download them. FungiDB needs an API key and sends no CORS headers, so only the CLI can download from it. ENA FASTQ is browser-readable. Ensembl: UNVERIFIED (no test run).
- qmutfind uses `ureq` for downloads in its CLI crate. Its core crate does no I/O and compiles for `wasm32`.
- Global lab rules: compress large text output (`zstd` for internal, `gzip` for portable); read compressed input in new scripts; use `$SCRATCH` for node-local job I/O; pixi manifests use `[workspace]`; never use `BASH_SOURCE[0]` in SLURM scripts.

## 4. Approaches considered

1. **Wrap the existing tools (chosen).** Rust workspace for genome, GFF3, sample sheet, transcript building, matrix and plots. Native CLI calls native `kallisto`, `hisat2`, `featureCounts`, `Rscript` (from pixi). The browser calls the genomics-wasm bundles. Trade-off: the browser build depends on a young third-party repository, and the "pure Rust, one build" goal of qmutfind does not hold here.
2. **Pure-Rust reimplementation** (kallisto-like pseudoaligner, `rammap` splice mode). Rejected: large effort, and results would need proof that they match the real tools.
3. **Hybrid:** approach 1 plus a Rust read counter instead of featureCounts. Deferred. Revisit if cross-origin isolation or the 4-thread cap causes trouble. Counting against a SAM stream is small, but it needs validation against featureCounts.

Decision by the user: approach 1.

## 5. Architecture

Cargo workspace. Library crates are written so that qmutfind can adopt them later.

| Crate | Role | I/O |
|---|---|---|
| `qseq-io` | Extracted from qmutfind-core: `read_fasta`, `FastqChunker`/`FastqReader`, GFF3 gene model (`Genome`, `Transcript`), gzip/zstd detection. | none, compiles for `wasm32` |
| `genome-catalog` | Owns the genome catalog schema (versioned) and the derived layer: per-genome FASTA and GFF3 URLs, file sizes, annotation facts (same fields as qmutfind `site/genomes/genomes.json`, curated menu). URL builders for NCBI, FungiDB (Ensembl only after spike 0), size checks. Holds a small built-in snapshot. The full list of fungal assemblies comes from PopSeq's `catalog/kingdom/fungi_assemblies.csv.gz` (26,826 rows on 2026-11-11, built from NCBI `assembly_summary.txt`) as upstream input. | none |
| `genome-fetch` | Native downloads: NCBI, Ensembl, FungiDB with API key (`EUPATHDB_KEY` or key file, key sent only to the key host, never printed; code taken from `qmutfind-cli/src/fetch.rs`). Writes `.gz`/`.zst`. | native only |
| `qrnamap-core` | Transcriptome builder (genome + GFF3 to transcript FASTA and `tx2gene`), sample sheet model and validation, result parsers (kallisto `abundance.tsv`, featureCounts table), count matrix, summary stats, SVG plots, DESeq2 R script generation. | none, compiles for `wasm32` |
| `qrnamap-cli` | Subcommands, native backends through `std::process`, threads, `$SCRATCH`. | native |
| `qrnamap-wasm` | `wasm-bindgen` wrapper for `qrnamap-core`. | wasm |
| `site/` | Web page. Worker orchestrates `ToolRunner` calls to the tool bundles and calls `qrnamap-wasm`. Reuses qmutfind's `gz.js`, `reads.js`, `reference-panel.js`, `persist.js`, `theme.js`. | browser |

Question answered: yes, the "get the genome" parts become a pair of crates (`genome-catalog` pure, `genome-fetch` native). The browser downloads with `fetch` in JavaScript, using the URLs from `genome-catalog`. FungiDB stays CLI-only because of the missing CORS headers.

Where the shared crates live (my call): start inside the qRNAmap workspace, as `crates/genome-catalog` and `crates/genome-fetch`, with no dependency on other qRNAmap crates, so a later move is a copy. When qmutfind is ready to adopt them, move `qseq-io`, `genome-catalog` and `genome-fetch` to their own repository (`stajichlab/genome-catalog` is free as a repo name and `genome-catalog` and `genome-fetch` are free on crates.io, checked 2026-11-19; use a Git dependency pinned by commit, no crates.io release needed) and pin by commit, as qmutfind pins `rammap`. Doing this first would slow the first working version.

A `Backend` trait in `qrnamap-core` describes a tool call as data: arguments, input files, expected output files. The CLI executes it natively. The worker executes it through `ToolRunner`. The same `Backend` description is used by both, which keeps the two paths identical in structure.

## 6. Pipeline

```
sample sheet + genome (FASTA, GFF3)
   |
   +-- build transcripts (Rust): transcripts.fa, tx2gene.tsv   [fast mode]
   |      kallisto index -> kallisto quant per sample -> sum est_counts to gene
   |
   +-- hisat2 build (with splice sites from GFF3) -> hisat2 align per sample -> featureCounts [genome mode]
   |
   v
gene count matrix (integers) + sample sheet
   |
   v
DESeq2 (R script generated from the sheet's design formula) -> results, normalized counts, PCA, volcano, MA
```

Notes
- Transcript building: for each mRNA/transcript, join exons (or CDS plus derived UTRs when exons are missing) in transcript order, reverse-complement on the minus strand. The result is checked against `gffread` on yeast and on *A. fumigatus* Af293 (identical sequences expected; any difference is reported, not hidden).
- Splice sites for HISAT2 come from the same GFF3 model (`hisat2_extract_splice_sites` equivalent in Rust), so no Python helper is needed.
- kallisto to gene counts: sum `est_counts` per gene from `tx2gene`, round to integers for DESeq2. `tximport` is not in the webR bundle. The rounding ignores per-sample transcript length changes. This is stated in the output notes. Replace later if `tximport` becomes available.
- featureCounts options come from the sample sheet (strandedness, paired-end, feature type, attribute).
- DESeq2 needs replicates. The tool refuses to run DE with fewer than 2 samples in any tested group and says why.

## 7. Sample sheet

A TSV. One row per sample.

| Column | Meaning |
|---|---|
| `sample` | unique name |
| `fastq1`, `fastq2` | files or URLs (`fastq2` empty for single-end) |
| `condition` | main group |
| other columns | any covariates, usable in the design (batch, replicate, time) |

A separate small block (CLI flags or a web form) sets: mode (`quant` or `align`), strandedness, design formula (default `~ condition`), contrast (`condition,treated,control`).
The web page edits the sheet in a table, validates it live (unique names, matching files, balanced groups, formula uses existing columns) and downloads it as TSV. The CLI reads the same TSV. Single-end reads need fragment length mean and SD for kallisto. The sheet validator asks for them.

## 8. CLI

```
qrnamap fetch      --source ncbi|ensembl|fungidb --id ...      # genome + annotation
qrnamap transcripts --reference g.fa.gz --gff g.gff3.gz --out tx.fa.gz --tx2gene tx2gene.tsv.gz
qrnamap quant      --samples sheet.tsv --transcripts tx.fa.gz --out dir     # kallisto
qrnamap align      --samples sheet.tsv --reference g.fa.gz --gff g.gff3.gz --out dir   # hisat2 + featureCounts
qrnamap matrix     --out dir                                    # gene count matrix
qrnamap de         --counts counts.tsv --samples sheet.tsv --design '~ condition' --contrast condition,treated,control
qrnamap run        --samples sheet.tsv ...                      # all steps
```

Per-sample steps are independent. On the cluster `qrnamap run` can emit a SLURM array. Each array task should hold about 1 to 1.5 hours of work (lab rule), so several samples go in one task for small data. The task size is set from a measured timing, not a guess. Intermediate large text files are compressed. Temporary files go in `$SCRATCH`. Scripts get paths as arguments and never use `BASH_SOURCE[0]`.

## 9. Web page

Same three-step feel as qmutfind. (1) Genome: reference panel with the NCBI menu, local file, or a URL; FungiDB shows the CLI command. (2) Samples: table editor for the sheet, file picker or ENA accession per sample, mode, design. (3) Results: mapping summary per sample, count table (sortable, download), PCA, sample distance, volcano and MA plots, DE table, all downloads. Progress is per sample and per step, and each can be cancelled. Results are saved in the browser as in qmutfind (`persist.js`). Limits are shown before a run: genome size, total read bytes, estimated memory. A run that would pass the limits is blocked with a reason.

Hosting: static files plus the two isolation headers, or the service-worker approach of genomics-wasm for hosts without custom headers (GitHub Pages). The page works in a reduced mode without isolation: kallisto, HISAT2 and DESeq2 run; featureCounts needs the headers. (UNVERIFIED; spike 0 tests it.)

## 10. Build order

0. **Spikes (measure before committing).** In Node and in one real browser:
   a. Get the kallisto, HISAT2, featureCounts and DESeq2 bundles (build with Emscripten 3.1.64, or use released `dist`). Run each on yeast data and compare with native tools from pixi.
   b. HISAT2 index of a 12 Mb (yeast) and a 29 Mb (Af293) genome: time and peak memory in WASM. Same for the kallisto index of the Af293 transcript set.
   c. Cross-origin isolation with NCBI and ENA downloads.
   d. Ensembl download and CORS test.
   e. Licence check of genomics-wasm and of bundling GPL tools.
   Results go in `PLAN.md`. If 0b fails for HISAT2, genome mode becomes CLI-only in version 1 and the page offers fast mode only.
1. Workspace; extract `qseq-io` and `genome-catalog` from qmutfind with their tests.
2. Transcriptome builder with the `gffread` comparison.
3. Native CLI: `transcripts`, `quant`, `align`, `matrix`, `de`. pixi environment with kallisto, hisat2, subread, R and DESeq2.
4. Sample sheet validation, shared by CLI and WASM.
5. `qrnamap-wasm` and the web worker for fast mode with kallisto, then genome mode, then DE, then the Salmon backend (CLI first, then browser).
6. Plots, saved results, hosting.
7. Validation: CLI against WASM on a small yeast set (counts identical), CLI against published kallisto/HISAT2 runs on the same data, DE on a public dataset with known results.

Each step has tests before the next one starts. Tests use real small datasets with checked accessions (ENA read counts compared, as in the lab download rules).

## 11. Decisions by the builder (my calls, open to change)

- Counting: featureCounts (decided by the user; the genomics-wasm port is used as it is). No Rust counter.
- Quantifier: both kallisto and Salmon (decided by the user). Build order: kallisto first, then Salmon behind the same `Backend` trait. `--quantifier kallisto|salmon` in the CLI and a choice on the page. Both write gene counts in the same format, so DESeq2 does not care which one ran.
- Salmon needs threads and shared memory in the browser, so it needs the isolation headers like featureCounts. Fast mode with kallisto works without them. Salmon index format and options (for example selective alignment, decoys) are not designed yet; spike 0 reads the Salmon package README and measures it.
- Shared crates start inside this workspace and move out later.
- Single-threaded WASM tools for the first release. Threads only for featureCounts.
- Bootstrap and `tximport`-style length handling: not in version 1.
- Default DE design `~ condition`, with the formula and contrast editable.

## 11a. Decisions by the user after the builder check (2026-11-19)

- **CLI first.** Version 1 commits to the CLI (build steps 1 to 4 and the CLI side of step 7). The browser build follows after spike 0 and only for tools that pass it. A tool that fails in the browser stays in the CLI.
- **Scale:** 8 to 12 samples per experiment, up to 48. Largest genome 50 Mb. 12,000 to 15,000 genes. Reads per sample were not given. Limits and tests use these numbers. The browser limit for sample count and read bytes is set from spike 0 measurements, not guessed. Up to 48 samples in the browser is not assumed to work.
- **Repository:** `stajichlab/qRNAmap` on GitHub, GitHub Actions allowed (including Emscripten builds).

- **Catalog ownership:** `genome-catalog` owns the schema and the derived layer (URLs, sizes, annotation). PopSeq's assembly table is upstream input for the full list. Python tools such as PopSeq read the published catalog file, not the Rust code. PopSeq is not changed without the user's say-so. UNCHECKED: whether PopSeq's columns fit the schema (they have no FASTA or GFF URL columns that I saw); compare before the schema is frozen.

## 11b. Calls added after the builder check (mine, open to change)

- qmutfind is the user's own code, local at `../qmutfind` (GPL-3.0-only, verified in its `Cargo.toml`). The "MIT" in section 3 is the genomics-wasm repository's own `LICENSE`, not qmutfind's. Extracting code from it is allowed. Take the current working tree and record the qmutfind commit in the first extraction commit.
- `Backend` models one tool call: arguments, input files, expected output files. Multi-step tools are two calls (`hisat2 align`, then `featureCounts`). In the CLI the SAM/BAM stays on disk (in `$SCRATCH` during a job, compressed copy kept only if asked). In the browser the bytes pass from one call to the next through the worker's memory, one sample at a time, and are dropped after counting. The page blocks a sample whose reads would pass the memory estimate.
- DESeq2: the CLI uses the same DESeq2 version as the webR bundle and the default Wald test with no LFC shrinkage, so CLI and browser tables match. Shrinkage (`apeglm`, `ashr`) is an optional CLI-only flag, labeled in the output. Volcano and MA plots use the unshrunk estimates. The design formula is checked against the sample sheet columns in Rust before R runs.
- Browsers: test Chromium and Firefox in CI, WebKit best effort. Kallisto uses `MEMORY64=2`, which the package says compiles to ordinary 32-bit WASM, so native memory64 support should not be needed (UNVERIFIED in practice; spike 0a tests it).
- Spike 0 pass criteria (proposed): the WASM output equals native output on the yeast fixtures (counts exactly equal for kallisto and featureCounts; SAM records equal for HISAT2 apart from `@PG`); the run completes in a real browser without an out-of-memory error; the 12 Mb (yeast) genome index builds in 10 minutes or less and under 2 GiB. A failure removes that tool from the browser build, not from the CLI.
- Hosting default: GitHub Pages with the service-worker headers of genomics-wasm. Vendor or CI-build of the tool bundles is decided after spike 0 from the bundle sizes and the licence check.
- Validation data: a small yeast RNA-seq set with public accessions (chosen in spike 0, read counts checked against ENA). "Known results" means agreement with native tools run on the same data, plus a DE run on a public set whose published gene list we compare. No accession is chosen yet.

## 12. Open items

- Result of spike 0 (may remove genome mode from the browser).
- Ensembl source (needs the CORS test).
- How to ship the tool bundles: vendor `dist` files in the repository or build in CI. Size and licence decide.
- Whether HISAT2 2.2.3 (upstream) differs in results from 2.2.2 (bundle). Not checked.
- A GitHub repository `stajichlab/qRNAmap` must be created by the user or with their approval.
