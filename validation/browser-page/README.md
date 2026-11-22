# Browser page checks

The workflow `.github/workflows/page.yml` does not download or convert reads from NCBI on each run.
It downloads static files from the release `page-v1-test-data-v1` (pre-release) of this repository.

- Default (push, or dispatch with all `counts` <= 1,000,000): the job `reads` downloads the two prefix files, the first 1,000,000 read pairs of each mate of SRR3396381 (83,401,396 and 83,784,618 bytes).
- Full path (the largest value in `counts` > 1,000,000): the job `reads` downloads the two complete mates (712,250,441 and 723,342,896 bytes; 8,405,534 pairs each).
- The SHA-256 of every file is checked against `measure/reads_prefix.tsv`, which is committed in this repository (not against the release). A replaced asset fails the job.
- The prefix files were made from the full files with `zcat SRR3396381_N.fastq.gz | head -n 4000000 | gzip -6 -n`. The first 4,000,000 decompressed lines are identical to those of the full files (SHA-256 compared).
- The full files were made by `measure/prepare_reads.sh` (NCBI SRA, fasterq-dump 3.2.0, gzip -6 -n) in run 44990656616. CI does not run this script. Use it to re-create the data.
- The script `measure/fetch_reads_release.sh <out dir> <prefix|full>` does the download. It also fetches the yeast transcripts from the spike 0 fixture release (SHA-256 from `validation/browser-spike0/manifest.tsv`).

## Genome step tests (`tests/`, `fixtures/`)

- `tests/genome.test.mjs`: browser tests of step 1 of `/analyze/` (Playwright 1.64.0, Node test runner). Workflow `.github/workflows/page-smoke.yml`, job `genome-step`, one job per engine.
- Run: `bash scripts/build-wasm-pkg.sh web pkg-web`, `bash scripts/build-site.sh`, `npm ci` in `tests/`, then
  `DIST=$PWD/dist-site ENGINE=chromium node --test validation/browser-page/tests/genome.test.mjs`.
- The tests do not contact NCBI. `page.route` answers the two catalog URLs of S. cerevisiae S288C (GCF_000146045.2) with the fixtures.
  The local test server can hold back or replace `vendor/kallisto/kallisto.mjs` (cancel during the index, and a simulated out-of-memory failure).
- The size-limit tests use `?genome_limit_bases=N`. This option can only lower the page limit.
- `REAL_NCBI=1` adds one test that downloads S288C from NCBI through the page. It is not run in CI.
- `fixtures/make_fixtures.py` makes the fixtures from the S288C RefSeq files: the mitochondrial genome NC_001224.1 only
  (`mito.fna`, `mito.gff3`, and gzip copies). `fixtures/PROVENANCE.tsv` has the SHA-256 of the sources and the outputs.
  `fixtures/expected.tsv` has the counts of `qrnamap transcripts` on the fixture (46 genes, 47 transcripts).

## Samples step tests (`tests/samples.test.mjs`, `tests/units.test.mjs`)

- `tests/samples.test.mjs`: browser tests of step 2 of `/analyze/` (sample table, sheet upload, checks, quantification, cancel, the abundance check).
  Each test uploads the fixture genome in step 1 first. Same job as the genome step tests (`genome-step`, one job per engine).
  Run: `DIST=$PWD/dist-site ENGINE=chromium node --test validation/browser-page/tests/samples.test.mjs`.
- The local test server can replace `vendor/kallisto/kallisto.mjs` after the index is built: a module that writes the largest `est_counts` token of
  `abundance.tsv` three times (as Firefox did in the measurements), for the next 1 or 2 runs, or it holds the request (cancel test). The page has no test switch for this.
- `tests/units.test.mjs`: Node tests (no browser) of `site/analyze/abundance-check.js`, `sheet.js` and `fastq.js`. Run: `node --test validation/browser-page/tests/units.test.mjs`.
- `fixtures/make_reads.py <fixtures dir> <qrnamap binary>` makes the read fixtures: reads simulated from the 47 transcripts of the fixture genome
  (fixed seed; `s1` 300 pairs gzip, `s2` 300 pairs plain, `s3` 300 single-end reads gzip), and the expected values from `qrnamap quant` (kallisto 0.52.0 on PATH):
  `fixtures/expected_quant.tsv` (gene counts per case) and `fixtures/expected_runinfo.tsv` (n_processed, n_pseudoaligned). Case `s1_first100` is the first 100 pairs of `s1`.
  The tests check that the page gives the same gene counts as the command-line tool for the same reads. This is a small check of the acceptance rule; Task 9 does the check on the real test sets.
