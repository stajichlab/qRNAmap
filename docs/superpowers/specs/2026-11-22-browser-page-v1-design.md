# qRNAmap browser page, version 1: design

Date: 2026-11-22. Status: draft for approval. Parent: `2026-11-19-qrnamap-design.md` (section 9, the page). Evidence: `docs/browser/spike0.md` (spike 0 results, 2026-11-21).

## Intent

Today `https://qrnamap.fungalgenomes.org/` is a placeholder. It links to no tool, because no tool exists. This design is for the first working page: a lab member or a collaborator brings their own RNA-seq reads, picks a genome, and gets gene counts and two quality plots. All computation runs in the visitor's browser. No server and no upload of reads.

Why: the lab and its collaborators need a first look at their data without installing tools or using a cluster. The command-line tool stays the route for full-size runs and for differential expression.

Users: the Stajich lab and collaborating labs. They bring their own reads. We also offer public test sets.

## Decisions (from the design conversation, 2026-11-22)

1. **Quantification mode only (kallisto 0.52.0).** HISAT2 plus featureCounts (alignment mode) is the next release.
2. **No differential expression in version 1.** The page ends at counts and quality plots. It gives the exact `qrnamap de` command for the next step.
3. **Genomes:** a menu of the built-in catalog (155 fungal genomes, NCBI RefSeq except three GenBank assemblies), downloaded from NCBI by the page. NCBI sets are preferred. Also: upload of the user's own FASTA and GFF3. Both paths feed the same transcript builder.
4. **Reads:** local FASTQ files (plain or gzip), plus a "Load test set" button. No free ENA accession box in version 1.
   **Sample info (added 2026-11-22, user request):** two ways in. (a) Enter it in a table in the page: sample name, condition and other columns, and a mate 1 and mate 2 file picker per row. (b) Upload a CLI-format sample sheet and pick the FASTQ files; the page matches each row's `fastq1`/`fastq2` to a picked file by file name, ignoring the folder part. Unmatched rows are flagged with a picker on the row; a duplicate file name is a clear error. The sheet the page writes has only file names, because the page cannot know disk paths; an optional "reads folder" box adds the user's cluster path to the downloaded sheet and to the `qrnamap de` command, and the page says so when it is empty.
5. **Large inputs:** (amended 2026-11-22, user decision: N depends on the read length; the page takes N from `previewLimit(readLength)`; 200,000 pairs for reads of up to 101 bases, the value for longer reads from the page quantification memory job, see docs/browser/page-v1-measurements.md) a preview of the first N read pairs per sample, labelled "preview of the first N pairs" on every table and plot, with the full command-line command shown. Above a hard limit the page stops with a message. N and the limit are measured, not guessed (see "To measure").
6. **Results step:**
   - mapping summary per sample (reads processed, reads pseudoaligned, fraction);
   - sortable, searchable gene count table (transcript counts summed per gene with the same transcript-to-gene table the CLI uses, rounded as the CLI does);
   - sample-to-sample correlation heatmap and PCA;
   - downloads: gene count matrix, sample sheet, the filled-in `qrnamap de` command;
   - results saved in the browser (cleared from the page, never sent anywhere).
   - Not in version 1: volcano and MA plots (they need DE).
7. **Test sets:**
   - Yeast amphotericin B (GSE80357): AmpB against control, 3 against 3. Native CLI results exist.
   - *C. neoformans* H99, Lkh1 study (BioProject PRJNA976751, GEO GSE233612): *lkh1* deletion against wild type, basal condition, 3 against 3. The paper reports 276 genes up and 133 down (P < 0.05, |log2 fold change| of at least 2; Supplementary Table S3).
   - *N. crassa*: later, when a collaborator names a set. The candidate studied (PRJDB15450) failed the rules: replicate counts could not be confirmed and the paper is not open access.
8. **Layout:** `/` is the landing page with an "Analyze reads" button, the test sets, the CLI install link and the spike 0 results. `/analyze/` is the tool: Genome, Samples, Results on one page. Hosting stays on GitHub Pages. kallisto needs no isolation headers.
9. **kallisto files:** commit the built files to the repo (about 1.7 MiB), with a provenance file.
11. **Gene aliases (added 2026-11-22, user request):** gene names stay the GFF `locus_tag` (fallback `Name`, then `ID`), as now; the count table's first column stays the locus name. Transcript names in the index stay the GFF transcript `ID` (NCBI files: the mRNA accession with an `rna-` prefix). A new file `gene_info.tsv`, one row per gene, carries the aliases: locus name, gene symbol (`Name`), `GeneID` (from `Dbxref`), `old_locus_tag`, transcript accessions, protein accessions (from the CDS `protein_id`), biotype and class. The page's count table search matches any alias, and the page offers `gene_info.tsv` as a download. The CLI writes the same file (page and CLI stay equal). Existing outputs, including `counts.tsv`, do not change.
10. **Page code:** copy `gz.js`, `reads.js`, `reference-panel.js`, `persist.js` and `theme.js` from `../qmutfind/site` and record the `qmutfind` commit in each copied file header. No shared package.

## Goals

- A visitor can go from the landing page to a count table for a test set, in a supported browser, without installing anything.
- The page's gene count table equals the command-line tool's table for the same reads and the same preview size (acceptance rule below).
- Every limit the page enforces comes from a measurement.

## Non-goals

- Alignment mode, featureCounts, HISAT2, DESeq2, Salmon in the browser.
- Differential expression, volcano and MA plots.
- A free ENA accession box. FungiDB downloads (needs a key, no CORS: CLI only).
- A server, accounts, or any upload of reads.
- A shared package with `qmutfind`.

## Anti-goals (a failure even if it technically works)

- A browser result that differs from the CLI result without a test that shows and explains it.
- A preview presented as a full analysis.
- A claim that "it works in browsers" beyond what was measured (headless Chromium, Firefox, WebKit on a 16 GB runner, one small dataset).
- An API key or other secret in page code.

## Approaches considered

**Where the kallisto files come from.**
1. **Commit the built files (chosen by the user).** About 1.7 MiB. Fast deploy, no Emscripten in the Pages workflow, the live page depends on nothing outside qRNAmap. Cost: a binary in git, which is checked by recording the SHA-256 and the upstream build run.
2. Build kallisto in the Pages workflow. Reproducible, no binary in git. Cost: a longer deploy (the upstream build job took 248 s) and a new way to break the site. Can replace option 1 later without changing the page.
3. Fetch from a release in `stajichlab/genomics-wasm` at page load. Small repo, but the live page depends on a second repository.

**Reads.** Local files (chosen) against a free ENA accession box (not in version 1) against both. ENA fetch is needed only for the test sets, so those can be verified in advance.

**Page code.** Copy the `qmutfind` modules with the commit recorded (chosen) against a shared package (rejected: ties the two release schedules together for little gain).

## Facts used (measured, from spike 0)

- `kallisto.wasm` is 1,669,837 bytes (512,352 gzipped); `kallisto.mjs` is 92,178 bytes. Built from kallisto 0.52.0 revision `4e9f29cf3b021260415430c057a22469ca081391`, zlib 1.3.1, Emscripten 3.1.64 (upstream `sources.lock.json` at the pin 0b015f6; our build run 44939577353).
- kallisto `abundance.tsv` from the WASM build is byte-identical to native in Node, Chromium, Firefox and WebKit on a 200,000-pair yeast subset.
- kallisto index build in the browser: yeast transcripts 19.34 to 37.03 s; *A. fumigatus* Af293 transcripts 43.96 to 84.41 s. Peak memory up to 1,176.5 MiB (WebKit, whole-browser process tree). The WASM heap limit is 2 GiB (a stated limit, not a measurement).
- kallisto runs without cross-origin isolation. The WASM build requires `-t 1`.
- Browsers can read NCBI genome FASTA and GFF3 files, and ENA FASTQ files including byte ranges (HTTP 206), from a page. Measured 2026-11-21 on one network.
- kallisto index bytes differ from run to run, and the cause is unknown. The quantification output was identical. The page must not compare or cache index bytes.
- A yeast sample, SRR3396381, has 8,405,534 read pairs; its mate 1 file is 669,296,728 bytes.

## To measure before fixing numbers

1. Preview size N and the hard limit: kallisto quantification (not only indexing) in the three browsers at several read counts, with time and peak memory, 3 runs each.
2. Whether a truncated gzip stream from ENA goes through kallisto in the browser. UNVERIFIED. Fallback if it fails: host small subset FASTQ files (public data) as release assets and have the test-set button fetch those.
3. Whether the *C. neoformans* paper's gene identifiers match the RefSeq H99 annotation in the catalog (GCF_000149245.1). UNVERIFIED.
4. The yeast set's other five samples: read counts and ENA file sizes. UNVERIFIED (only SRR3396381 was checked).
5. Behaviour on a phone-sized screen. Not measured.
6. Genome menu limits: download size and transcript-building time for the larger genomes in the catalog.

Results (2026-11-22): `docs/browser/page-v1-measurements.md`. In short (headless Playwright browsers on a GitHub runner, SRR3396381, yeast transcripts):
item 1: N = 200,000 pairs, hard limit = 200,000 pairs (for 101-base reads; amended 2026-11-22: N depends on the read length, see below), reads given to kallisto as plain FASTQ; item 2: a cut gzip goes through kallisto (exit 0, it stops at the
shorter mate), and cutting both mates at N records works in all three browsers, but ENA failed in 11 of 30 sessions and the release-asset fallback is blocked by CORS
(a same-origin fallback is UNVERIFIED; the CORS test used a release asset of stajichlab/genomics-wasm); item 3: the 409 genes of Supplementary Table S3
match the RefSeq H99 locus tags and the qRNAmap gene names (column-to-run mapping UNVERIFIED); item 4: measured (ENA portal API); item 5: not measured;
item 6: the three largest catalog genomes finished in all engines.
Firefox wrote a damaged `abundance.tsv` in 4 of 130 kallisto runs (about 3%, no error, cause UNVERIFIED). So a byte-identical comparison of Firefox output
with the CLI (acceptance rule) can fail without a code fault: Task 9 must re-run or check the file before it reports a mismatch. Needed page behaviour:
check each `abundance.tsv` and re-run the sample once when the check fails (untested).

## Acceptance rule

For each test set, with the same reads and the same preview size, the page's gene count table equals the CLI's table. Gene summing, rounding and matrix building use `qrnamap-core` compiled to WASM in a `qrnamap-wasm` wrapper, so the CLI and the page share one code path. A difference blocks the release until it is explained.

## Constraints

- `qseq-io`, `genome-catalog` and `qrnamap-core` stay pure (no I/O) and must keep building for `wasm32-unknown-unknown`.
- Licences: qRNAmap is GPL-3.0-only; the copied `qmutfind` modules are GPL-3.0-only (same author). kallisto is BSD-2-Clause, with bundled zlib, Bifrost (BSD-2-Clause), CRoaring (Apache-2.0), BooPHF (MIT), libdivide, kseq and wyhash. The page ships the upstream notices next to the kallisto files. See `docs/browser/spike0.md`, licence section. No legal conclusion is drawn here.
- No secrets in any committed file. Reads never leave the browser.
- Docs in Simplified Technical English. Numbers on the page come from measurements.

## Left to the builder

The worker design and message protocol; progress and cancel behaviour; the layout of the three steps; how the Pages workflow builds `qrnamap-wasm` (the `qmutfind` workflow is the template: `pixi run build-wasm-web`); the form of the provenance file; error text for bad uploads (FASTA names that do not match the GFF3, annotations with no usable transcripts); test layout, including a headless-browser end-to-end test on the test sets.

## Risks

- A truncated gzip may not work through kallisto. The fallback is in "To measure".
- Memory on small devices. All measurements used a 16 GB runner.
- A committed binary can drift from its recorded source. The provenance file and a CI check of the SHA-256 reduce this.
- Test-set data on ENA can change or move. The page must show a clear error if a fetch fails.
