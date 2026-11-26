# Browser Page Version 1 Results and Documentation

Date: 2026-11-23. Branch: `page-v1`. Status: Complete.

## Summary

`qRNAmap` version 1 provides a fully client-side browser application at `/analyze/`.
The tool runs bulk RNA-seq pseudo-alignment and quality control inside the user's browser using WebAssembly (WASM).
No sequencing reads, genome data, or user files are uploaded to any server.

## Architecture

The browser application comprises three main components:

1. **`qrnamap-wasm`** (`crates/qrnamap-wasm`):
   - Rust core compiled to `wasm32-unknown-unknown` with `wasm-bindgen`.
   - Handles genome catalog listing, transcriptome construction from FASTA and GFF3 (`build_transcriptome`), alias table extraction (`gene_info_tsv`), per-sample transcript-to-gene summation (`gene_counts`), count matrix assembly (`assemble_matrix`), sample correlation matrix calculation, and PCA coordinates (`qc`).
2. **`kallisto 0.52.0` WASM**:
   - Vendored C++ binary compiled to WebAssembly via Emscripten.
   - Executed through an isolated Web Worker using `ToolRunner` with a virtual in-memory filesystem (MEMFS).
3. **Dedicated Prep Worker & Memory Isolation**:
   - Transcript extraction and tx2gene parsing run in a transient Web Worker that terminates upon completion.
   - This releases JavaScript memory before the kallisto index builder executes in the main page worker, keeping peak memory safely below browser limits.
4. **Streaming Read Fetcher (`streamCutFastq`)**:
   - For public test sets, FASTQ files stream over HTTPS from the European Nucleotide Archive (ENA).
   - Decompressed chunks pass through `VerifyingGunzip` (fflate with strict CRC-32 and length validation).
   - The stream counts newlines and cancels the HTTP request immediately upon reading $N \times 4$ lines ($N$ pairs).

## Limits and Margins

Measurements in Task 3 established the following operational boundaries:

| Parameter | Value | Rationale |
|---|---|---|
| **Max read pairs ($N$), reads $\le 101$ bp** | 200,000 | Tested on *S. cerevisiae* HiSeq data. Peak memory remains safely under 1 GiB. |
| **Max read pairs ($N$), reads $> 101$ bp** | 75,000 | Determined from 150 bp *C. neoformans* NovaSeq measurements to stay below WebKit's 1.5 GiB worker ceiling. |
| **Max genome size (custom upload)** | 50,000,000 bp | Enforced by client validation to prevent browser out-of-memory errors. |
| **Catalog genomes tested** | Up to 29.4 Mb | Verified on *S. cerevisiae* (12.2 Mb), *C. neoformans* H99 (18.9 Mb), and *A. fumigatus* Af293 (29.4 Mb). |

## Browser Support and Verification

Tests ran across Playwright headless browsers:

| Engine | Genome Step | Samples & WASM Quant | Results View & Downloads | Test Sets Acceptance |
|---|---|---|---|---|
| **Chromium** (156.0) | Pass | Pass | Pass | Pass (byte-identical) |
| **Firefox** (157.0) | Pass | Pass | Pass | Pass (byte-identical) |
| **WebKit** (27.2) | Pass (CI) | Pass (CI) | Pass (CI) | Pass (CI)* |

\* WebKit runs in Ubuntu 24.04 CI. Local host execution requires desktop GTK/WebKit shared libraries.

### Firefox Output Validation Check
During determinism testing in Task 3, Firefox exhibited a rare formatting anomaly where number tokens were duplicated in the text output of `abundance.tsv` (e.g. `9` written as `999`).
To prevent corrupted counts, `abundance-check.js` validates that:
1. Every row parses into 4 numeric values.
2. The target count matches the transcript catalog count.
3. The sum of `est_counts` matches `n_pseudoaligned` within numeric tolerance.
If a check fails, the sample automatically re-runs once.

## Test Sets and CLI Parity

Two public test sets were evaluated against the native `qrnamap` CLI reference tables:

### 1. Yeast Amphotericin B (`yeast-ampb`)
- **Source**: GSE80357 / PRJNA318684 (*Saccharomyces cerevisiae* S288C).
- **Design**: 3 control vs 3 amphotericin B-treated samples (matched HiSeq 2000 101 bp paired-end).
- **Preview Size**: $N = 200,000$ pairs per sample.
- **CLI Reference Hash**: `32fa3cb829b0ea1dfd2054fb36f72ecdd3f1a4b44dbcfde2f2b61017406bd11c` (`counts.tsv`).
- **Browser Result**: 100% byte-identical to CLI reference in Chromium and Firefox.

### 2. *Cryptococcus neoformans* H99 Lkh1 Deletion (`cneo-lkh1`)
- **Source**: GSE233612 / PRJNA976751 (*Cryptococcus neoformans* var. *grubii* H99).
- **Design**: 3 wild-type vs 3 *lkh1* deletion basal samples (NovaSeq 6000 150 bp paired-end).
- **Preview Size**: $N = 75,000$ pairs per sample.
- **CLI Reference Hash**: `4de7a7582b5d00d61131826dea7120213d62a050f64e5d22efa2d4d6d200b1e2` (`counts.tsv`).
- **Browser Result**: 100% byte-identical to CLI reference in Chromium and Firefox.

### Differential Expression Comparison with Paper Table S3
We compared differential expression calls from the preview depth ($N = 75,000$ pairs) against the published Supplementary Table S3 (Kim et al., *Front Cell Infect Microbiol* 2024;14:1367098):

| Metric | Value |
|---|---|
| Total genes tested | 8,338 |
| Paper Table S3 basal DE genes ($p_{\text{adj}} < 0.05, |\log_2 \text{FC}| \ge 2$) | 409 |
| CLI preview DE genes ($N = 75,000$, $p_{\text{adj}} < 0.05, |\log_2 \text{FC}| \ge 2$) | 0 |
| Overlap | 0 |
| Minimum $p_{\text{adj}}$ observed in preview | 0.3202 |

**Conclusion**: At preview depth (75,000 pairs, ~0.4% of full library sequencing depth), statistical power is insufficient to call genome-wide significant differential expression after Benjamini-Hochberg adjustment across 8,338 genes. The preview provides quality control, library sizing, mapping rates, and exploratory PCA/heatmaps; full-depth differential expression requires complete sequencing depth executed via the exported CLI script.

## Known Gaps

1. **Mobile devices**: High-memory requirements of kallisto index building (250–800 MiB) target desktop browsers. Mobile browsers with strict tab memory caps (e.g. iOS Safari < 1 GiB) may fail during indexing.
2. **Single-end reads**: Requires manual entry of fragment length mean and standard deviation, mirroring the CLI requirement.
3. **Full-depth differential expression**: The browser performs exploratory PCA on $\log_2(\text{CPM} + 1)$ of variable genes. Full DESeq2 differential expression is executed via the downloadable `de-command.txt` and `samples.tsv` on native infrastructure.

## Addendum (2026-11-24): Volcano link-outs and `dbxref` in `gene_info.tsv`

The gene features of the parsed GFF3 carry `Dbxref` entries; `gene_info_tsv` now writes them to a new final column `dbxref` (comma-separated, for example `GeneID:854598,SGD:S000007260`), so the page can link genes to their source databases. The header gains the column; all other columns and outputs stay byte-identical.

In the results step, the volcano plot's labeled genes (top five up, top five down) now:

- show a mouseover tooltip with the locus, symbol, log2FC, padj and the linked databases;
- link out on click (label text and the point itself, both `target="_blank"`): `SGD:` xrefs open `https://www.yeastgenome.org/locus/<SGDID>`; `FungiDB:` xrefs and `AFU*_` Aspergillus fumigatus locus tags (Af293 `AFUA_1G00100`, A1163 `AFUB_102960`) open `https://fungidb.org/gene/<ID>`; a numeric `GeneID` opens `https://www.ncbi.nlm.nih.gov/gene/<GeneID>`. SGD takes precedence, then FungiDB; genes with neither stay unlinked. The builders live in `site/analyze/gene-links.js` (pure, unit-tested in `tests/units.test.mjs`).

Note: A1163's GenBank annotation (GCA_000150145.1) carries no `Dbxref` at all, so its genes link only through the `AFUB_` locus-tag rule; Af293's RefSeq GFF (GCF_000002655.1) carries both `FungiDB:` and `GeneID:` xrefs.

The same fix defines the previously missing `downloadBlob` helper, so the `differential_expression.tsv` button works (it threw `ReferenceError: downloadBlob is not defined`).

## Addendum (2026-11-26): the in-page DE test was on the wrong scale

The volcano's test divided a log2 fold change by a standard error computed from normalized counts on the linear
scale (`sqrt(var/n + 0.05)`). The SE grows with read depth, so deeper runs found less: a full-depth
*A. fumigatus* run gave a median SE of about 148 log2 units, a largest |stat| of 3.2 and no gene at padj < 0.05.
On simulated 3-vs-3 data with 20 genes up fourfold, the old test found none of them at any depth. The
"insufficient power" conclusion of the Table S3 comparison above was measured with this test and should be re-run.

The test now lives in `site/analyze/de-stats.js` (pure, unit-tested in `tests/units.test.mjs`):

- DESeq2 median-of-ratios size factors, then log2(normalized + 0.5) per sample;
- a two-sample t-test on that scale with a limma-style moderated variance
  `s2_post = (d0*s0^2 + d*s2) / (d0 + d)`, `s0^2` the median gene variance among genes with mean ≥ 10, `d0 = 4`;
- p-values from Student's t on `d + d0` degrees of freedom, then Benjamini-Hochberg;
- no replicate in either group: p = 1 (nothing to test against).

On the simulation it finds all 20 changed genes at 1x, 100x and 10,000x depth, with at most 5 of 380 false calls. It is
a screen, not DESeq2: `qrnamap de` remains the analysis of record. The padj cutoff is now inclusive (`padj ≤ cutoff`),
so a cutoff of 1 labels the top genes even when none is significant. The correlation heatmap's column labels are
rotated -45° and the margins grow with the longest sample name, so long names no longer overlap.

## Addendum (2026-11-26): gene counts table, gene function and popup

`gene_info.tsv` gains a final column `description`: the gene's function text from the GFF3. Sources, in order: the gene
feature's `description`, `Note` or `product`; else the `product` of its first transcript that has one; else that
transcript's CDS `product` (NCBI RefSeq puts the product on the mRNA and CDS). Tabs and line breaks become spaces. The
page reads columns by header name, so older tables (without `dbxref` or `description`) still load. Pre-built index sets
made before this change have no `description` column: their genes show no function text until the set is rebuilt.

The gene counts table:

- shows 25, 50, 100 or 250 rows per page, with first / previous / next / last buttons and a page number box;
- sits in a scroll box with a sticky header that can be dragged taller, or switched between short and tall (80% of the window);
- has an optional "Gene info" column (symbol and function, wrapped to two lines; full text on hover);
- searches the function text as well as the IDs;
- opens a popup on a gene name with the function, IDs, location and the same link-outs as the volcano
  (`genePopupHtml` in `site/analyze/gene-links.js`, unit-tested).

Rows per page, the info column and the table height are remembered in the browser (`localStorage`, when available).
