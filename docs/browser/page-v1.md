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
