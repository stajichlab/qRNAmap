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
