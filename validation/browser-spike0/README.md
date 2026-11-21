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
