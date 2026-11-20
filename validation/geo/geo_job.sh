#!/usr/bin/env bash
# geo_job.sh: SLURM job body of the GSE80357 comparison (DESeq2 runs and compare_geo.R).
#
# Usage: sbatch -p short -c 1 --mem=4G --time=00:30:00 --output=LOG \
#            --wrap "bash /ABS/validation/geo/geo_job.sh REPO DATA QRNAMAP"
#   REPO     absolute path of the worktree (scripts, pixi env, results/)
#   DATA     absolute path of the AMB data directory (run_quant, run_align, samples.tsv, geo/)
#   QRNAMAP  absolute path of the qrnamap binary
# No BASH_SOURCE; no hardcoded /scratch: temporary files go under "${SCRATCH:?}".
set -euo pipefail
REPO="${1:?repo}"; DATA="${2:?data}"; QRNAMAP="${3:?qrnamap binary}"
export PATH="$REPO/.pixi/envs/default/bin:$PATH"
TMP="${SCRATCH:?SCRATCH is not set: submit from a login shell}/geo_compare"
mkdir -p "$TMP"
GEO="$DATA/geo"
echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -Is)"
command -v Rscript zstd

# 4. DESeq2 on the authors' counts (all authors' genes, our sample names, the AMB sheet).
mkdir -p "$GEO/de_authors"
"$QRNAMAP" de --counts "$GEO/prep/authors_counts_amb.tsv" --samples "$DATA/samples.tsv" \
    --design '~ condition' --contrast condition,AmpB,Control \
    --out "$GEO/de_authors/deseq2" --tmpdir "$TMP"
# Check: the same step on our quant counts with this binary (vs run_quant/deseq2.results.csv).
mkdir -p "$GEO/de_quant_rerun"
"$QRNAMAP" de --counts "$DATA/run_quant/counts.tsv.zst" --samples "$DATA/samples.tsv" \
    --design '~ condition' --contrast condition,AmpB,Control \
    --out "$GEO/de_quant_rerun/deseq2" --tmpdir "$TMP"

Rscript "$REPO/validation/geo/compare_geo.R" \
    "$GEO/prep/authors_counts_amb.tsv" \
    "$DATA/run_quant/counts.tsv.zst" "$DATA/run_align/counts.tsv.zst" \
    "$GEO/download/GSE80357_log_fold_change_and_FDR_values.csv.gz" \
    "$DATA/run_quant/deseq2.results.csv" "$DATA/run_align/deseq2.results.csv" \
    "$GEO/de_authors/deseq2.results.csv" "$GEO/de_quant_rerun/deseq2.results.csv" \
    "$GEO/download/GSE80357_average_cpm_per_treatment_group.csv.gz" \
    "$REPO/validation/geo/results"
rm -rf "$TMP"
echo "end $(date -Is)"
