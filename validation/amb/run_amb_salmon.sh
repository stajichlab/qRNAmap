#!/usr/bin/env bash
# run_amb_salmon.sh: the AMB data with `qrnamap run --quantifier salmon` on the UCR HPCC, and the
# comparison with the earlier kallisto run ($DATA/run_quant) and the authors' GSE80357 tables.
# Adapted from run_amb.sh (stages slurm and submit-login) and validation/geo/geo_job.sh.
#
# Usage: REPO=/ABS/PATH bash /ABS/PATH/validation/amb/run_amb_salmon.sh STAGE [ARGS]
# REPO is required: the absolute path of the qRNAmap checkout.
# All paths are absolute and come from the variables below (override them in the environment).
# This script never uses BASH_SOURCE (the compare stage runs inside a SLURM job) and never
# hardcodes /scratch: the jobs use "${SCRATCH:?}".
#
# Stages:
#   gen        qrnamap run --quantifier salmon --slurm (default --mem-gb 8, no --tmpdir)
#   submit     run $DATA/run_quant_salmon/slurm/submit.sh from a fresh login shell (as
#              run_amb.sh submit-login), so that the jobs get $SCRATCH
#   compare    SLURM job body: compare_modes.R (Salmon vs kallisto), qrnamap de on the Salmon
#              counts (check), compare_geo.R (authors vs Salmon vs kallisto)
set -euo pipefail

REPO="${REPO:?set REPO to the absolute path of the qRNAmap checkout (scripts, pixi env, results)}"
DATA="${DATA:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb}"
PIXI_BIN="${PIXI_BIN:-$REPO/.pixi/envs/default/bin}"
QRNAMAP="${QRNAMAP:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-salmon/release/qrnamap}"
GENOME_DIR="${GENOME_DIR:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/genome}"
REF="$GENOME_DIR/GCF_000146045.2_R64_genomic.fna.gz"
GFF="$GENOME_DIR/GCF_000146045.2_R64_genomic.gff.gz"
SHEET="$DATA/samples.tsv"
CONTRAST="condition,AmpB,Control"
SAL="$DATA/run_quant_salmon"
KAL="$DATA/run_quant"
RES="$REPO/validation/amb/results_salmon"

stage="${1:?stage}"; shift || true

case "$stage" in
gen)
    export PATH="$PIXI_BIN:$PATH"
    "$QRNAMAP" run --samples "$SHEET" --reference "$REF" --gff "$GFF" \
        --mode quant --quantifier salmon --out "$SAL" --strandedness unstranded --threads 1 \
        --contrast "$CONTRAST" --slurm --samples-per-task 6 \
        --partition short --time 02:00:00
    ;;
submit)
    env -i HOME="$HOME" USER="$USER" PATH="/usr/bin:/bin" bash -l -c '
        echo "submit shell PATH=$PATH"
        for t in salmon kallisto Rscript; do
            printf "  %s: %s\n" "$t" "$(command -v "$t" || echo NOT-FOUND)"
        done
        bash "$0"' "$SAL/slurm/submit.sh"
    ;;
compare)
    export PATH="$PIXI_BIN:$PATH"
    TMP="${SCRATCH:?SCRATCH is not set: submit from a login shell}/salmon_compare"
    mkdir -p "$TMP" "$RES/geo"
    echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -Is)"
    command -v Rscript zstd
    # (a) Salmon vs kallisto. In the output tables "quant" is Salmon and "align" is kallisto.
    Rscript "$REPO/validation/amb/compare_modes.R" "$SAL" "$KAL" "$RES"
    # Check of the DESeq2 step: qrnamap de again on the Salmon counts with this binary.
    mkdir -p "$DATA/geo/de_salmon_rerun"
    "$QRNAMAP" de --counts "$SAL/counts.tsv.zst" --samples "$SHEET" \
        --design '~ condition' --contrast "$CONTRAST" \
        --out "$DATA/geo/de_salmon_rerun/deseq2" --tmpdir "$TMP"
    # (b) The GEO comparison script with Salmon in the QUANT slot and kallisto in the ALIGN
    # slot. In the output tables "quant" is Salmon and "align" is kallisto.
    Rscript "$REPO/validation/geo/compare_geo.R" \
        "$DATA/geo/prep/authors_counts_amb.tsv" \
        "$SAL/counts.tsv.zst" "$KAL/counts.tsv.zst" \
        "$DATA/geo/download/GSE80357_log_fold_change_and_FDR_values.csv.gz" \
        "$SAL/deseq2.results.csv" "$KAL/deseq2.results.csv" \
        "$DATA/geo/de_authors/deseq2.results.csv" "$DATA/geo/de_salmon_rerun/deseq2.results.csv" \
        "$DATA/geo/download/GSE80357_average_cpm_per_treatment_group.csv.gz" \
        "$RES/geo"
    rm -rf "$TMP"
    echo "end $(date -Is)"
    ;;
*)
    echo "unknown stage: $stage" >&2; exit 2 ;;
esac
