#!/usr/bin/env bash
# run_amb_salmon_overhang.sh: the AMB data with `qrnamap run --quantifier salmon` and the two Salmon
# overhang settings, on the UCR HPCC. Same stages and comparisons as run_amb_salmon.sh.
#
# Usage: REPO=/ABS/PATH SETTING=A|B bash /ABS/PATH/validation/amb/run_amb_salmon_overhang.sh STAGE
# REPO is required: the absolute path of the qRNAmap checkout. SETTING is required:
#   A  --salmon-softclip-overhangs
#   B  --salmon-softclip-overhangs --salmon-min-score-fraction 0.5
# All paths are absolute and come from the variables below (override them in the environment).
# This script never uses BASH_SOURCE (the compare stage runs inside a SLURM job) and never
# hardcodes /scratch: the jobs use "${SCRATCH:?}".
#
# Salmon 2.8.0 does not record --softclipOverhangs or --minScoreFraction in any output file.
# The record of the setting is: the `gen` log (it prints the qrnamap argv) and the generated
# $SAL/slurm/samples.sbatch (it holds the qrnamap command with the options).
#
# Stages:
#   gen        qrnamap run --quantifier salmon --slurm + setting options (default --mem-gb 8,
#              no --tmpdir)
#   submit     run $SAL/slurm/submit.sh from a fresh login shell, so that the jobs get $SCRATCH
#   compare    SLURM job body: compare_modes.R (setting vs kallisto; setting vs default Salmon),
#              qrnamap de on the setting counts (check), compare_geo.R (authors vs setting vs
#              kallisto)
#   collect    login node: collect_salmon.py and zero_genes_by_length.py (needs TASK_START, the
#              sacct Start of the sample array task)
set -euo pipefail

REPO="${REPO:?set REPO to the absolute path of the qRNAmap checkout (scripts, pixi env, results)}"
SETTING="${SETTING:?set SETTING to A or B}"
DATA="${DATA:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb}"
PIXI_BIN="${PIXI_BIN:-$REPO/.pixi/envs/default/bin}"
QRNAMAP="${QRNAMAP:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-overhang/release/qrnamap}"
GENOME_DIR="${GENOME_DIR:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/genome}"
REF="$GENOME_DIR/GCF_000146045.2_R64_genomic.fna.gz"
GFF="$GENOME_DIR/GCF_000146045.2_R64_genomic.gff.gz"
SHEET="$DATA/samples.tsv"
CONTRAST="condition,AmpB,Control"
KAL="$DATA/run_quant"
SAL0="$DATA/run_quant_salmon"

case "$SETTING" in
A) OPTS=(--salmon-softclip-overhangs) ;;
B) OPTS=(--salmon-softclip-overhangs --salmon-min-score-fraction 0.5) ;;
*) echo "SETTING must be A or B" >&2; exit 2 ;;
esac
SAL="$DATA/run_quant_salmon_overhang_$SETTING"
RES="$REPO/validation/amb/results_salmon_overhang/$SETTING"

stage="${1:?stage}"; shift || true

case "$stage" in
gen)
    export PATH="$PIXI_BIN:$PATH"
    cmd=("$QRNAMAP" run --samples "$SHEET" --reference "$REF" --gff "$GFF"
        --mode quant --quantifier salmon "${OPTS[@]}" --out "$SAL" --strandedness unstranded
        --threads 1 --contrast "$CONTRAST" --slurm --samples-per-task 6
        --partition short --time 02:00:00)
    echo "setting $SETTING argv: ${cmd[*]}"
    "${cmd[@]}"
    echo "options in samples.sbatch:"
    grep -n -e '--salmon-' "$SAL/slurm/samples.sbatch"
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
    TMP="${SCRATCH:?SCRATCH is not set: submit from a login shell}/salmon_overhang_compare_$SETTING"
    mkdir -p "$TMP" "$RES/geo" "$RES/vs_salmon_default"
    echo "host $(hostname) job ${SLURM_JOB_ID:-none} setting $SETTING start $(date -Is)"
    command -v Rscript zstd
    # (a) setting vs kallisto. In the output tables "quant" is this setting and "align" is kallisto.
    Rscript "$REPO/validation/amb/compare_modes.R" "$SAL" "$KAL" "$RES"
    # (a2) setting vs default Salmon. "quant" is this setting, "align" is default Salmon.
    Rscript "$REPO/validation/amb/compare_modes.R" "$SAL" "$SAL0" "$RES/vs_salmon_default"
    # Check of the DESeq2 step: qrnamap de again on the setting counts with this binary.
    RERUN="$DATA/geo/de_salmon_overhang_${SETTING}_rerun"
    mkdir -p "$RERUN"
    "$QRNAMAP" de --counts "$SAL/counts.tsv.zst" --samples "$SHEET" \
        --design '~ condition' --contrast "$CONTRAST" \
        --out "$RERUN/deseq2" --tmpdir "$TMP"
    # (b) GEO comparison: setting in the QUANT slot, kallisto in the ALIGN slot.
    Rscript "$REPO/validation/geo/compare_geo.R" \
        "$DATA/geo/prep/authors_counts_amb.tsv" \
        "$SAL/counts.tsv.zst" "$KAL/counts.tsv.zst" \
        "$DATA/geo/download/GSE80357_log_fold_change_and_FDR_values.csv.gz" \
        "$SAL/deseq2.results.csv" "$KAL/deseq2.results.csv" \
        "$DATA/geo/de_authors/deseq2.results.csv" "$RERUN/deseq2.results.csv" \
        "$DATA/geo/download/GSE80357_average_cpm_per_treatment_group.csv.gz" \
        "$RES/geo"
    rm -rf "$TMP"
    echo "end $(date -Is)"
    ;;
collect)
    export PATH="$PIXI_BIN:$PATH"
    TASK_START="${TASK_START:?set TASK_START to the sacct Start of the sample array task}"
    mkdir -p "$RES"
    python -I "$REPO/validation/amb/collect_salmon.py" --run-dir "$SAL" --sheet "$SHEET" \
        --task-start "$TASK_START" --out "$RES/per_sample.tsv"
    python -I "$REPO/validation/amb/zero_genes_by_length.py" --sheet "$SHEET" --salmon-dir "$SAL" \
        --kallisto-dir "$KAL" --authors "$DATA/geo/prep/authors_counts_amb.tsv" --out-dir "$RES"
    zstd -dc "$SAL/mapping_summary.tsv.zst" > "$RES/mapping_summary.tsv"
    ;;
*)
    echo "unknown stage: $stage" >&2; exit 2 ;;
esac
