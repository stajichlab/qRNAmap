#!/usr/bin/env bash
# salmon_subset.sh: SLURM job body. Salmon option variants on the first N pairs of one AMB sample,
# against kallisto, for the short transcripts RPL41A (YDL184C), RPL41B (YDL133C-A) and YOR302W
# (78 nt each, reads 101 nt). Uses the indexes of the earlier runs (read only). Diagnostic only:
# qrnamap passes none of these tuning flags.
#
# Usage: sbatch -p short -c 2 --mem=8G --time=00:30:00 --output=LOG \
#            --wrap "REPO=/ABS/worktree bash /ABS/worktree/validation/amb/salmon_subset.sh SAMPLE PAIRS"
# REPO must be set (no BASH_SOURCE). Temporary files go under "${SCRATCH:?}". Writes
# $REPO/validation/amb/results_salmon/subset_SAMPLE.tsv.
set -euo pipefail
REPO="${REPO:?set REPO to the absolute path of the qRNAmap checkout}"
DATA="${DATA:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb}"
SAMPLE="${1:?sample}"; PAIRS="${2:?pairs}"
export PATH="$REPO/.pixi/envs/default/bin:$PATH"
T="${SCRATCH:?SCRATCH is not set: submit from a login shell}/salmon_subset"
mkdir -p "$T"
echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -Is)"
salmon --version; kallisto version
f1=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2}' "$DATA/samples.tsv")
f2=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $3}' "$DATA/samples.tsv")
lines=$((PAIRS * 4))
# zcat | head ends zcat with SIGPIPE; pipefail would stop the script, so allow it here.
(zcat "$f1" | head -n "$lines" > "$T/r1.fq") || true
(zcat "$f2" | head -n "$lines" > "$T/r2.fq") || true
echo "pairs in subset: $(( $(wc -l < "$T/r1.fq") / 4 )) $(( $(wc -l < "$T/r2.fq") / 4 ))"
TX="rna-NM_001180244.1 rna-NM_001180194.1 rna-NM_001183721.1"   # RPL41A RPL41B YOR302W

out="$REPO/validation/amb/results_salmon/subset_${SAMPLE}.tsv"
printf 'variant\toptions\tpairs\tpercent_mapped\tRPL41A\tRPL41B\tYOR302W\n' > "$out"
run_salmon() {
    local name="$1"; shift
    salmon quant -i "$DATA/run_quant_salmon/salmon_index" -l IU -1 "$T/r1.fq" -2 "$T/r2.fq" \
        -p 2 -o "$T/$name" "$@" 2> "$T/$name.log"
    local pm np
    pm=$(grep -o '"percent_mapped": [0-9.]*' "$T/$name/aux_info/meta_info.json" | grep -o '[0-9.]*$')
    np=$(grep -o '"num_processed": [0-9]*' "$T/$name/aux_info/meta_info.json" | grep -o '[0-9]*$')
    local vals=()
    for t in $TX; do vals+=("$(awk -F'\t' -v t="$t" '$1==t{print $5}' "$T/$name/quant.sf")"); done
    printf '%s\tsalmon %s\t%s\t%s\t%s\t%s\t%s\n' "$name" "${*:-(defaults)}" "$np" "$pm" "${vals[@]}" >> "$out"
}
run_salmon defaults
run_salmon softclipOverhangs --softclipOverhangs
run_salmon softclip --softclip
run_salmon softclipOverhangs_msf0.5 --softclipOverhangs --minScoreFraction 0.5
run_salmon useEM --useEM
kallisto quant -i "$DATA/run_quant/transcripts.kidx" -o "$T/kallisto" --plaintext -t 2 \
    "$T/r1.fq" "$T/r2.fq" 2> "$T/kallisto.log"
vals=()
for t in $TX; do vals+=("$(awk -F'\t' -v t="$t" '$1==t{print $4}' "$T/kallisto/abundance.tsv")"); done
np=$(grep -o '"n_processed": [0-9]*' "$T/kallisto/run_info.json" | grep -o '[0-9]*$')
pm=$(grep -o '"p_pseudoaligned": [0-9.]*' "$T/kallisto/run_info.json" | grep -o '[0-9.]*$')
printf 'kallisto\tkallisto 0.52.0 (unstranded)\t%s\t%s\t%s\t%s\t%s\n' "$np" "$pm" "${vals[@]}" >> "$out"
cat "$out"
rm -rf "$T"
echo "end $(date -Is)"
