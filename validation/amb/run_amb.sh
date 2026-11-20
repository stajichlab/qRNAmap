#!/usr/bin/env bash
# run_amb.sh: driver of the AMB experiments (E1 to E3) on the UCR HPCC. Adapted from
# validation/run_yeast.sh.
#
# Usage: bash /ABS/PATH/validation/amb/run_amb.sh STAGE [ARGS]
# All paths are absolute and come from the variables below (override them in the environment).
# This script never uses BASH_SOURCE (some stages run inside SLURM jobs) and never hardcodes
# /scratch: jobs use "${SCRATCH:?}".
#
# Stages:
#   e1-gen                 qrnamap run --mode quant --slurm on a sheet whose FASTQ files do not exist
#   submit-clean RUN_DIR   run RUN_DIR/slurm/submit.sh from `env -i` with only sbatch on PATH (E1a;
#                          on this cluster the jobs then have no $SCRATCH, see README)
#   submit-login RUN_DIR   run RUN_DIR/slurm/submit.sh from a fresh login shell without the pixi env
#   strand                 SLURM job body: qrnamap strand-check on the 6 samples (default --pairs)
#   slurm MODE SPT STRAND MEM_GB TIME
#                          qrnamap run --slurm for MODE (pixi env on PATH of this shell only)
#   compare                SLURM job body: compare_modes.R on run_quant and run_align
#   e3-manual SAMPLE       SLURM job body: the align steps of one sample by hand, with /usr/bin/time -v
#   e3-align SAMPLE KIND   SLURM job body: qrnamap align --task-id 0 for one sample (index copied
#                          from run_align); KIND scratch|gpfs sets --tmpdir
set -euo pipefail

REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/qRNAmap/.claude/worktrees/harden-cli}"
DATA="${DATA:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb}"
PIXI_BIN="${PIXI_BIN:-$REPO/.pixi/envs/default/bin}"
QRNAMAP="${QRNAMAP:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/target-harden/release/qrnamap}"
GENOME_DIR="${GENOME_DIR:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data/genome}"
REF="$GENOME_DIR/GCF_000146045.2_R64_genomic.fna.gz"
GFF="$GENOME_DIR/GCF_000146045.2_R64_genomic.gff.gz"
SHEET="$DATA/samples.tsv"
CONTRAST="condition,AmpB,Control"

stage="${1:?stage}"; shift || true

case "$stage" in
e1-gen)
    export PATH="$PIXI_BIN:$PATH"
    "$QRNAMAP" run --samples "$DATA/e1/samples_missing.tsv" --reference "$REF" --gff "$GFF" \
        --mode quant --out "$DATA/e1/run_quant" --strandedness unstranded --threads 1 \
        --contrast "$CONTRAST" --slurm --samples-per-task 6 \
        --partition short --mem-gb 4 --time 00:10:00
    ;;
submit-clean)
    RUN_DIR="${1:?run dir}"
    # A shell with only the sbatch directory and the system directories on PATH.
    SB_DIR="$(dirname "$(command -v sbatch)")"
    env -i HOME="$HOME" USER="$USER" PATH="$SB_DIR:/usr/bin:/bin" \
        bash -c 'echo "submit shell PATH=$PATH"; command -v kallisto hisat2-align-s featureCounts Rscript || echo "(no pipeline tool on PATH)"; bash "$0"' \
        "$RUN_DIR/slurm/submit.sh"
    ;;
submit-login)
    # A fresh login shell (env -i, then bash -l): the site profile loads the default modules
    # (among them workspace/scratch, which the SLURM task prolog needs to set $SCRATCH in the
    # job), but the pixi env is not on PATH: no kallisto, hisat2 or featureCounts, and Rscript is
    # the cluster R module, not the pixi one.
    RUN_DIR="${1:?run dir}"
    env -i HOME="$HOME" USER="$USER" PATH="/usr/bin:/bin" bash -l -c '
        echo "submit shell PATH=$PATH"
        for t in kallisto hisat2-align-s hisat2-build-s featureCounts Rscript; do
            printf "  %s: %s\n" "$t" "$(command -v "$t" || echo NOT-FOUND)"
        done
        bash "$0"' "$RUN_DIR/slurm/submit.sh"
    ;;
strand)
    export PATH="$PIXI_BIN:$PATH"
    echo "node $(hostname) cpus ${SLURM_CPUS_ON_NODE:-?} start $(date -Is)"
    "$QRNAMAP" strand-check --samples "$SHEET" --reference "$REF" --gff "$GFF" \
        --out "$DATA/strand" --threads 4 --tmpdir "${SCRATCH:?}"
    echo "end $(date -Is)"
    ;;
slurm)
    MODE="${1:?mode}"; SPT="${2:?samples per task}"; STRAND="${3:?strandedness}"
    MEM="${4:?mem gb}"; TIME="${5:?time}"
    export PATH="$PIXI_BIN:$PATH"
    "$QRNAMAP" run --samples "$SHEET" --reference "$REF" --gff "$GFF" \
        --mode "$MODE" --out "$DATA/run_$MODE" --strandedness "$STRAND" --threads 1 \
        --contrast "$CONTRAST" --slurm --samples-per-task "$SPT" \
        --partition short --mem-gb "$MEM" --time "$TIME"
    ;;
compare)
    export PATH="$PIXI_BIN:$PATH"
    "$PIXI_BIN/Rscript" "$REPO/validation/amb/compare_modes.R" "$DATA/run_quant" "$DATA/run_align" \
        "$REPO/validation/amb/results"
    ;;
e3-manual)
    # The align steps of the CLI (crates/qrnamap-core/src/toolcall.rs hisat2_align and
    # feature_counts; crates/qrnamap-cli/src/pipeline.rs align_sample) run by hand, each tool
    # under /usr/bin/time -v. Uses the HISAT2 index, splice sites and SAF of $DATA/run_align.
    SAMPLE="${1:?sample}"; STRAND="${STRAND:-unstranded}"
    export PATH="$PIXI_BIN:$PATH"
    R="$DATA/run_align"
    W="${SCRATCH:?}/e3manual"; mkdir -p "$W"
    echo "node $(hostname) cpus ${SLURM_CPUS_ON_NODE:-?} job ${SLURM_JOB_ID:-?} start $(date -Is)"
    echo "== df -T SCRATCH"; df -T "${SCRATCH:?}"
    echo "== findmnt -T SCRATCH"; findmnt -T "${SCRATCH:?}"
    echo "== /proc/self/cgroup"; cat /proc/self/cgroup
    f1=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $2}' "$SHEET")
    f2=$(awk -F'\t' -v s="$SAMPLE" '$1==s{print $3}' "$SHEET")
    # Staging as the CLI does it: decompress each mate to a plain file in scratch.
    /usr/bin/time -v -o "$W/time_zcat1.txt" bash -c 'zcat "$0" > "$1"' "$f1" "$W/lane1_R1.fq"
    /usr/bin/time -v -o "$W/time_zcat2.txt" bash -c 'zcat "$0" > "$1"' "$f2" "$W/lane1_R2.fq"
    ss=()
    case "$STRAND" in
        reverse) ss=(--rna-strandness RF); fcs=2 ;;
        forward) ss=(--rna-strandness FR); fcs=1 ;;
        *) fcs=0 ;;
    esac
    /usr/bin/time -v -o "$W/time_hisat2.txt" hisat2-align-s -p 1 -x "$R/index/genome" \
        --known-splicesite-infile "$R/splice_sites.txt" --no-temp-splicesite --max-intronlen 5000 \
        "${ss[@]}" -1 "$W/lane1_R1.fq" -2 "$W/lane1_R2.fq" -S "$W/$SAMPLE.sam" 2> "$W/hisat2.stderr"
    cat "$W/hisat2.stderr"
    /usr/bin/time -v -o "$W/time_fc.txt" featureCounts -T 1 -F SAF -a "$R/genes.saf" -p --countReadPairs \
        -s "$fcs" -o "$W/featurecounts.tsv" "$W/$SAMPLE.sam" > "$W/fc.stdout" 2>&1
    tail -5 "$W/featurecounts.tsv.summary" || true
    echo "== byte sizes"
    stat -c '%s %n' "$W/lane1_R1.fq" "$W/lane1_R2.fq" "$W/$SAMPLE.sam"
    for t in zcat1 zcat2 hisat2 fc; do
        echo "== /usr/bin/time -v $t"; grep -E 'Maximum resident|Elapsed|Exit status|File system outputs' "$W/time_$t.txt"
    done
    mkdir -p "$DATA/e3/manual_${SLURM_JOB_ID:-x}"
    cp "$W"/time_*.txt "$W/hisat2.stderr" "$W/featurecounts.tsv.summary" "$DATA/e3/manual_${SLURM_JOB_ID:-x}/"
    echo "== memory accounting at the end of the job (files still in scratch)"
    cg=$(awk -F: '$2 ~ /memory/ || $1 == "0" {print $3; exit}' /proc/self/cgroup)
    for base in /sys/fs/cgroup/memory /sys/fs/cgroup; do
        if [ -r "$base$cg/memory.stat" ]; then
            echo "memory.stat: $base$cg/memory.stat"; cat "$base$cg/memory.stat"
            for f in memory.max_usage_in_bytes memory.usage_in_bytes memory.peak memory.current memory.limit_in_bytes memory.max; do
                [ -r "$base$cg/$f" ] && echo "$f $(cat "$base$cg/$f")"
            done
        fi
    done
    # Also the job-level cgroup (v1: .../job_ID), when the step cgroup is below it.
    jc=$(echo "$cg" | sed -n 's#^\(.*/job_[0-9]*\)/.*#\1#p')
    if [ -n "$jc" ] && [ -r "/sys/fs/cgroup/memory$jc/memory.stat" ]; then
        echo "job cgroup memory.stat: /sys/fs/cgroup/memory$jc/memory.stat"; cat "/sys/fs/cgroup/memory$jc/memory.stat"
        for f in memory.max_usage_in_bytes memory.usage_in_bytes memory.limit_in_bytes; do
            [ -r "/sys/fs/cgroup/memory$jc/$f" ] && echo "job $f $(cat "/sys/fs/cgroup/memory$jc/$f")"
        done
    fi
    echo "end $(date -Is)"
    ;;
e3-align)
    # The whole `qrnamap align` for one sample, with --tmpdir on node-local scratch or on GPFS.
    SAMPLE="${1:?sample}"; KIND="${2:?scratch|gpfs}"; STRAND="${STRAND:-unstranded}"
    export PATH="$PIXI_BIN:$PATH"
    O="$DATA/e3/align_${KIND}_${SLURM_JOB_ID:?}"; mkdir -p "$O"
    awk -F'\t' -v s="$SAMPLE" 'NR==1 || $1==s' "$SHEET" > "$O/sheet.tsv"
    case "$KIND" in
        scratch) T="${SCRATCH:?}/e3tmp" ;;
        gpfs) T="$DATA/e3/tmp_gpfs_${SLURM_JOB_ID:?}" ;;
        *) echo "KIND must be scratch or gpfs" >&2; exit 2 ;;
    esac
    mkdir -p "$T"
    echo "node $(hostname) cpus ${SLURM_CPUS_ON_NODE:-?} job $SLURM_JOB_ID kind $KIND tmpdir $T start $(date -Is)"
    df -T "$T"
    # Same work as one array task: reuse the prepared index of $DATA/run_align (copied, with
    # its .prepared marker), so no index build enters this job's MaxRSS.
    cp -rp "$DATA/run_align/index" "$DATA/run_align/splice_sites.txt" "$DATA/run_align/genes.saf" \
        "$DATA/run_align/.prepared" "$O/"
    "$QRNAMAP" align --samples "$O/sheet.tsv" --reference "$REF" --gff "$GFF" --out "$O" \
        --strandedness "$STRAND" --threads 1 --tmpdir "$T" --task-id 0 --samples-per-task 1
    rm -rf "$T" "$O/index"
    echo "end $(date -Is)"
    ;;
*)
    echo "unknown stage: $stage" >&2; exit 2 ;;
esac
