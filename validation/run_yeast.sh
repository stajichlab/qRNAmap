#!/usr/bin/env bash
# run_yeast.sh: Task 14 validation driver (qRNAmap CLI on public yeast RNA-seq, UCR HPCC).
#
# Usage: bash /ABS/PATH/run_yeast.sh STAGE [ARGS]
# All paths are absolute and come from the variables below (or the environment).
# This script never uses BASH_SOURCE (it runs inside SLURM jobs).
#
# Stages (in order):
#   prefetch            login node: prefetch every run in $DATA/ena/chosen.tsv to $DATA/sra
#   dump RUN PAIRS      SLURM job body: parallel-fastq-dump, check both mates == PAIRS, delete the .sra
#   submit-dump         login node: sbatch one dump job per run
#   genome              login node: qrnamap fetch (NCBI download only)
#   parity              SLURM job body: qrnamap transcripts + gffread -w + compare_gffread.py
#   sheet               write $DATA/samples.tsv
#   slurm MODE [K] [STRAND]  qrnamap run --slurm for MODE (quant|align), K samples per task
#   submit MODE         login node: run the generated submit.sh with the pixi env on PATH
#   strand              SLURM job body: qrnamap strand-check on every sample
#   compare             SLURM job body: compare_modes (Spearman, DE overlap) with R from the pixi env
set -euo pipefail

REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/qRNAmap/.claude/worktrees/cli-core}"
DATA="${DATA:-/bigdata/stajichlab/jstajich/projects/qRNAmap_data}"
PIXI_BIN="$REPO/.pixi/envs/default/bin"
QRNAMAP="$REPO/target/release/qrnamap"
GENOME_DIR="$DATA/genome"
REF="$GENOME_DIR/GCF_000146045.2_R64_genomic.fna.gz"
GFF="$GENOME_DIR/GCF_000146045.2_R64_genomic.gff.gz"
PARTITION="${PARTITION:-short}"
CONTRAST="condition,ups1KO,WT"

stage="${1:?stage}"; shift || true

case "$stage" in
prefetch)
    source /etc/profile.d/modules.sh
    module load sratoolkit/3.2.0
    mkdir -p "$DATA/sra" "$DATA/logs"
    for r in $(cut -f1 "$DATA/ena/chosen.tsv"); do
        prefetch "$r" -O "$DATA/sra" --max-size 30G > "$DATA/logs/prefetch_$r.log" 2>&1
    done
    ;;
dump)
    RUN="${1:?run}"; EXPECT="${2:?expected pairs}"
    source /etc/profile.d/modules.sh
    module load sratoolkit/3.2.0
    module load parallel-fastq-dump/0.6.7
    mkdir -p "$DATA/fastq"
    # $SCRATCH is node-local: temporary files only. The FASTQ output goes to /bigdata.
    parallel-fastq-dump --sra-id "$DATA/sra/$RUN/$RUN.sra" --threads 8 --outdir "$DATA/fastq" \
        --split-files --gzip --tmpdir "${SCRATCH:?}"
    for m in 1 2; do
        n=$(( $(zcat "$DATA/fastq/${RUN}_$m.fastq.gz" | wc -l) / 4 ))
        echo "${RUN}_$m: $n reads (expected $EXPECT)"
        [ "$n" -eq "$EXPECT" ] || { echo "ERROR: read count mismatch for ${RUN}_$m" >&2; exit 1; }
    done
    rm -f "$DATA/sra/$RUN/$RUN.sra"
    rmdir "$DATA/sra/$RUN" 2>/dev/null || true
    echo "DUMP_OK $RUN"
    ;;
submit-dump)
    mkdir -p "$DATA/logs"
    while IFS=$'\t' read -r run pairs _; do
        sbatch -p "$PARTITION" -N 1 -n 1 -c 8 --mem 16G -t 02:00:00 --job-name "dump_$run" \
            -o "$DATA/logs/dump_$run.%j.out" -e "$DATA/logs/dump_$run.%j.err" \
            --wrap "bash $REPO/validation/run_yeast.sh dump $run $pairs"
    done < "$DATA/ena/chosen.tsv"
    ;;
genome)
    "$QRNAMAP" fetch --source ncbi --id GCF_000146045.2 --out "$GENOME_DIR"
    ;;
parity)
    export PATH="$PIXI_BIN:$PATH"
    P="$DATA/parity"; mkdir -p "$P"
    T="${SCRATCH:?}/parity"; mkdir -p "$T"
    "$QRNAMAP" transcripts --reference "$REF" --gff "$GFF" \
        --out "$P/qrnamap_transcripts.fa.zst" --tx2gene "$P/tx2gene.tsv.zst"
    # gffread needs a plain, indexable FASTA and reads the GFF3 from a plain file.
    zcat "$REF" > "$T/genome.fa"
    zcat "$GFF" > "$T/genes.gff3"
    gffread -w "$T/gffread_transcripts.fa" -g "$T/genome.fa" "$T/genes.gff3"
    zstd -q -f -T0 "$T/gffread_transcripts.fa" -o "$P/gffread_transcripts.fa.zst"
    "$PIXI_BIN/python" -I "$REPO/validation/compare_gffread.py" \
        --ours "$P/qrnamap_transcripts.fa.zst" --gffread "$P/gffread_transcripts.fa.zst" \
        --gff "$GFF" --out "$P/parity_diff.tsv"
    ;;
sheet)
    S="$DATA/samples.tsv"
    printf 'sample\tfastq1\tfastq2\tcondition\treplicate\n' > "$S"
    while IFS=$'\t' read -r run _ _ _ title _; do
        case "$title" in
            *WT-___*) cond=WT ;;
            *KO-___-ups1*) cond=ups1KO ;;
            *) echo "unknown title $title" >&2; exit 1 ;;
        esac
        rep="${title##*REP_}"
        printf '%s\t%s\t%s\t%s\tr%s\n' "${cond}_$rep" "$DATA/fastq/${run}_1.fastq.gz" \
            "$DATA/fastq/${run}_2.fastq.gz" "$cond" "$rep" >> "$S"
    done < "$DATA/ena/chosen.tsv"
    cat "$S"
    ;;
slurm)
    MODE="${1:?mode}"; SPT="${2:-1}"; STRAND="${3:-reverse}"
    export PATH="$PIXI_BIN:$PATH"
    "$QRNAMAP" run --samples "$DATA/samples.tsv" --reference "$REF" --gff "$GFF" \
        --mode "$MODE" --out "$DATA/run_$MODE" --strandedness "$STRAND" --threads 1 \
        --contrast "$CONTRAST" --slurm --samples-per-task "$SPT" \
        --partition "$PARTITION" --mem-gb 16 --time 02:00:00
    ;;
submit)
    # The generated scripts call kallisto, hisat2, featureCounts and Rscript by name. sbatch
    # exports the submitting shell's environment, so the pixi env must be on PATH here.
    MODE="${1:?mode}"
    export PATH="$PIXI_BIN:$PATH"
    bash "$DATA/run_$MODE/slurm/submit.sh"
    ;;
strand)
    export PATH="$PIXI_BIN:$PATH"
    "$QRNAMAP" strand-check --samples "$DATA/samples.tsv" --reference "$REF" --gff "$GFF" \
        --out "$DATA/strand" --threads 4
    ;;
compare)
    export PATH="$PIXI_BIN:$PATH"
    "$PIXI_BIN/Rscript" "$REPO/validation/compare_modes.R" "$DATA/run_quant" "$DATA/run_align" \
        "$REPO/validation/results"
    ;;
*)
    echo "unknown stage: $stage" >&2; exit 2 ;;
esac
