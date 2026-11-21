#!/usr/bin/env bash
# Browser spike 0, Task 2: build the fixtures. No BASH_SOURCE; all paths come from the environment.
# Usage: fixtures.sh STAGE      STAGE = fetch (needs internet: run on the login node)
#                                     | build (reads, transcripts, count matrix: run under SLURM)
# Required environment (absolute paths):
#   FIX         fixture directory (outputs)
#   QRNAMAP     path of the qrnamap binary (built from main)
#   PIXI_BIN    <worktree>/.pixi/envs/default/bin
#   SCRIPTS     directory with subset_reads.py and merge_counts.py (validation/browser-spike0)
#   AMB         /bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb
set -euo pipefail
: "${FIX:?}" "${QRNAMAP:?}" "${PIXI_BIN:?}" "${SCRIPTS:?}" "${AMB:?}"
export PATH="$PIXI_BIN:$PATH"
YEAST=GCF_000146045.2
AFUM=GCF_000002655.1
SRC_SAMPLE=AmpB_1
SRC_RUN=SRR3396381
SRC_PAIRS=8405534      # ENA read_count of SRR3396381 (amb/results/ena_runs.tsv); checked by subset_reads.py
N_PAIRS=200000
SEED=20261121
READS_DIR=/bigdata/gen220/shared/data-examples/rnaseq/kallisto/data

case "${1:?stage: fetch|build}" in
fetch)
    mkdir -p "$FIX/genomes/yeast" "$FIX/genomes/afum"
    "$QRNAMAP" fetch --source ncbi --id "$YEAST" --out "$FIX/genomes/yeast"
    "$QRNAMAP" fetch --source ncbi --id "$AFUM" --out "$FIX/genomes/afum"
    ls -l "$FIX/genomes/yeast" "$FIX/genomes/afum"
    ;;
build)
    mkdir -p "$FIX/reads" "$FIX/transcripts" "$FIX/counts"
    echo "ENA check of $SRC_RUN:"
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRC_RUN&result=read_run&fields=run_accession,read_count,base_count,library_layout&format=tsv" \
        | tee "$FIX/reads/ena_$SRC_RUN.tsv" || echo "(ENA not reachable from this node; the count in amb/results/ena_runs.tsv is used)"
    python3 -I "$SCRIPTS/subset_reads.py" "$READS_DIR/${SRC_RUN}_1.fastq.gz" "$READS_DIR/${SRC_RUN}_2.fastq.gz" \
        "$SRC_PAIRS" "$N_PAIRS" "$SEED" "$FIX/reads/yeast_sub_1.fastq.gz" "$FIX/reads/yeast_sub_2.fastq.gz"
    for g in yeast afum; do
        fa=$(ls "$FIX/genomes/$g"/*_genomic.fna.gz); gff=$(ls "$FIX/genomes/$g"/*_genomic.gff.gz)
        "$QRNAMAP" transcripts --reference "$fa" --gff "$gff" \
            --out "$FIX/transcripts/$g.transcripts.fa.zst" --tx2gene "$FIX/transcripts/$g.tx2gene.tsv.zst" \
            --saf "$FIX/transcripts/$g.genes.saf" --splice-sites "$FIX/transcripts/$g.splice_sites.txt"
    done
    # AMB gene-count matrix: the six kallisto gene-count tables of the AMB run (amb/run_quant).
    python3 -I "$SCRIPTS/merge_counts.py" "$AMB/run_quant" "$FIX/counts/amb_gene_counts.tsv" \
        Control_1 Control_2 Control_3 AmpB_1 AmpB_2 AmpB_3
    ;;
esac
