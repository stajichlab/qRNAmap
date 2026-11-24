#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=qrnamap_expected
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/qRNAmap_data/page-v1-cli-expected/slurm_%j.out

set -euo pipefail

PROJ_ROOT="/bigdata/stajichlab/jstajich/projects/qRNAmap/.claude/worktrees/page-v1"
DATA_ROOT="/bigdata/stajichlab/jstajich/projects/qRNAmap_data"
BIN="$PROJ_ROOT/target/debug/qrnamap"

export PATH="$PROJ_ROOT/.pixi/envs/default/bin:$PATH"
export R_LIBS_USER="$PROJ_ROOT/.pixi/envs/default/lib/R/library"

# 1. Yeast AmpB
YEAST_OUT="$DATA_ROOT/page-v1-cli-expected/yeast-ampb"
mkdir -p "$YEAST_OUT"
"$BIN" run \
  --samples "$PROJ_ROOT/validation/browser-page/expected/yeast-ampb/samples.tsv" \
  --reference "$DATA_ROOT/spike0/fixtures/genomes/yeast/GCF_000146045.2_R64_genomic.fna.gz" \
  --gff "$DATA_ROOT/spike0/fixtures/genomes/yeast/GCF_000146045.2_R64_genomic.gff.gz" \
  --mode quant \
  --quantifier kallisto \
  --threads 1 \
  --design "~ condition" \
  --contrast "condition,AmpB,control" \
  --out "$YEAST_OUT"

# 2. C. neoformans H99 Lkh1
CNEO_OUT="$DATA_ROOT/page-v1-cli-expected/cneo-lkh1"
mkdir -p "$CNEO_OUT"
"$BIN" run \
  --samples "$PROJ_ROOT/validation/browser-page/expected/cneo-lkh1/samples.tsv" \
  --reference "$DATA_ROOT/page-v1/dl_h99/GCF_000149245.1_CNA3_genomic.fna.gz" \
  --gff "$DATA_ROOT/page-v1/dl_h99/GCF_000149245.1_CNA3_genomic.gff.gz" \
  --mode quant \
  --quantifier kallisto \
  --threads 1 \
  --design "~ genotype" \
  --contrast "genotype,lkh1_delta,wild_type" \
  --out "$CNEO_OUT"

# Decompress counts.tsv.zst to counts.tsv and copy to validation/browser-page/expected/
zstd -dc "$YEAST_OUT/counts.tsv.zst" > "$PROJ_ROOT/validation/browser-page/expected/yeast-ampb/counts.tsv"
zstd -dc "$CNEO_OUT/counts.tsv.zst" > "$PROJ_ROOT/validation/browser-page/expected/cneo-lkh1/counts.tsv"

# Copy mapping summaries and DESeq2 outputs
cp "$YEAST_OUT/mapping_summary.tsv.zst" "$PROJ_ROOT/validation/browser-page/expected/yeast-ampb/"
cp "$CNEO_OUT/mapping_summary.tsv.zst" "$PROJ_ROOT/validation/browser-page/expected/cneo-lkh1/"
if [ -f "$CNEO_OUT/deseq2.results.csv" ]; then
  cp "$CNEO_OUT/deseq2.results.csv" "$PROJ_ROOT/validation/browser-page/expected/cneo-lkh1/"
fi

# Hash counts.tsv
(cd "$PROJ_ROOT/validation/browser-page/expected" && sha256sum yeast-ampb/counts.tsv cneo-lkh1/counts.tsv > counts_sha256.tsv)

echo "ALL DONE"
