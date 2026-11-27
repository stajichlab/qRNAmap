#!/usr/bin/env bash
# Page DE vs DESeq2 on one count matrix (HPCC; pixi environment of the repo for R/DESeq2, node for the page code).
#   bash validation/de/run_compare.sh COUNTS COLDATA FACTOR LEVEL REF OUTDIR ["~ condition" "~ strain + condition" ...]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
counts=$1 coldata=$2 factor=$3 level=$4 ref=$5 out=$6; shift 6
designs=("$@"); [[ ${#designs[@]} -gt 0 ]] || designs=("~ $factor")
mkdir -p "$out"
node "$HERE/page_de.mjs" "$counts" "$coldata" "$factor" "$level" "$ref" "$out/page.tsv"
i=0
for d in "${designs[@]}"; do
  i=$((i + 1))
  echo "== DESeq2 design: $d"
  pixi run --manifest-path "$HERE/../../pixi.toml" -e de Rscript "$HERE/deseq2_reference.R" "$counts" "$coldata" "$d" "$factor" "$level" "$ref" "$out/deseq2_$i"
  echo "== page vs DESeq2 ($d)"
  python3 "$HERE/compare_de.py" "$out/deseq2_$i.results.tsv" "$out/page.tsv"
done
