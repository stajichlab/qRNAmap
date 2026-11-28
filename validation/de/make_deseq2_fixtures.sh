#!/usr/bin/env bash
# Phase 1 of the DESeq2 port (docs/superpowers/specs/2026-11-27-deseq2-rust-design.md): fetch the three datasets, build
# the fixture inputs, run R DESeq2 on the full data, add the genes it flagged to the fixtures, run DESeq2 with every intermediate saved (deseq2_dump.R), make the m - p <= 3 prior table and
# check it against DESeq2, then copy fixtures and table into the crate. Runs on the HPCC (in a job).
#   bash validation/de/make_deseq2_fixtures.sh WORK_DIR
# WORK_DIR (one-off, purgeable): raw/ downloads, inputs/, dumps/fixtures/ and dumps/full/ (full-data reference runs,
# not committed), prior_hist.tsv, summary.tsv. Writes into the repo: crates/qrnamap-core/tests/data/deseq/<fixture>/
# and crates/qrnamap-core/src/deseq/prior_hist.tsv.
# RSCRIPT: the Rscript to use (default: pixi environment `de` of the repo).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
WORK="${1:?usage: make_deseq2_fixtures.sh WORK_DIR}"
read -r -a R <<< "${RSCRIPT:-pixi run --manifest-path $ROOT/pixi.toml -e de Rscript}"
RAW="$WORK/raw"
mkdir -p "$RAW/yeast" "$RAW/pasilla" "$RAW/afum"

# Downloads (kept if present). Yeast: Schurch et al. 2016, MIT; pasilla: Bioconductor data package, LGPL.
G=https://github.com/bartongroup/profDGE48/raw/master
if [[ -z "$(ls "$RAW"/yeast/*.gbgout 2>/dev/null)" ]]; then
  for f in WT_countdata.tar.gz Snf2_countdata.tar.gz; do
    curl -sSfL --retry 3 -o "$RAW/$f" "$G/Preprocessed_data/$f"
    tar -xzf "$RAW/$f" -C "$RAW/yeast"
  done
fi
[[ -s "$RAW/exclude.lst" ]] || curl -sSfL --retry 3 -o "$RAW/exclude.lst" "$G/Bad_replicate_identification/exclude.lst"
if [[ ! -s "$RAW/pasilla/pasilla_gene_counts.tsv" ]]; then
  curl -sSfL --retry 3 -o "$RAW/pasilla.tar.gz" \
    https://bioconductor.org/packages/release/data/experiment/src/contrib/pasilla_1.42.0.tar.gz
  tar -xzf "$RAW/pasilla.tar.gz" -C "$RAW/pasilla" --strip-components=3 \
    pasilla/inst/extdata/pasilla_gene_counts.tsv pasilla/inst/extdata/pasilla_sample_annotation.csv
fi
gzip -dc "$HERE/data/afum_azole_gmm.counts.tsv.gz" > "$RAW/afum/counts.full.tsv"
cp "$HERE/data/afum_azole_gmm.coldata.tsv" "$RAW/afum/coldata.tsv"

"${R[@]}" "$HERE/build_deseq2_fixtures.R" "$RAW" "$WORK/inputs"

dump() { # INPUT_DIR OUT_DIR
  local spec; spec="$1/spec.tsv"
  val() { awk -F'\t' -v k="$1" '$1 == k { print $2 }' "$spec"; }
  "${R[@]}" "$HERE/deseq2_dump.R" "$1/counts.tsv" "$1/coldata.tsv" "$(val design)" "$(val contrast)" "$(val ref)" "$2"
}
for d in "$WORK"/inputs/full/*; do dump "$d" "$WORK/dumps/full/$(basename "$d")"; done
for d in "$WORK"/inputs/fixtures/*; do
  n=$(basename "$d")
  if [[ -d "$WORK/dumps/full/$n" ]]; then
    "${R[@]}" "$HERE/add_flagged_genes.R" "$WORK/dumps/full/$n" "$WORK/inputs/full/$n/counts.tsv" "$d/counts.tsv"
  fi
  dump "$d" "$WORK/dumps/fixtures/$n"
done

"${R[@]}" "$HERE/deseq2_prior_table.R" "$WORK/prior_hist.tsv"
"${R[@]}" "$HERE/check_prior_table.R" "$WORK/prior_hist.tsv" "$WORK"/dumps/fixtures/* "$WORK"/dumps/full/*

# One line per run: the cases each dump covers.
keys=(m p trend dispPriorVar_path genes genes_all_zero n_disp_outlier n_gene_grid n_optim_rows n_cooks_flagged n_replaced n_filtered filter_theta n_padj_0.05)
{ printf 'kind\tname'; printf '\t%s' "${keys[@]}"; echo
  for d in "$WORK"/dumps/*/*; do
    printf '%s\t%s' "$(basename "$(dirname "$d")")" "$(basename "$d")"
    for k in "${keys[@]}"; do printf '\t%s' "$(awk -F'\t' -v k="$k" '$1 == k { print $2 }' "$d/run.tsv")"; done; echo
  done; } > "$WORK/summary.tsv"
column -t -s $'\t' "$WORK/summary.tsv"

DEST="$ROOT/crates/qrnamap-core/tests/data/deseq"
for d in "$WORK"/dumps/fixtures/*; do
  n=$(basename "$d"); mkdir -p "$DEST/$n"
  cp "$WORK/inputs/fixtures/$n"/{counts.tsv,coldata.tsv,spec.tsv} "$d"/* "$DEST/$n/"
done
mkdir -p "$ROOT/crates/qrnamap-core/src/deseq"
cp "$WORK/prior_hist.tsv" "$ROOT/crates/qrnamap-core/src/deseq/prior_hist.tsv"
du -sh "$DEST"
