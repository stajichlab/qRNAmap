#!/usr/bin/env bash
# Build pre-computed kallisto index + annotations for a genome accession.
# Outputs:
#   OUTDIR/<acc>.idx.gz
#   OUTDIR/<acc>.tx2gene.tsv.gz
#   OUTDIR/<acc>.gene_info.tsv.gz
#   OUTDIR/<acc>.genes.tsv.gz
# INDEX_FROM=<URL of an existing <acc>.idx.gz> keeps that index instead of building one (to refresh the annotation tables
# of a published set, for example after gene_info.tsv gained its description column). The new tx2gene must equal the one
# next to that index (<acc>.tx2gene.tsv.gz at the same place), else the index and the tables would not match: the run stops.
set -euo pipefail

ACCESSION="${1:?Usage: build-one-index.sh ACCESSION FNA_URL GFF_URL OUTDIR}"
FNA_URL="${2:?missing FNA_URL}"
GFF_URL="${3:?missing GFF_URL}"
OUTDIR="${4:?missing OUTDIR}"

TMPDIR="${SCRATCH:-/tmp}/qidx_${ACCESSION}_$$"
mkdir -p "$TMPDIR" "$OUTDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== Processing $ACCESSION ==="
echo "Downloading genome FASTA..."
curl -sSL --retry 3 --max-time 300 "$FNA_URL" -o "$TMPDIR/genome.fna.gz"
echo "Downloading genome GFF3..."
curl -sSL --retry 3 --max-time 300 "$GFF_URL" -o "$TMPDIR/genome.gff.gz"

echo "Building transcripts with qrnamap..."
qrnamap transcripts \
  --reference "$TMPDIR/genome.fna.gz" \
  --gff "$TMPDIR/genome.gff.gz" \
  --out "$TMPDIR/transcripts.fa" \
  --tx2gene "$TMPDIR/tx2gene.tsv" \
  --gene-info "$TMPDIR/gene_info.tsv" \
  --genes-table "$TMPDIR/genes.tsv"

# The index is built with the kallisto.wasm build of the page (scripts/build-index-wasm.mjs), not a native kallisto:
# index files differ between kallisto builds of the same version (bench 2026-11-26_index-compat: an index from the Linux
# conda kallisto 0.52.0 is "Corrupted" for kallisto.wasm and for macOS Homebrew kallisto 0.52.0).
# INDEXER=native uses `kallisto index` instead (for a native-only set).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${INDEX_FROM:-}" ]]; then
  echo "Keeping the index $INDEX_FROM (no index build)..."
  curl -sSL --fail --retry 3 --max-time 600 "$INDEX_FROM" -o "$TMPDIR/kept.idx.gz"
  gzip -t "$TMPDIR/kept.idx.gz"
  curl -sSL --fail --retry 3 --max-time 300 "${INDEX_FROM%.idx.gz}.tx2gene.tsv.gz" -o "$TMPDIR/kept.tx2gene.tsv.gz"
  if ! gzip -dc "$TMPDIR/kept.tx2gene.tsv.gz" | cmp -s - "$TMPDIR/tx2gene.tsv"; then
    echo "the new tx2gene differs from the one of the kept index: rebuild the index instead (unset INDEX_FROM)" >&2
    exit 1
  fi
  echo "tx2gene equals the kept index's"
elif [[ "${INDEXER:-wasm}" == native ]]; then
  echo "Building kallisto index (k=31, native)..."
  kallisto index -i "$TMPDIR/transcripts.idx" -k 31 "$TMPDIR/transcripts.fa"
else
  echo "Building kallisto index (k=31, wasm)..."
  node "$HERE/build-index-wasm.mjs" "$TMPDIR/transcripts.fa" "$TMPDIR/transcripts.idx"
  # Self-check in kallisto.wasm: 2,000 single-end reads of 100 bases cut from the transcripts themselves must load and
  # mostly pseudoalign (a genome-independent check; the page's fixture reads are yeast).
  awk 'BEGIN { q = sprintf("%100s", ""); gsub(/ /, "I", q) }
       function emit() { if (length(seq) >= 150 && n < 2000) { print "@r" n; print substr(seq, 26, 100); print "+"; print q; n++ } }
       /^>/ { emit(); seq = ""; next } { seq = seq $0 } END { emit() }' "$TMPDIR/transcripts.fa" > "$TMPDIR/selfcheck.fq"
  check=$(node "$HERE/check-index-wasm.mjs" "$TMPDIR/transcripts.idx" "$TMPDIR/selfcheck.fq")
  echo "self-check: $check"
  read -r status processed aligned <<<"$check"
  [[ "$status" == ok && "$processed" -gt 0 && $(( aligned * 10 )) -ge $(( processed * 9 )) ]] || { echo "index self-check failed" >&2; exit 1; }
fi

echo "Compressing deliverables with gzip -9..."
if [[ -n "${INDEX_FROM:-}" ]]; then
  cp "$TMPDIR/kept.idx.gz" "$OUTDIR/${ACCESSION}.idx.gz" # byte for byte: the published index stays as it is
else
  gzip -9 -c "$TMPDIR/transcripts.idx" > "$OUTDIR/${ACCESSION}.idx.gz"
fi
gzip -9 -c "$TMPDIR/tx2gene.tsv" > "$OUTDIR/${ACCESSION}.tx2gene.tsv.gz"
gzip -9 -c "$TMPDIR/gene_info.tsv" > "$OUTDIR/${ACCESSION}.gene_info.tsv.gz"
gzip -9 -c "$TMPDIR/genes.tsv" > "$OUTDIR/${ACCESSION}.genes.tsv.gz"

# Compute sha256 checksums (one file per accession: the release step joins them; a shared name was overwritten when
# the build artifacts were merged, which left only the last genome in checksums.sha256 of indexes-v1)
(cd "$OUTDIR" && sha256sum \
  "${ACCESSION}.idx.gz" \
  "${ACCESSION}.tx2gene.tsv.gz" \
  "${ACCESSION}.gene_info.tsv.gz" \
  "${ACCESSION}.genes.tsv.gz" > "${ACCESSION}.sha256")

echo "=== Successfully built index for $ACCESSION ==="
ls -lh "$OUTDIR/${ACCESSION}"*
