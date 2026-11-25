#!/usr/bin/env bash
# Build pre-computed kallisto index + annotations for a genome accession.
# Outputs:
#   OUTDIR/<acc>.idx.gz
#   OUTDIR/<acc>.tx2gene.tsv.gz
#   OUTDIR/<acc>.gene_info.tsv.gz
#   OUTDIR/<acc>.genes.tsv.gz
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

echo "Building kallisto index (k=31)..."
kallisto index -i "$TMPDIR/transcripts.idx" -k 31 "$TMPDIR/transcripts.fa"

echo "Compressing deliverables with gzip -9..."
gzip -9 -c "$TMPDIR/transcripts.idx" > "$OUTDIR/${ACCESSION}.idx.gz"
gzip -9 -c "$TMPDIR/tx2gene.tsv" > "$OUTDIR/${ACCESSION}.tx2gene.tsv.gz"
gzip -9 -c "$TMPDIR/gene_info.tsv" > "$OUTDIR/${ACCESSION}.gene_info.tsv.gz"
gzip -9 -c "$TMPDIR/genes.tsv" > "$OUTDIR/${ACCESSION}.genes.tsv.gz"

# Compute sha256 checksums
(cd "$OUTDIR" && sha256sum \
  "${ACCESSION}.idx.gz" \
  "${ACCESSION}.tx2gene.tsv.gz" \
  "${ACCESSION}.gene_info.tsv.gz" \
  "${ACCESSION}.genes.tsv.gz" >> checksums.sha256)

echo "=== Successfully built index for $ACCESSION ==="
ls -lh "$OUTDIR/${ACCESSION}"*
