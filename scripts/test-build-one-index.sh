#!/usr/bin/env bash
# Fixture test of build-one-index.sh with INDEX_FROM (keep a published index, refresh the annotation tables). No network:
# the genome and the "published" set are file:// URLs to the yeast mitochondrial fixtures.
#   bash scripts/test-build-one-index.sh        (needs target/release/qrnamap: cargo build --release -p qrnamap-cli)
# Cases: a published set whose tx2gene matches keeps its index byte for byte and gets a gene_info.tsv with the
# description column; a published tx2gene that differs stops the run with no output for that genome.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
FIX="$ROOT/validation/browser-page/fixtures"
export PATH="$ROOT/target/release:$PATH"
command -v qrnamap >/dev/null || { echo "qrnamap not built: cargo build --release -p qrnamap-cli" >&2; exit 1; }
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
export SCRATCH="$T/scratch"
mkdir -p "$SCRATCH" "$T/pub"
fail() { echo "FAIL: $*" >&2; exit 1; }

# The "published" set: the fixture's tx2gene and a stand-in index.
qrnamap transcripts --reference "$FIX/mito.fna.gz" --gff "$FIX/mito.gff3.gz" --out "$T/tx.fa" --tx2gene "$T/tx2gene.tsv" >/dev/null 2>&1
gzip -c "$T/tx2gene.tsv" > "$T/pub/GCF_T.tx2gene.tsv.gz"
echo "published index bytes" | gzip -c > "$T/pub/GCF_T.idx.gz"

# 1. Matching tx2gene: the index is kept as it is, the tables are new.
INDEX_FROM="file://$T/pub/GCF_T.idx.gz" bash "$HERE/build-one-index.sh" GCF_T "file://$FIX/mito.fna.gz" "file://$FIX/mito.gff3.gz" "$T/out" > "$T/log1" 2>&1 \
  || { cat "$T/log1"; fail "the refresh with a matching tx2gene failed"; }
cmp -s "$T/pub/GCF_T.idx.gz" "$T/out/GCF_T.idx.gz" || fail "the kept index changed"
gzip -dc "$T/out/GCF_T.gene_info.tsv.gz" | head -1 | grep -q $'\tdescription$' || fail "gene_info.tsv has no description column"
for f in GCF_T.idx.gz GCF_T.tx2gene.tsv.gz GCF_T.gene_info.tsv.gz GCF_T.genes.tsv.gz; do
  grep -q "  $f$" "$T/out/GCF_T.sha256" || fail "$f is not in GCF_T.sha256"
done
(cd "$T/out" && sha256sum --check --quiet GCF_T.sha256) || fail "GCF_T.sha256 does not match the files"
grep -q "tx2gene equals the kept index's" "$T/log1" || fail "no tx2gene check in the log"
echo "ok 1 - a matching published tx2gene keeps the index and refreshes the tables"

# 2. A published tx2gene that differs: the run stops and writes nothing for the genome.
printf 'rna-x\tgene-x\n' | gzip -c > "$T/pub/GCF_T.tx2gene.tsv.gz"
if INDEX_FROM="file://$T/pub/GCF_T.idx.gz" bash "$HERE/build-one-index.sh" GCF_T "file://$FIX/mito.fna.gz" "file://$FIX/mito.gff3.gz" "$T/out2" > "$T/log2" 2>&1; then
  fail "a different published tx2gene was accepted"
fi
grep -q "differs from the one of the kept index" "$T/log2" || { cat "$T/log2"; fail "no tx2gene message"; }
[[ -z "$(ls -A "$T/out2" 2>/dev/null)" ]] || fail "files were written for a refused genome"
echo "ok 2 - a different published tx2gene stops the run"
echo "all passed"
