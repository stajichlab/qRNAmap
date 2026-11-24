#!/usr/bin/env bash
# Page-v1 measurements, CI only: fetch the reads and the yeast transcripts as static files from GitHub releases.
# Usage: fetch_reads_release.sh <out dir> <prefix|full>. Run from the repository root. No BASH_SOURCE is used.
#   prefix: the first 1,000,000 read pairs of each mate of SRR3396381 (for counts <= 1,000,000).
#   full:   both complete mates of SRR3396381 (8,405,534 pairs each).
# The files come from release page-v1-test-data-v1 of stajichlab/qRNAmap. They are saved as SRR3396381_1.fastq.gz and
# SRR3396381_2.fastq.gz in <out dir>. The SHA-256 of each file is checked against validation/browser-page/measure/reads_prefix.tsv
# (committed in the repository), not against a file of the release. A replaced asset fails the check.
# The transcripts come from the spike 0 fixture release; the SHA-256 is checked against validation/browser-spike0/manifest.tsv.
# prepare_reads.sh is the documented way to re-create the reads from NCBI SRA.
set -euo pipefail
out="${1:?usage: fetch_reads_release.sh <out dir> <prefix|full>}"
kind="${2:?usage: fetch_reads_release.sh <out dir> <prefix|full>}"
case "$kind" in prefix|full) ;; *) echo "kind must be prefix or full" >&2; exit 2 ;; esac
run=SRR3396381
rel=https://github.com/stajichlab/qRNAmap/releases/download/page-v1-test-data-v1
table=validation/browser-page/measure/reads_prefix.tsv
mkdir -p "$out"
command -v zstd >/dev/null || { sudo apt-get update; sudo apt-get install -y zstd; }

printf 'file\tbytes\trecords\tsource\tdownload_s\n' > "$out/reads_check.tsv"
for m in 1 2; do
  if [ "$kind" = prefix ]; then asset="${run}_${m}.first1000000.fastq.gz"; else asset="${run}_${m}.fastq.gz"; fi
  want=$(awk -F'\t' -v f="$asset" '$1==f{print $3}' "$table")
  records=$(awk -F'\t' -v f="$asset" '$1==f{print $4}' "$table")
  test -n "$want"
  start=$(date +%s.%N)
  curl -fsSL --retry 3 --retry-delay 10 -o "$out/${run}_${m}.fastq.gz" "$rel/$asset"
  secs=$(awk -v a="$start" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f", b-a}')
  got=$(sha256sum "$out/${run}_${m}.fastq.gz" | cut -d' ' -f1)
  echo "$asset table=$want got=$got seconds=$secs"
  test "$want" = "$got"
  printf '%s\t%s\t%s\t%s\t%s\n' "${run}_${m}.fastq.gz" "$(stat -c %s "$out/${run}_${m}.fastq.gz")" "$records" \
    "release page-v1-test-data-v1 asset $asset ($kind), sha256 checked against $table" "$secs" >> "$out/reads_check.tsv"
done

curl -fsSL --retry 3 -o "$out/yeast.transcripts.fa.zst" \
  https://github.com/stajichlab/genomics-wasm/releases/download/qrnamap-spike0-fixtures-v1/yeast.transcripts.fa.zst
tw=$(awk -F'\t' '$1=="transcripts/yeast.transcripts.fa.zst"{print $3}' validation/browser-spike0/manifest.tsv)
tg=$(sha256sum "$out/yeast.transcripts.fa.zst" | cut -d' ' -f1)
echo "yeast.transcripts.fa.zst manifest=$tw got=$tg"
test -n "$tw" && test "$tw" = "$tg"
zstd -q -d -f "$out/yeast.transcripts.fa.zst" -o "$out/yeast.transcripts.fa"
rm "$out/yeast.transcripts.fa.zst"
cat "$out/reads_check.tsv"
(cd "$out" && sha256sum ./*.gz ./*.fa | tee sha256.txt)
