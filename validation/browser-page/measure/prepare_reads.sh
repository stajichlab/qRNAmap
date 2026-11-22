#!/usr/bin/env bash
# Page-v1 measurements, CI only: download SRR3396381 (both mates) from ENA and the yeast transcripts of the spike 0
# fixture release, and check them. Usage: prepare_reads.sh <out dir>
# Checks: MD5 of each mate equals the ENA filereport fastq_md5; records per mate equal the ENA read_count;
# the transcript FASTA (zstd) SHA-256 equals validation/browser-spike0/manifest.tsv. Writes <out dir>/reads_check.tsv.
# Run from the repository root. No BASH_SOURCE is used.
set -euo pipefail
out="${1:?usage: prepare_reads.sh <out dir>}"
run=SRR3396381
mkdir -p "$out"
report="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,read_count,fastq_bytes,fastq_md5,fastq_ftp&format=tsv"
curl -fsS --retry 3 --max-time 120 "$report" -o "$out/ena_${run}.tsv"
cat "$out/ena_${run}.tsv"
read_count=$(awk -F'\t' 'NR==2{print $2}' "$out/ena_${run}.tsv")
IFS=';' read -r -a bytes <<< "$(awk -F'\t' 'NR==2{print $3}' "$out/ena_${run}.tsv")"
IFS=';' read -r -a md5s <<< "$(awk -F'\t' 'NR==2{print $4}' "$out/ena_${run}.tsv")"
IFS=';' read -r -a ftps <<< "$(awk -F'\t' 'NR==2{print $5}' "$out/ena_${run}.tsv")"
printf 'file\tbytes\tena_bytes\tmd5\tena_md5\trecords\tena_read_count\tdownload_s\n' > "$out/reads_check.tsv"
for i in 0 1; do
  f="$out/${run}_$((i + 1)).fastq.gz"
  start=$(date +%s.%N)
  curl -fsS --retry 5 --retry-delay 10 --max-time 5400 -o "$f" "https://${ftps[$i]}"
  end=$(date +%s.%N)
  got_bytes=$(stat -c %s "$f")
  got_md5=$(md5sum "$f" | cut -d' ' -f1)
  records=$(( $(gzip -dc "$f" | wc -l) / 4 ))
  secs=$(awk -v a="$start" -v b="$end" 'BEGIN{printf "%.1f", b-a}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(basename "$f")" "$got_bytes" "${bytes[$i]}" "$got_md5" "${md5s[$i]}" "$records" "$read_count" "$secs" >> "$out/reads_check.tsv"
  test "$got_bytes" = "${bytes[$i]}"
  test "$got_md5" = "${md5s[$i]}"
  test "$records" = "$read_count"
done
curl -fsSL --retry 3 -o "$out/yeast.transcripts.fa.zst" \
  https://github.com/stajichlab/genomics-wasm/releases/download/qrnamap-spike0-fixtures-v1/yeast.transcripts.fa.zst
want=$(awk -F'\t' '$1=="transcripts/yeast.transcripts.fa.zst"{print $3}' validation/browser-spike0/manifest.tsv)
got=$(sha256sum "$out/yeast.transcripts.fa.zst" | cut -d' ' -f1)
echo "yeast.transcripts.fa.zst manifest=$want got=$got"
test -n "$want" && test "$want" = "$got"
zstd -q -d -f "$out/yeast.transcripts.fa.zst" -o "$out/yeast.transcripts.fa"
rm "$out/yeast.transcripts.fa.zst"
cat "$out/reads_check.tsv"
(cd "$out" && sha256sum ./* | tee sha256.txt)
