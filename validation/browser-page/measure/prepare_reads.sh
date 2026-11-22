#!/usr/bin/env bash
# Page-v1 measurements, CI only: SRR3396381 reads (both mates) and the yeast transcripts of the spike 0 fixture release.
# Usage: prepare_reads.sh <out dir>. Run from the repository root. No BASH_SOURCE is used.
# Note: CI no longer runs this script. CI downloads the files made by it from release page-v1-test-data-v1
# (fetch_reads_release.sh). Use this script to re-create the data.
#
# Source of the reads: the NCBI SRA file of SRR3396381 (AWS open data URL from the NCBI locate API), converted with
# fasterq-dump of the NCBI SRA Toolkit 3.2.0, then gzip -6. Reason: in run 44984006776 (2026-11-22) a full download of the
# ENA FASTQ from a GitHub runner failed (curl error 18 after about 21.6 MB in 480 s).
# Checks: MD5 of the .sra file equals the locate API value; records per mate equal the ENA read_count; the sequence and
# quality lines of the records in the first 1 MiB of each ENA FASTQ (Range request) equal the first records here (headers are
# not compared; this check is recorded, and a failed ENA request is recorded as not checked);
# the transcript FASTA (zstd) SHA-256 equals validation/browser-spike0/manifest.tsv. Writes <out dir>/reads_check.tsv.
set -euo pipefail
out="${1:?usage: prepare_reads.sh <out dir>}"
run=SRR3396381
work="${RUNNER_TEMP:-/tmp}/sra-work"
mkdir -p "$out" "$work"
command -v zstd >/dev/null || { sudo apt-get update; sudo apt-get install -y zstd; }

curl -fsS --retry 3 --max-time 120 -o "$out/ena_${run}.tsv" \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp&format=tsv"
cat "$out/ena_${run}.tsv"
read_count=$(awk -F'\t' 'NR==2{print $2}' "$out/ena_${run}.tsv")
IFS=';' read -r -a ftps <<< "$(awk -F'\t' 'NR==2{print $6}' "$out/ena_${run}.tsv")"

curl -fsS --retry 3 --max-time 120 -o "$work/locate.json" "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=${run}&accept-alternate-locations=yes"
read -r sra_url sra_md5 sra_size < <(python3 -I -c '
import json, sys
d = json.load(open(sys.argv[1]))
for f in d["result"][0]["files"]:
    if f.get("type") == "sra" and f.get("name") == sys.argv[2]:
        for loc in f["locations"]:
            if loc.get("link") and loc.get("service") == "s3" and not loc.get("payRequired"):
                print(loc["link"], f["md5"], f["size"]); sys.exit(0)
sys.exit("no free s3 link for the .sra file")' "$work/locate.json" "$run")
echo "sra $sra_url md5 $sra_md5 size $sra_size"
start=$(date +%s.%N)
curl -fsS --retry 5 --retry-delay 10 --max-time 3600 -o "$work/${run}.sra" "$sra_url"
sra_s=$(awk -v a="$start" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f", b-a}')
test "$(md5sum "$work/${run}.sra" | cut -d' ' -f1)" = "$sra_md5"
test "$(stat -c %s "$work/${run}.sra")" = "$sra_size"

curl -fsSL --retry 3 -o "$work/sratoolkit.tar.gz" https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.0/sratoolkit.3.2.0-ubuntu64.tar.gz
sha256sum "$work/sratoolkit.tar.gz" | tee "$out/sratoolkit.sha256"
tar -xzf "$work/sratoolkit.tar.gz" -C "$work"
bin="$(ls -d "$work"/sratoolkit.3.2.0-ubuntu64/bin)"
mkdir -p "$HOME/.ncbi"
printf '/LIBS/GUID = "%s"\n/libs/cloud/report_instance_identity = "false"\n' "$(cat /proc/sys/kernel/random/uuid)" > "$HOME/.ncbi/user-settings.mkfg"
"$bin/fasterq-dump" --version | tee "$out/fasterq-dump.version"
start=$(date +%s.%N)
"$bin/fasterq-dump" --split-files --threads 4 --temp "$work" --outdir "$work" "$work/${run}.sra"
dump_s=$(awk -v a="$start" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f", b-a}')
rm "$work/${run}.sra"

printf 'file\tbytes\trecords\tena_read_count\tsource\tsra_download_s\tfasterq_dump_s\tena_head_records_compared\tena_head_seq_qual_equal\n' > "$out/reads_check.tsv"
for m in 1 2; do
  fq="$work/${run}_${m}.fastq"
  records=$(( $(wc -l < "$fq") / 4 ))
  # ENA head check: complete records in the first 1 MiB of the ENA file.
  compared=0; equal=not_checked
  if curl -fsS --max-time 600 -r 0-1048575 -o "$work/ena_head_${m}.gz" "https://${ftps[$((m - 1))]}"; then
    gzip -dc "$work/ena_head_${m}.gz" 2>/dev/null > "$work/ena_head_${m}.fastq" || true
    compared=$(( $(wc -l < "$work/ena_head_${m}.fastq") / 4 - 1 ))
    if [ "$compared" -gt 0 ]; then
      a=$(head -n $((compared * 4)) "$work/ena_head_${m}.fastq" | awk 'NR%4==2||NR%4==0' | sha256sum)
      b=$(head -n $((compared * 4)) "$fq" | awk 'NR%4==2||NR%4==0' | sha256sum)
      if [ "$a" = "$b" ]; then equal=yes; else equal=no; fi
      head -n 1 "$work/ena_head_${m}.fastq" > "$out/ena_head_${m}.header.txt"; head -n 1 "$fq" > "$out/sra_head_${m}.header.txt"
    fi
  fi
  gzip -6 -n -c "$fq" > "$out/${run}_${m}.fastq.gz"
  rm "$fq"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${run}_${m}.fastq.gz" "$(stat -c %s "$out/${run}_${m}.fastq.gz")" "$records" "$read_count" \
    "NCBI SRA ${run}.sra (md5 $sra_md5), fasterq-dump 3.2.0 --split-files, gzip -6 -n" "$sra_s" "$dump_s" "$compared" "$equal" >> "$out/reads_check.tsv"
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
(cd "$out" && sha256sum ./*.gz ./*.fa | tee sha256.txt)
