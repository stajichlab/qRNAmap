#!/usr/bin/env bash
# Reads for quant-memory.mjs (CI). Usage: quant-memory-reads.sh <out dir> <yeast|h99>. Run from the repository root. No BASH_SOURCE.
#   yeast: the first 1,000,000 pairs of SRR3396381 (2 x 101 bases) from release page-v1-test-data-v1, SHA-256 checked against
#          validation/browser-page/measure/reads_prefix.tsv. Writes yeast_{1,2}.fastq.gz (as downloaded) and yeast_{1,2}.fastq
#          (the same records, plain). The page cuts both at N = 200,000 pairs itself.
#   h99:   the first 200,000 pairs of C. neoformans H99 run SRR23037513 (BioProject PRJNA922535) from NCBI SRA with fastq-dump -X
#          (SRA Toolkit 3.2.0, as prepare_reads.sh). Checks: 200,000 records per mate. The read lengths are written, not assumed.
# Appends one row per file to <out dir>/reads_check.tsv (file, bytes, records, read lengths, source).
set -euo pipefail
out="${1:?usage: quant-memory-reads.sh <out dir> <yeast|h99>}"
kind="${2:?usage: quant-memory-reads.sh <out dir> <yeast|h99>}"
mkdir -p "$out"
check="$out/reads_check.tsv"
[ -f "$check" ] || printf 'file\tbytes\trecords\tread_lengths\tsource\n' > "$check"
stats() { # <file> <source>: records and read length counts of a plain or gzip FASTQ
  local f="$1"
  local s
  s=$( (if [[ $f == *.gz ]]; then gzip -dc "$f"; else cat "$f"; fi) | awk 'NR%4==2{n++; l[length($0)]++} END{printf "%d\t", n; sep=""; for(k in l){printf "%s%s:%d", sep, k, l[k]; sep=","}}')
  printf '%s\t%s\t%s\t%s\n' "$(basename "$f")" "$(stat -c %s "$f")" "$s" "$2" >> "$check"
  echo "$(basename "$f") $s"
}
case "$kind" in
  yeast)
    rel=https://github.com/stajichlab/qRNAmap/releases/download/page-v1-test-data-v1
    table=validation/browser-page/measure/reads_prefix.tsv
    for m in 1 2; do
      asset="SRR3396381_${m}.first1000000.fastq.gz"
      want=$(awk -F'\t' -v f="$asset" '$1==f{print $3}' "$table")
      test -n "$want"
      curl -fsSL --retry 3 --retry-delay 10 -o "$out/yeast_${m}.fastq.gz" "$rel/$asset"
      got=$(sha256sum "$out/yeast_${m}.fastq.gz" | cut -d' ' -f1)
      echo "$asset table=$want got=$got"
      test "$want" = "$got"
      gzip -dc "$out/yeast_${m}.fastq.gz" > "$out/yeast_${m}.fastq"
      stats "$out/yeast_${m}.fastq.gz" "release page-v1-test-data-v1 $asset, sha256 checked against $table"
      stats "$out/yeast_${m}.fastq" "gzip -dc of yeast_${m}.fastq.gz"
    done
    ;;
  h99)
    run=SRR23037513
    work="${RUNNER_TEMP:-/tmp}/sra-work"
    mkdir -p "$work"
    curl -fsSL --retry 3 -o "$work/sratoolkit.tar.gz" https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.0/sratoolkit.3.2.0-ubuntu64.tar.gz
    sha256sum "$work/sratoolkit.tar.gz"
    tar -xzf "$work/sratoolkit.tar.gz" -C "$work"
    bin="$(ls -d "$work"/sratoolkit.3.2.0-ubuntu64/bin)"
    mkdir -p "$HOME/.ncbi"
    printf '/LIBS/GUID = "%s"\n/libs/cloud/report_instance_identity = "false"\n' "$(cat /proc/sys/kernel/random/uuid)" > "$HOME/.ncbi/user-settings.mkfg"
    "$bin/fastq-dump" --version
    "$bin/fastq-dump" -X 200000 --split-files --outdir "$work" "$run"
    for m in 1 2; do
      mv "$work/${run}_${m}.fastq" "$out/h99_${m}.fastq"
      stats "$out/h99_${m}.fastq" "NCBI SRA $run, fastq-dump 3.2.0 -X 200000 --split-files"
      test "$(awk -F'\t' -v f="h99_${m}.fastq" '$1==f{print $3}' "$check")" = 200000
    done
    ;;
  *) echo "kind must be yeast or h99" >&2; exit 2 ;;
esac
cat "$check"
