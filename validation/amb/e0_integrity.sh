#!/usr/bin/env bash
# e0_integrity.sh: E0 job body. gzip -t and record count of every FASTQ of the AMB set.
# Usage (sbatch): sbatch ... --wrap "bash /ABS/validation/amb/e0_integrity.sh SHARED_DIR OUT_TSV"
# Reads SHARED_DIR read-only. Writes one row per file to OUT_TSV:
#   file  gzip_t_exit  lines  records  seconds
# No BASH_SOURCE, no hardcoded /scratch.
set -euo pipefail
SHARED="${1:?shared data dir}"
OUT="${2:?output tsv}"
RUNS="SRR3396385 SRR3396386 SRR3396387 SRR3396381 SRR3396382 SRR3396384"
echo "node: $(hostname)  cpus: ${SLURM_CPUS_ON_NODE:-?}  date: $(date -Is)"
printf 'file\tgzip_t_exit\tlines\trecords\tseconds\n' > "$OUT"
for r in $RUNS; do
    for m in 1 2; do
        f="$SHARED/${r}_$m.fastq.gz"
        t0=$(date +%s)
        set +e
        gzip -t "$f"; gt=$?
        set -e
        lines=$(zcat "$f" | wc -l)
        t1=$(date +%s)
        rec=$(( lines / 4 ))
        printf '%s\t%s\t%s\t%s\t%s\n' "${r}_$m.fastq.gz" "$gt" "$lines" "$rec" "$(( t1 - t0 ))" | tee -a "$OUT"
    done
done
echo "E0_DONE $(date -Is)"
