#!/usr/bin/env bash
# Browser spike 0, Task 2: native reference outputs. Run under SLURM. No BASH_SOURCE.
# Required environment (absolute paths): FIX, OUT (reference dir), PIXI_BIN, SCRIPTS.
# Every tool runs with one thread and under /usr/bin/time -v. Temp files go to node-local $SCRATCH.
set -euo pipefail
: "${FIX:?}" "${OUT:?}" "${PIXI_BIN:?}" "${SCRIPTS:?}"
export PATH="$PIXI_BIN:$PATH"
W="${SCRATCH:?}/s0native"; mkdir -p "$W" "$OUT"
cd "$W"
TIMES="$OUT/native_times.tsv"
printf 'step\tcommand\texit\twall_s\tpeak_rss_kb\n' > "$TIMES"
timed() { # timed STEP COMMAND... : run under /usr/bin/time -v, append one row to native_times.tsv
    local step=$1; shift
    local tf="$W/time_$step.txt"
    set +e
    /usr/bin/time -v -o "$tf" "$@" > "$W/$step.stdout" 2> "$W/$step.stderr"
    local rc=$?
    set -e
    local wall rss
    wall=$(awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; print s}' "$tf")
    rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tf")
    printf '%s\t%s\t%s\t%s\t%s\n' "$step" "$*" "$rc" "$wall" "$rss" >> "$TIMES"
    cp "$W/$step.stderr" "$OUT/$step.stderr.txt"; cp "$tf" "$OUT/$step.time.txt"
    [ "$rc" -eq 0 ] || { echo "step $step failed ($rc)"; exit "$rc"; }
}
{ echo "kallisto $(kallisto version)"; hisat2-align-s --version | head -1; featureCounts -v 2>&1 | head -2
  Rscript -e 'cat("DESeq2", as.character(packageVersion("DESeq2")), R.version.string, "\n")'; zstd -V; } > "$OUT/versions.txt" 2>&1

R1="$FIX/reads/yeast_sub_1.fastq.gz"; R2="$FIX/reads/yeast_sub_2.fastq.gz"
# kallisto: index and quant against the yeast transcripts
zstd -dc "$FIX/transcripts/yeast.transcripts.fa.zst" > yeast.transcripts.fa
timed kallisto_index kallisto index -i yeast.idx yeast.transcripts.fa
timed kallisto_quant kallisto quant -i yeast.idx -o kq -t 1 --plaintext "$R1" "$R2"
mkdir -p "$OUT/kallisto"
for f in abundance.tsv run_info.json; do cp "kq/$f" "$OUT/kallisto/$f"; done
zstd -q -19 -f -o "$OUT/kallisto/yeast.idx.zst" yeast.idx
# HISAT2: both genomes, same reads (the Af293 run is cross-species, to measure build cost)
for g in yeast afum; do
    fa=$(ls "$FIX/genomes/$g"/*_genomic.fna.gz)
    zcat "$fa" > "$g.fa"
    timed "hisat2_build_$g" hisat2-build-s -p 1 "$g.fa" "$g.hs"
    timed "hisat2_align_$g" hisat2-align-s -p 1 -x "$g.hs" -1 <(zcat "$R1") -2 <(zcat "$R2") -S "$g.sam"
    mkdir -p "$OUT/hisat2"
    zstd -q -19 -f -o "$OUT/hisat2/$g.sam.zst" "$g.sam"
    cp "$W/hisat2_align_$g.stderr" "$OUT/hisat2/$g.align_summary.txt"
    cat "$g".hs.*.ht2 | wc -c > "$OUT/hisat2/$g.index_bytes.txt"
done
# featureCounts on the yeast HISAT2 SAM: SAF from qrnamap transcripts, pairs, unstranded (as the CLI does)
mkdir -p "$OUT/featurecounts"
timed featurecounts featureCounts -T 1 -F SAF -a "$FIX/transcripts/yeast.genes.saf" -p --countReadPairs -s 0 -o fc.tsv yeast.sam
zstd -q -19 -f -o "$OUT/featurecounts/yeast.counts.tsv.zst" fc.tsv
cp fc.tsv.summary "$OUT/featurecounts/yeast.counts.summary"
# DESeq2 on the AMB gene-count matrix
mkdir -p "$OUT/deseq2"
timed deseq2 Rscript "$SCRIPTS/deseq2_ref.R" "$FIX/counts/amb_gene_counts.tsv" "$OUT/deseq2/amb"
echo "done $(date -Is)"
