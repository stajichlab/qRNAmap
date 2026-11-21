#!/usr/bin/env python3
"""zero_genes_by_length.py: Salmon zero-count genes compared with kallisto, by transcript length.

Usage (pixi python, all paths absolute):
  python -I zero_genes_by_length.py --sheet SAMPLES.tsv --salmon-dir RUN_SALMON --kallisto-dir RUN_KALLISTO \
      --authors AUTHORS_COUNTS.tsv --out-dir RESULTS_DIR [--reads-checked 10000]

Inputs (read only):
  SAMPLES.tsv                 qrnamap sample sheet (columns sample, fastq1, fastq2)
  RUN_SALMON/tx2gene.tsv.zst  transcript -> gene
  RUN_SALMON/<s>/quant.sf     Salmon transcript table (Length, NumReads)
  RUN_SALMON/counts.tsv.zst, RUN_KALLISTO/counts.tsv.zst   qrnamap gene matrices
  AUTHORS_COUNTS.tsv          authors' counts with our sample names (validation/geo/prepare_geo.py)
Outputs (small TSVs in RESULTS_DIR):
  zero_genes_by_length.tsv    per sample: Salmon-zero/kallisto-nonzero genes, how many have every
                              transcript no longer than the read length, and a length breakdown
  short_transcripts.tsv       per sample: Salmon transcripts shorter than 150 nt and how many have NumReads 0
  zero_genes_examples.tsv     RPL41A (YDL184C), RPL41B (YDL133C-A), YOR302W: counts per sample in Salmon,
                              kallisto and the authors' table, and transcript lengths
Definitions:
  read length: the maximum sequence length in the first --reads-checked records of each mate file.
  gene length class: the length of the gene's LONGEST transcript (from quant.sf Length).
  "all transcripts <= read length" means the longest transcript is <= the read length.
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys

EXAMPLES = [("YDL184C", "RPL41A"), ("YDL133C-A", "RPL41B"), ("YOR302W", "YOR302W")]
BINS = [(0, 101, "<=101"), (102, 150, "102-150"), (151, 300, "151-300"), (301, 10**9, ">300")]


def zstd_lines(path):
    out = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True, text=True).stdout
    return out.splitlines()


def matrix(path):
    lines = zstd_lines(path) if path.endswith(".zst") else open(path).read().splitlines()
    head = lines[0].split("\t")
    m = {}
    for l in lines[1:]:
        f = l.split("\t")
        m[f[0]] = dict(zip(head[1:], (float(x) for x in f[1:])))
    return m


def read_len(path, n):
    op = gzip.open if path.endswith(".gz") else open
    best = 0
    with op(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i // 4 >= n:
                break
            if i % 4 == 1:
                best = max(best, len(line.rstrip("\n")))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--salmon-dir", required=True)
    ap.add_argument("--kallisto-dir", required=True)
    ap.add_argument("--authors", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--reads-checked", type=int, default=10000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    with open(a.sheet) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    t2g = dict(l.split("\t")[:2] for l in zstd_lines(os.path.join(a.salmon_dir, "tx2gene.tsv.zst")))
    sal = matrix(os.path.join(a.salmon_dir, "counts.tsv.zst"))
    kal = matrix(os.path.join(a.kallisto_dir, "counts.tsv.zst"))
    auth = matrix(a.authors)

    zrows, srows, erows = [], [], []
    for r in rows:
        s = r["sample"]
        rl = max(read_len(r["fastq1"], a.reads_checked), read_len(r["fastq2"], a.reads_checked))
        tx = {}  # gene -> list of (length, numreads)
        with open(os.path.join(a.salmon_dir, s, "quant.sf")) as fh:
            for q in csv.DictReader(fh, delimiter="\t"):
                tx.setdefault(t2g[q["Name"]], []).append((int(q["Length"]), float(q["NumReads"])))
        zero = [g for g in sal if sal[g][s] == 0 and kal[g][s] > 0]
        longest = {g: max(L for L, _ in tx[g]) for g in zero}
        short = sum(1 for g in zero if longest[g] <= rl)
        row = {"sample": s, "read_length": rl, "salmon_zero_kallisto_nonzero": len(zero),
               "all_transcripts_le_read_length": short,
               "kallisto_reads_in_these_genes": int(sum(kal[g][s] for g in zero))}
        for lo, hi, lab in BINS:
            row[f"longest_tx_{lab}"] = sum(1 for g in zero if lo <= longest[g] <= hi)
        zrows.append(row)
        alltx = [t for v in tx.values() for t in v]
        lt150 = [t for t in alltx if t[0] < 150]
        lerl = [t for t in alltx if t[0] <= rl]
        srows.append({"sample": s, "transcripts": len(alltx), "tx_lt_150": len(lt150),
                      "tx_lt_150_numreads_0": sum(1 for t in lt150 if t[1] == 0),
                      "tx_le_read_length": len(lerl),
                      "tx_le_read_length_numreads_0": sum(1 for t in lerl if t[1] == 0)})
        for g, name in EXAMPLES:
            erows.append({"sample": s, "gene": g, "name": name,
                          "transcript_lengths": ",".join(str(L) for L, _ in tx[g]),
                          "salmon": int(sal[g][s]), "kallisto": int(kal[g][s]),
                          "authors": int(auth[g][s]) if g in auth else "NA"})

    for name, data in [("zero_genes_by_length.tsv", zrows), ("short_transcripts.tsv", srows),
                       ("zero_genes_examples.tsv", erows)]:
        with open(os.path.join(a.out_dir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
        print("==", name)
        print(open(os.path.join(a.out_dir, name)).read(), end="")


if __name__ == "__main__":
    sys.exit(main())
