#!/usr/bin/env python3
"""merge_counts.py DIR OUT SAMPLE... : merge DIR/<sample>.genecounts.tsv.zst (gene, <sample>) into one table.
Genes are kept in the order of the first sample; every sample must list the same genes."""
import subprocess, sys

def read(path):
    txt = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True, text=True).stdout
    rows = [l.split("\t") for l in txt.splitlines()]
    return rows[0], rows[1:]

def main():
    d, out, samples = sys.argv[1], sys.argv[2], sys.argv[3:]
    genes, cols = None, []
    for s in samples:
        hdr, rows = read(f"{d}/{s}.genecounts.tsv.zst")
        if hdr != ["gene", s]:
            sys.exit(f"error: unexpected header {hdr} for {s}")
        g = [r[0] for r in rows]
        if genes is None:
            genes = g
        elif g != genes:
            sys.exit(f"error: gene list of {s} differs")
        cols.append([r[1] for r in rows])
    with open(out, "w") as f:
        f.write("gene\t" + "\t".join(samples) + "\n")
        for i, g in enumerate(genes):
            f.write(g + "\t" + "\t".join(c[i] for c in cols) + "\n")
    print(f"{len(genes)} genes x {len(samples)} samples -> {out}", file=sys.stderr)

if __name__ == "__main__":
    main()
