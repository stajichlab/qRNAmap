#!/usr/bin/env python3
"""Compare two DE result tables (gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj; NA allowed).
  python3 compare_de.py REFERENCE.tsv OTHER.tsv [--alpha 0.05]
Prints, over genes with a reference padj (DESeq2 leaves NA for filtered genes): the call table at padj < alpha,
the Spearman and Pearson correlation of log2FoldChange and of -log10 pvalue, and the median absolute LFC difference.
Standard library only."""
import argparse, csv, math, statistics

def load(p):
    with open(p) as f:
        return {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}

def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None

def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r

def pearson(a, b):
    ma, mb = statistics.fmean(a), statistics.fmean(b)
    num_ = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num_ / den if den else float("nan")

ap = argparse.ArgumentParser()
ap.add_argument("reference"); ap.add_argument("other"); ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args()
ref, oth = load(a.reference), load(a.other)
genes = [g for g in ref if g in oth and num(ref[g]["padj"]) is not None and num(oth[g]["padj"]) is not None]
both = sum(num(ref[g]["padj"]) < a.alpha and num(oth[g]["padj"]) < a.alpha for g in genes)
ref_only = sum(num(ref[g]["padj"]) < a.alpha <= num(oth[g]["padj"]) for g in genes)
oth_only = sum(num(oth[g]["padj"]) < a.alpha <= num(ref[g]["padj"]) for g in genes)
print(f"genes compared: {len(genes)} (reference padj not NA)")
print(f"padj < {a.alpha}: both {both}, reference only {ref_only}, other only {oth_only}")
for col, tf in (("log2FoldChange", lambda v: v), ("pvalue", lambda v: -math.log10(max(v, 1e-300)))):
    pairs = [(tf(num(ref[g][col])), tf(num(oth[g][col]))) for g in genes if num(ref[g][col]) is not None and num(oth[g][col]) is not None]
    x, y = [p[0] for p in pairs], [p[1] for p in pairs]
    name = col if col != "pvalue" else "-log10 pvalue"
    print(f"{name}: Pearson {pearson(x, y):.4f}, Spearman {pearson(ranks(x), ranks(y)):.4f} (n={len(pairs)})")
d = [abs(num(ref[g]["log2FoldChange"]) - num(oth[g]["log2FoldChange"])) for g in genes
     if num(ref[g]["log2FoldChange"]) is not None and num(oth[g]["log2FoldChange"]) is not None]
print(f"|LFC difference|: median {statistics.median(d):.4f}, 95th percentile {sorted(d)[int(0.95 * (len(d) - 1))]:.4f}")
