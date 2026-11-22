#!/usr/bin/env python3
"""Compare kallisto abundance.tsv files with a reference file, line by line, and check each file on its own.

Usage: check_abundance.py REFERENCE FILE...
Per file: rows, rows that differ from the reference (with the first difference shown as old -> new), whether every number
parses, and the sum of est_counts. Stdlib only. The files are data.
"""
import sys

ref_path, paths = sys.argv[1], sys.argv[2:]
ref = open(ref_path, encoding='utf-8').read().splitlines()
print('file\trows\trows_differing\tall_numbers_parse\test_counts_sum\tfirst_difference')
for p in paths:
    lines = open(p, encoding='utf-8').read().splitlines()
    diffs = [(i, a, b) for i, (a, b) in enumerate(zip(ref, lines)) if a != b]
    diffs += [(i, '', '') for i in range(min(len(ref), len(lines)), max(len(ref), len(lines)))]
    ok, total = True, 0.0
    for line in lines[1:]:
        f = line.split('\t')
        try:
            vals = [float(x) for x in f[1:]]
            if len(vals) != 4:
                ok = False
            else:
                total += vals[2]
        except ValueError:
            ok = False
    first = ''
    if diffs:
        i, a, b = diffs[0]
        first = f'line {i + 1}: {a!r} -> {b!r}'
    print(f'{p}\t{len(lines)}\t{len(diffs)}\t{ok}\t{total:.3f}\t{first}')
