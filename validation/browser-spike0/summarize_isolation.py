#!/usr/bin/env python3
"""Merge the per-browser TSV files of isolation.mjs (one download directory per browser) into results/*.tsv.

Usage: summarize_isolation.py RUN_DIR OUT_DIR
RUN_DIR holds spike0-isolation-chromium, spike0-isolation-firefox, spike0-isolation-webkit (gh run download).
The node probe (no browser) is identical in intent for each browser job; it is written once per job with a column 'job'.
"""
import csv
import sys
from pathlib import Path

run, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
engines = ["chromium", "firefox", "webkit"]
for name in ["isolation_pages.tsv", "isolation_tools.tsv", "cors_results.tsv", "cors_node_probe.tsv"]:
    rows, header = [], None
    for engine in engines:
        with open(run / f"spike0-isolation-{engine}" / name, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            for row in reader:
                if name == "cors_node_probe.tsv":
                    row = {"job": engine, **row}
                rows.append(row)
    cols = (["job"] + header) if name == "cors_node_probe.tsv" else header
    with open(out / name, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(name, len(rows), "rows")
