#!/usr/bin/env python3
"""collect_salmon.py: per-sample table of the Salmon AMB run.

Usage (pixi python): python -I collect_salmon.py --run-dir RUN --sheet SAMPLES.tsv \
    --task-start 2026-11-20T17:13:47 --out per_sample.tsv
Columns: sample, salmon_version, num_processed, num_mapped, percent_mapped (all from
RUN/<s>/aux_info/meta_info.json), genecounts_mtime, seconds (difference of the mtimes of
RUN/<s>.genecounts.tsv.zst to the previous sample in sheet order, or to --task-start, the sacct Start
of the array task, for the first sample; local time of this machine), s_per_million_pairs.
"""
import argparse
import csv
import datetime
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--task-start", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.sheet) as fh:
        samples = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
    prev = datetime.datetime.fromisoformat(a.task_start).timestamp()
    rows = []
    for s in samples:
        with open(os.path.join(a.run_dir, s, "aux_info", "meta_info.json")) as fh:
            m = json.load(fh)
        t = os.stat(os.path.join(a.run_dir, f"{s}.genecounts.tsv.zst")).st_mtime
        sec = t - prev
        prev = t
        rows.append({"sample": s, "salmon_version": m["salmon_version"],
                     "num_processed": m["num_processed"], "num_mapped": m["num_mapped"],
                     "percent_mapped": m["percent_mapped"],
                     "genecounts_mtime": datetime.datetime.fromtimestamp(t).isoformat(timespec="seconds"),
                     "seconds": round(sec), "s_per_million_pairs": round(sec / m["num_processed"] * 1e6, 1)})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(open(a.out).read(), end="")


if __name__ == "__main__":
    main()
