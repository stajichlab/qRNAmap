"""collect_amb.py: per-sample numbers of the AMB runs (E2) into one small table.

Usage: python -I collect_amb.py --data DATA_DIR --out OUT_TSV --quant-spt K --align-spt K
       [--task-start MODE:TASK=ISO8601 ...]

Reads, for each sample of DATA_DIR/samples.tsv:
  quant: DATA_DIR/run_quant/SAMPLE/run_info.json (n_processed, n_pseudoaligned, p_pseudoaligned)
  align: the line "SAMPLE<TAB>reads_processed=N<TAB>mapped=X%<TAB>assigned_fraction=F" from the
         array task logs in DATA_DIR/run_align/slurm/logs/, and the featureCounts summary
         DATA_DIR/run_align/SAMPLE/featurecounts.summary (Assigned and the sum of all rows).
  both:  the mtime of DATA_DIR/run_MODE/SAMPLE.genecounts.tsv.zst (UTC, ISO 8601). Within one
         array task the samples run in sheet order, so the difference between consecutive mtimes
         is the time of one sample (a derived value, not a timer). For the first sample of a
         task the time is counted from the task's sacct Start (--task-start).
Assigned fraction (align) = featureCounts Assigned / read pairs processed by HISAT2.
"""
import argparse
import csv
import datetime
import glob
import json
import os
import re


def summary(path):
    rows = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            k, v = line.rstrip("\n").split("\t")[:2]
            rows[k] = int(v)
    return rows


def mtime(path):
    t = os.stat(path).st_mtime
    return datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z", t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--quant-spt", type=int, required=True, help="samples per array task, quant")
    ap.add_argument("--align-spt", type=int, required=True, help="samples per array task, align")
    ap.add_argument("--task-start", action="append", default=[],
                    help="MODE:TASK=ISO8601 with UTC offset (sacct Start), e.g. quant:0=2026-11-20T09:53:15-08:00")
    a = ap.parse_args()
    with open(os.path.join(a.data, "samples.tsv")) as fh:
        samples = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
    log_line = re.compile(r"^(\S+)\treads_processed=(\d+)\tmapped=([0-9.]+)%\tassigned_fraction=([0-9.]+)")
    align_log = {}
    for f in sorted(glob.glob(os.path.join(a.data, "run_align", "slurm", "logs", "qrnamap-samples_*.out"))):
        with open(f) as fh:
            for line in fh:
                m = log_line.match(line)
                if m:
                    align_log[m.group(1)] = (int(m.group(2)), float(m.group(3)), float(m.group(4)), os.path.basename(f))
    out = []
    starts = {}
    for item in a.task_start:
        key, iso = item.split("=", 1)
        mode, task = key.split(":")
        starts[(mode, int(task))] = datetime.datetime.fromisoformat(iso).timestamp()
    spt = {"quant": a.quant_spt, "align": a.align_spt}
    prev = {"quant": None, "align": None}
    for i, s in enumerate(samples):
        with open(os.path.join(a.data, "run_quant", s, "run_info.json")) as fh:
            ri = json.load(fh)
        fc = summary(os.path.join(a.data, "run_align", s, "featurecounts.summary"))
        n, mapped, frac, logname = align_log[s]
        row = {
            "sample": s,
            "quant_n_processed": ri["n_processed"],
            "quant_n_pseudoaligned": ri["n_pseudoaligned"],
            "quant_p_pseudoaligned": ri["p_pseudoaligned"],
            "align_pairs_processed": n,
            "align_hisat2_overall_pct": mapped,
            "align_fc_assigned": fc["Assigned"],
            "align_assigned_fraction": round(fc["Assigned"] / n, 4),
            "align_assigned_fraction_cli_log": frac,
            "align_fc_summary_sum": sum(fc.values()),
            "align_fc_unassigned_multimapping": fc.get("Unassigned_MultiMapping", 0),
            "align_fc_unassigned_nofeatures": fc.get("Unassigned_NoFeatures", 0),
            "align_fc_unassigned_ambiguity": fc.get("Unassigned_Ambiguity", 0),
            "align_fc_unassigned_unmapped": fc.get("Unassigned_Unmapped", 0),
            "align_log": logname,
        }
        for mode in ("quant", "align"):
            iso, t = mtime(os.path.join(a.data, f"run_{mode}", f"{s}.genecounts.tsv.zst"))
            row[f"{mode}_genecounts_mtime"] = iso
            task = i // spt[mode]
            row[f"{mode}_task"] = task
            if i % spt[mode] == 0:
                # First sample of a task: from the task's sacct Start (whole seconds).
                prev[mode] = starts.get((mode, task))
            row[f"{mode}_sample_seconds"] = "" if prev[mode] is None else round(t - prev[mode], 1)
            prev[mode] = t
        out.append(row)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    with open(a.out) as fh:
        print(fh.read(), end="")


if __name__ == "__main__":
    main()
