#!/usr/bin/env python3
"""Browser spike 0, Task 4: merge the per-engine index_size.<engine>.tsv files and apply the pass criteria.

Usage: summarize_index_size.py NATIVE_TIMES_TSV OUT_DIR ENGINE_TSV [ENGINE_TSV ...]

Writes OUT_DIR/index_size_runs.tsv (every run, all columns, plus run_id) and OUT_DIR/index_size_summary.tsv
(one row per case and engine: the three values of time and memory, outcome, native reference, criterion).

Criterion (design section 15): the 12 Mb yeast HISAT2 index builds in 600 s or less and under 2 GiB.
"Under 2 GiB" is checked on every memory probe the engine has: WASM linear memory (all engines),
node_maxrss_kb (Node), tree_pss_peak_kb (browsers). A run that does not finish fails the criterion.
The same limits are reported for the other cases as information, not as a criterion.
"""
import csv
import sys
from pathlib import Path

LIMIT_BYTES = 2 * 2**30  # 2 GiB
LIMIT_S = 600.0


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def run_ok(row):
    """Return (ok, reasons) for one run against 600 s and 2 GiB on every probe present."""
    reasons = []
    if row["outcome"] != "finished":
        reasons.append(f"outcome {row['outcome']}")
    run_s = num(row["run_s"])
    if run_s is None or run_s > LIMIT_S:
        reasons.append(f"run_s {row['run_s']}")
    probes = {
        "wasm_heap_final_bytes": num(row["wasm_heap_final_bytes"]),
        "node_maxrss_kb": None if num(row["node_maxrss_kb"]) is None else num(row["node_maxrss_kb"]) * 1024,
        "tree_pss_peak_kb": None if num(row["tree_pss_peak_kb"]) in (None, 0.0) else num(row["tree_pss_peak_kb"]) * 1024,
    }
    if probes["wasm_heap_final_bytes"] is None:
        reasons.append("no WASM heap value")
    for name, value in probes.items():
        if value is not None and value >= LIMIT_BYTES:
            reasons.append(f"{name} >= 2 GiB")
    return (not reasons, reasons)


def main(argv):
    if len(argv) < 4:
        sys.exit(__doc__)
    native = {row["step"]: row for row in read_tsv(argv[1])}
    out = Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    runs = []
    for path in argv[3:]:
        run_id = next((part.rsplit("run-", 1)[1] for part in Path(path).parts if "run-" in part), "")
        for row in read_tsv(path):
            row["run_id"] = run_id
            runs.append(row)
    columns = ["run_id"] + [c for c in runs[0].keys() if c != "run_id"]
    with open(out / "index_size_runs.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(runs)

    native_step = {"hisat2_yeast": "hisat2_build_yeast", "hisat2_afum": "hisat2_build_afum",
                   "kallisto_yeast": "kallisto_index", "kallisto_afum": ""}
    order_case = ["hisat2_yeast", "hisat2_afum", "kallisto_yeast", "kallisto_afum"]
    order_engine = ["node", "chromium", "firefox", "webkit"]
    groups = {}
    for row in runs:
        groups.setdefault((row["case"], row["engine"]), []).append(row)
    header = ["case", "engine", "runs", "outcomes", "run_s", "callmain_s", "wasm_heap_final_MiB", "node_maxrss_MiB",
              "tree_pss_peak_MiB", "tree_pss_baseline_MiB", "cdp_js_heap_used_peak_MiB", "index_bytes",
              "native_wall_s", "native_peak_rss_MiB", "criterion", "criterion_detail"]
    mib = lambda v, scale: "" if num(v) is None else f"{num(v) * scale / 2**20:.1f}"
    rows = []
    for case in order_case:
        for engine in order_engine:
            group = sorted(groups.get((case, engine), []), key=lambda r: int(r["rep"]))
            if not group:
                continue
            nat = native.get(native_step[case]) if native_step[case] else None
            joined = lambda key, scale=None: ",".join(r[key] if scale is None else mib(r[key], scale) for r in group)
            checks = [run_ok(r) for r in group]
            if case == "hisat2_yeast":
                criterion = "pass" if all(ok for ok, _ in checks) else "fail"
            else:
                criterion = "info: within 600 s and 2 GiB" if all(ok for ok, _ in checks) else "info: not within 600 s and 2 GiB"
            detail = "; ".join(f"rep {r['rep']}: {', '.join(reasons)}" for r, (ok, reasons) in zip(group, checks) if not ok)
            rows.append({
                "case": case, "engine": engine, "runs": len(group), "outcomes": joined("outcome"),
                "run_s": joined("run_s"), "callmain_s": joined("callmain_s"),
                "wasm_heap_final_MiB": joined("wasm_heap_final_bytes", 1),
                "node_maxrss_MiB": joined("node_maxrss_kb", 1024) if engine == "node" else "",
                "tree_pss_peak_MiB": joined("tree_pss_peak_kb", 1024) if engine != "node" else "",
                "tree_pss_baseline_MiB": joined("tree_pss_baseline_kb", 1024) if engine != "node" else "",
                "cdp_js_heap_used_peak_MiB": joined("cdp_js_heap_used_peak_bytes", 1) if engine == "chromium" else "",
                "index_bytes": joined("index_bytes"),
                "native_wall_s": nat["wall_s"] if nat else "not measured",
                "native_peak_rss_MiB": mib(nat["peak_rss_kb"], 1024) if nat else "not measured",
                "criterion": criterion, "criterion_detail": detail,
            })
    with open(out / "index_size_summary.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print("\t".join(str(row[h]) for h in header))


if __name__ == "__main__":
    main(sys.argv)
