#!/usr/bin/env python3
"""Browser spike 0, Task 3: compare WASM outputs with the native references.

Usage:
  compare_parity.py --native FIXTURES/native --parity PARITY_DIR --out results.tsv

PARITY_DIR holds the downloaded workflow artifacts:
  spike0-parity-<engine>/<tool>/...        engine = node, chromium, firefox, webkit
  spike0-parity-deseq2-node/deseq2/...

Checks (plan Task 3, step 3):
  kallisto       abundance.tsv: same target_id order; est_counts and tpm text identical.
                 run_info.json is not compared (start time, paths).
  featureCounts  count table identical after the '#' program/command line; summary identical.
  HISAT2         SAM records identical after removing @PG header lines.
  DESeq2         same genes, same columns, same NA pattern; every numeric value within
                 DESEQ2_TOL (stated before the run); max absolute and relative differences reported.
Plain, .gz and .zst inputs are read transparently.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import math
import subprocess
import sys
from pathlib import Path

ENGINES = ["node", "chromium", "firefox", "webkit"]
TOOLS = ["kallisto", "hisat2", "featurecounts", "deseq2"]
# Stated before the first parity run. Native: DESeq2 1.50.2 with R 4.5.3 (pixi).
# WASM: release deseq2-package.zip (DESeq2 1.52.0, webR, R 4.6.0). The versions differ,
# so equality is not expected. A value passes when abs diff <= ABS or rel diff <= REL.
DESEQ2_TOL = {"absolute": 1e-12, "relative": 1e-6}
DESEQ2_NUMERIC = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def open_text(path: Path) -> io.TextIOBase:
    """Open plain, .gz or .zst text for reading."""
    path = Path(path)
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    if path.suffix == ".zst":
        proc = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE)
        return io.TextIOWrapper(proc.stdout, encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def find(base: Path, *names: str) -> Path | None:
    for name in names:
        for cand in (base / name, base / (name + ".zst"), base / (name + ".gz")):
            if cand.exists():
                return cand
    return None


def sha256(path: Path | None) -> str:
    if path is None:
        return ""
    h = hashlib.sha256()
    with open_text(path) as fh:
        for line in fh:
            h.update(line.encode("utf-8"))
    return h.hexdigest()


def short(text: str, n: int = 160) -> str:
    text = text.rstrip("\r\n").replace("\t", "\\t")
    return text if len(text) <= n else text[:n] + "..."


def compare_lines(native_iter, wasm_iter, label: str) -> tuple[bool, str, int]:
    """Compare two line streams in order. Return (equal, first difference, lines compared)."""
    n = 0
    for n_line, w_line in zip_longest_lines(native_iter, wasm_iter):
        n += 1
        if n_line != w_line:
            if n_line is None:
                return False, f"{label}: WASM has extra line {n}: {short(w_line)}", n
            if w_line is None:
                return False, f"{label}: WASM ends early; native line {n}: {short(n_line)}", n
            return False, f"{label} line {n}: native '{short(n_line)}' vs WASM '{short(w_line)}'", n
    return True, "", n


def zip_longest_lines(a, b):
    a, b = iter(a), iter(b)
    while True:
        x, y = next(a, None), next(b, None)
        if x is None and y is None:
            return
        yield x, y


def compare_kallisto(native: Path, wasm: Path) -> tuple[bool, str]:
    def rows(p):
        with open_text(p) as fh:
            return list(csv.reader(fh, delimiter="\t"))
    a, b = rows(native), rows(wasm)
    if not a or not b:
        return False, "empty abundance.tsv"
    if a[0] != b[0]:
        return False, f"header differs: native {a[0]} vs WASM {b[0]}"
    hdr = a[0]
    try:
        it, ic, itpm = hdr.index("target_id"), hdr.index("est_counts"), hdr.index("tpm")
    except ValueError as e:
        return False, f"missing column: {e}"
    if len(a) != len(b):
        return False, f"row count differs: native {len(a) - 1} vs WASM {len(b) - 1}"
    for i, (x, y) in enumerate(zip(a[1:], b[1:]), start=2):
        if x[it] != y[it]:
            return False, f"line {i}: target_id native {x[it]} vs WASM {y[it]}"
        for col, j in (("est_counts", ic), ("tpm", itpm)):
            if x[j] != y[j]:
                return False, f"line {i} {x[it]}: {col} native {x[j]} vs WASM {y[j]}"
    same_file = [r for r in a] == [r for r in b]
    return True, f"{len(a) - 1} targets; est_counts and tpm identical; whole file identical: {same_file}"


def compare_featurecounts(native_counts: Path, wasm_counts: Path,
                          native_summary: Path | None, wasm_summary: Path | None) -> tuple[bool, str]:
    def body(p):
        with open_text(p) as fh:
            for line in fh:
                if not line.startswith("#"):
                    yield line
    ok, diff, n = compare_lines(body(native_counts), body(wasm_counts), "counts")
    if not ok:
        return False, diff
    if native_summary is None or wasm_summary is None:
        return False, "summary file missing"
    with open_text(native_summary) as fa, open_text(wasm_summary) as fb:
        ok, diff, m = compare_lines(fa, fb, "summary")
    if not ok:
        return False, diff
    return True, f"{n - 1} genes identical (comment line excluded); summary identical ({m} lines)"


def compare_sam(native: Path, wasm: Path) -> tuple[bool, str]:
    def records(p, counter):
        with open_text(p) as fh:
            for line in fh:
                if line.startswith("@PG"):
                    continue
                if not line.startswith("@"):
                    counter[0] += 1
                yield line
    cn, cw = [0], [0]
    ok, diff, n = compare_lines(records(native, cn), records(wasm, cw), "SAM (no @PG)")
    if not ok:
        return False, diff
    return True, f"{cn[0]} alignment records and header (without @PG) identical"


def _num(s: str) -> float | None:
    if s in ("NA", ""):
        return None
    return float(s)


def compare_deseq2(native: Path, wasm: Path, tol: dict = DESEQ2_TOL) -> tuple[bool, str]:
    def rows(p):
        with open_text(p) as fh:
            return list(csv.DictReader(fh))
    a, b = rows(native), rows(wasm)
    if len(a) != len(b):
        return False, f"row count differs: native {len(a)} vs WASM {len(b)}"
    if not a:
        return False, "empty results"
    if list(a[0].keys()) != list(b[0].keys()):
        return False, f"columns differ: native {list(a[0].keys())} vs WASM {list(b[0].keys())}"
    stats = {c: {"max_abs": 0.0, "max_rel": 0.0, "fail": 0} for c in DESEQ2_NUMERIC}
    first = ""
    for i, (x, y) in enumerate(zip(a, b), start=2):
        if x["gene"] != y["gene"]:
            return False, f"line {i}: gene native {x['gene']} vs WASM {y['gene']}"
        for c in DESEQ2_NUMERIC:
            u, v = _num(x[c]), _num(y[c])
            if (u is None) != (v is None):
                return False, f"line {i} {x['gene']} {c}: NA pattern differs (native {x[c]}, WASM {y[c]})"
            if u is None:
                continue
            d = abs(u - v)
            m = max(abs(u), abs(v))
            r = d / m if m > 0 else 0.0
            s = stats[c]
            s["max_abs"] = max(s["max_abs"], d)
            s["max_rel"] = max(s["max_rel"], r)
            if not (d <= tol["absolute"] or r <= tol["relative"]) or math.isnan(d):
                s["fail"] += 1
                if not first:
                    first = f"line {i} {x['gene']} {c}: native {x[c]} vs WASM {y[c]} (abs {d:.3g}, rel {r:.3g})"
    sig = lambda rs: sum(1 for r in rs if _num(r["padj"]) is not None and _num(r["padj"]) < 0.05)
    parts = [f"{c} max_abs={stats[c]['max_abs']:.3g} max_rel={stats[c]['max_rel']:.3g} n_out={stats[c]['fail']}"
             for c in DESEQ2_NUMERIC]
    summary = (f"{len(a)} genes; tol abs<={tol['absolute']:g} or rel<={tol['relative']:g}; "
               + "; ".join(parts) + f"; padj<0.05 native {sig(a)} WASM {sig(b)}")
    bad = sum(s["fail"] for s in stats.values())
    if bad:
        return False, f"{bad} values out of tolerance; first: {first}; {summary}"
    return True, summary


def step_status(engine_dir: Path, tool: str) -> tuple[str, str]:
    """Read <tool>.steps.tsv written by parity.mjs. Return (worst status, first failure detail)."""
    tsv = engine_dir / f"{tool}.steps.tsv"
    if not tsv.exists():
        return "missing", f"no {tsv.name}"
    with open(tsv, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        if r["status"] != "ok":
            return r["status"], f"step {r['step']} {r['status']}: exit {r['exit_code']}; {r['detail'][:300]}"
    return "ok", ""


def evaluate(native: Path, parity: Path, engine: str, tool: str) -> dict:
    row = {"tool": tool, "engine": engine, "result": "", "detail": "", "wasm_sha256": "",
           "native_sha256": "", "log_ref": ""}
    if tool == "deseq2" and engine != "node":
        row.update(result="skipped", detail="DESeq2 runs in Node only in Task 3 (plan step 2)")
        return row
    art = parity / (f"spike0-parity-deseq2-node" if tool == "deseq2" else f"spike0-parity-{engine}")
    tdir = art / tool
    row["log_ref"] = str(art.relative_to(parity) / f"{tool}.steps.tsv")
    if not art.exists():
        row.update(result="error", detail=f"artifact {art.name} not found (job did not upload)")
        return row
    status, why = step_status(art, tool)
    if status != "ok":
        row.update(result="fail" if status == "fail" else "error", detail=why)
        return row
    if tool == "kallisto":
        n, w = find(native / "kallisto", "abundance.tsv"), find(tdir / "kq", "abundance.tsv")
        res = compare_kallisto(n, w) if n and w else (False, "abundance.tsv missing")
    elif tool == "hisat2":
        n, w = find(native / "hisat2", "yeast.sam"), find(tdir, "yeast.sam")
        res = compare_sam(n, w) if n and w else (False, "yeast.sam missing")
    elif tool == "featurecounts":
        n = find(native / "featurecounts", "yeast.counts.tsv")
        w = find(tdir, "fc.tsv")
        ns = find(native / "featurecounts", "yeast.counts.summary")
        ws = find(tdir, "fc.tsv.summary")
        res = compare_featurecounts(n, w, ns, ws) if n and w else (False, "count table missing")
    else:
        n, w = find(native / "deseq2", "amb.results.csv"), find(tdir, "amb.results.csv")
        res = compare_deseq2(n, w) if n and w else (False, "amb.results.csv missing")
    row["native_sha256"], row["wasm_sha256"] = sha256(n), sha256(w)
    row.update(result="pass" if res[0] else "fail", detail=res[1])
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native", type=Path, required=True)
    ap.add_argument("--parity", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    rows = [evaluate(args.native, args.parity, e, t) for t in TOOLS for e in ENGINES]
    cols = ["tool", "engine", "result", "detail", "wasm_sha256", "native_sha256", "log_ref"]
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: str(r[k]).replace("\t", " ").replace("\n", " ") for k in cols})
    for r in rows:
        print(f"{r['tool']:14s} {r['engine']:9s} {r['result']:8s} {r['detail'][:200]}")
    return 0 if all(r["result"] in ("pass", "skipped") for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
