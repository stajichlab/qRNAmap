#!/usr/bin/env python3
"""Cut preview reads for testsets exactly matching the browser's 4-line boundary cut rule.

Usage: pixi run python validation/browser-page/make_preview_reads.py <out_dir>

The <run>_<mate>.fastq.gz files it writes are published as GitHub release testsets-v1 (flat, with checksums.sha256)
and mirrored to s3://stajichlab/qRNAmap/testsets-v1/ with scripts/mirror-indexes-to-s3.sh testsets-v1. New reads
get a new tag (testsets-v2) and a matching TESTSET_READS_URL in site/analyze/testsets.js.
"""
import hashlib
import json
import os
import sys
import urllib.request
import zlib
import gzip

LOCAL_READS_DIR = "/bigdata/gen220/shared/data-examples/rnaseq/kallisto/data"

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def check_fastq(path: str, expected_records: int):
    records = 0
    with open(path, "rb") as f:
        line_no = 0
        seq_len = 0
        for line in f:
            line_no += 1
            kind = (line_no - 1) % 4
            clean_line = line.rstrip(b"\r\n")
            if kind == 0:
                assert clean_line.startswith(b"@"), f"Line {line_no} does not start with @"
                records += 1
            elif kind == 1:
                seq_len = len(clean_line)
            elif kind == 2:
                assert clean_line.startswith(b"+"), f"Line {line_no} does not start with +"
            elif kind == 3:
                assert len(clean_line) == seq_len, f"Line {line_no} qual len {len(clean_line)} != seq len {seq_len}"
    assert records == expected_records, f"Got {records} records, expected {expected_records}"
    assert line_no == expected_records * 4, f"Got {line_no} lines, expected {expected_records * 4}"

def cut_local_gzip(in_gz: str, out_plain: str, target_lines: int):
    with gzip.open(in_gz, "rb") as fin, open(out_plain, "wb") as fout:
        lines = 0
        for line in fin:
            fout.write(line)
            lines += 1
            if lines == target_lines:
                break
        assert lines == target_lines, f"Expected {target_lines} lines, read {lines}"

def cut_remote_stream(url: str, out_plain: str, target_lines: int):
    req = urllib.request.Request(url, headers={"User-Agent": "qRNAmap-preview-fetch/1.0"})
    decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
    lines = 0
    with urllib.request.urlopen(req, timeout=60) as resp, open(out_plain, "wb") as fout:
        buf = b""
        while lines < target_lines:
            chunk = resp.read(65536)
            if not chunk:
                # Flush remaining
                plain = decomp.flush()
                if plain:
                    buf += plain
                break
            plain = decomp.decompress(chunk)
            if not plain:
                continue
            buf += plain
            # Count newlines in buf
            while b"\n" in buf and lines < target_lines:
                nl_idx = buf.index(b"\n")
                line = buf[:nl_idx + 1]
                buf = buf[nl_idx + 1:]
                fout.write(line)
                lines += 1
                if lines == target_lines:
                    break
    assert lines == target_lines, f"Streamed {lines} lines, expected {target_lines}"

def main():
    if len(sys.argv) < 2:
        out_base = "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/page-v1-preview-reads"
    else:
        out_base = sys.argv[1]

    manifest_path = "site/testsets.json"
    with open(manifest_path) as f:
        testsets = json.load(f)

    os.makedirs(out_base, exist_ok=True)
    summary_rows = []

    for ts in testsets:
        ts_id = ts["id"]
        n_pairs = ts["n_pairs"]
        target_lines = n_pairs * 4
        ts_dir = os.path.join(out_base, ts_id)
        os.makedirs(ts_dir, exist_ok=True)

        print(f"=== Testset: {ts_id} (N={n_pairs} pairs, {target_lines} lines) ===")
        for s in ts["samples"]:
            sample = s["sample"]
            run = s["run"]
            for mate, url_key in [(1, "fastq1_url"), (2, "fastq2_url")]:
                plain_name = f"{run}_{mate}.fastq"
                plain_path = os.path.join(ts_dir, plain_name)
                gz_name = f"{run}_{mate}.fastq.gz"
                gz_path = os.path.join(ts_dir, gz_name)

                already_done = False
                if os.path.isfile(plain_path) and os.path.isfile(gz_path) and os.path.getsize(gz_path) > 0:
                    try:
                        check_fastq(plain_path, n_pairs)
                        already_done = True
                    except Exception:
                        already_done = False

                if not already_done:
                    # Check if local file exists in LOCAL_READS_DIR
                    local_gz = os.path.join(LOCAL_READS_DIR, f"{run}_{mate}.fastq.gz")
                    if os.path.isfile(local_gz):
                        print(f"Cutting from local file: {local_gz} -> {plain_name}")
                        cut_local_gzip(local_gz, plain_path, target_lines)
                    else:
                        url = s[url_key]
                        print(f"Streaming from ENA: {url} -> {plain_name}")
                        cut_remote_stream(url, plain_path, target_lines)

                    check_fastq(plain_path, n_pairs)
                    # Also create .fastq.gz with reproducible gzip
                    with open(plain_path, "rb") as fin, gzip.GzipFile(gz_path, "wb", mtime=0) as fout:
                        fout.writelines(fin)
                else:
                    print(f"Already done and verified: {plain_name}")

                plain_size = os.path.getsize(plain_path)
                plain_hash = sha256_file(plain_path)
                gz_size = os.path.getsize(gz_path)
                gz_hash = sha256_file(gz_path)

                summary_rows.append({
                    "set_id": ts_id, "sample": sample, "run": run, "mate": mate,
                    "filename": plain_name, "bytes": plain_size, "records": n_pairs, "sha256": plain_hash
                })
                summary_rows.append({
                    "set_id": ts_id, "sample": sample, "run": run, "mate": mate,
                    "filename": gz_name, "bytes": gz_size, "records": n_pairs, "sha256": gz_hash
                })

            # Check mate 1 and 2 first read stems match
            p1 = os.path.join(ts_dir, f"{run}_1.fastq")
            p2 = os.path.join(ts_dir, f"{run}_2.fastq")
            with open(p1, "rb") as f1, open(p2, "rb") as f2:
                h1 = f1.readline().decode("latin1").split()[0]
                h2 = f2.readline().decode("latin1").split()[0]
                if h1.endswith("/1"): h1 = h1[:-2]
                if h2.endswith("/2"): h2 = h2[:-2]
                assert h1 == h2, f"Stem mismatch between mates for {run}: {h1} vs {h2}"

    tsv_out = "validation/browser-page/expected/reads_sha256.tsv"
    with open(tsv_out, "w") as f:
        f.write("set_id\tsample\trun\tmate\tfilename\tbytes\trecords\tsha256\n")
        for r in summary_rows:
            f.write(f"{r['set_id']}\t{r['sample']}\t{r['run']}\t{r['mate']}\t{r['filename']}\t{r['bytes']}\t{r['records']}\t{r['sha256']}\n")

    print(f"Wrote {len(summary_rows)} rows to {tsv_out}")

if __name__ == "__main__":
    main()
