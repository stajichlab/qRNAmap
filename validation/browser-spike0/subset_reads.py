#!/usr/bin/env python3
"""Deterministic paired FASTQ subset. Usage: subset_reads.py R1 R2 TOTAL_PAIRS N SEED OUT1 OUT2
Picks N record indexes with random.Random(SEED).sample(range(TOTAL_PAIRS), N), writes them in file order.
Fails if a file has a different record count than TOTAL_PAIRS."""
import gzip, random, sys

def main():
    r1, r2, total, n, seed, o1, o2 = sys.argv[1:8]
    total, n, seed = int(total), int(n), int(seed)
    keep = sorted(random.Random(seed).sample(range(total), n))
    for src, dst in ((r1, o1), (r2, o2)):
        # gzip with mtime=0 so the output bytes are reproducible
        with gzip.open(src, "rt") as f, open(dst, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as out:
            i = 0
            k = 0
            while True:
                rec = [f.readline() for _ in range(4)]
                if not rec[0]:
                    break
                if k < n and i == keep[k]:
                    out.write("".join(rec).encode())
                    k += 1
                i += 1
        if i != total or k != n:
            sys.exit(f"error: {src} has {i} records (expected {total}), kept {k}/{n}")
        print(f"{src}: {i} records read, {k} kept -> {dst}", file=sys.stderr)

if __name__ == "__main__":
    main()
