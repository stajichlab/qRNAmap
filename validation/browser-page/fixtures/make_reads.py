#!/usr/bin/env python3
"""Make the small read fixtures for the browser tests of the samples step, and the expected CLI gene counts.

Usage: python3 make_reads.py <fixtures dir> <qrnamap binary>

The reads are simulated from the transcripts that `qrnamap transcripts` makes from mito.fna and mito.gff3 (the yeast
mitochondrial genome of make_fixtures.py). Real yeast reads were not used: almost all of them come from nuclear genes,
so they would give almost no counts against this genome.
- Random generator: Python `random.Random(SEED + k)` for sample k. The output is the same on every run.
- Transcripts of at least 150 bases are drawn with a weight per sample (log-normal, so the samples differ).
  A fragment of N(200, 20) bases (cut to 120 bases or more, and to the transcript length) is taken at a random place.
  Paired: mate 1 = the first 50 bases of the fragment, mate 2 = the reverse complement of the last 50 bases. Half of the
  pairs are swapped (an unstranded library). Single-end: one 50-base read from either strand.
  One read of 10 is random sequence (it does not pseudoalign). Quality is `I` for every base. No sequencing errors.
- s1: 300 pairs, gzip (s1_R1.fastq.gz, s1_R2.fastq.gz). s2: 300 pairs, plain (s2_R1.fastq, s2_R2.fastq).
  s3: 300 single-end reads, gzip (s3.fastq.gz). gzip: one member, mtime 0, level 9.
- expected_quant.tsv: `case<TAB>gene<TAB>count`, the gene counts (`<sample>.genecounts.tsv.zst` without its header) of `qrnamap quant` (kallisto 0.52.0 from PATH) for
  case s1, s2, s3 (all reads) and s1_first100 (the first 100 pairs of s1, the page cut with N = 100).
  s3 runs with --fragment-length 200,20.
- expected_runinfo.tsv: `case<TAB>pairs<TAB>n_processed<TAB>n_pseudoaligned` from run_info.json of the same runs.
- PROVENANCE.tsv: the rows of the files above are replaced (other rows are kept).
Python standard library, plus `zstd` on PATH to read the CLI output. No network.
"""
import gzip
import hashlib
import json
import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path

SEED = 20261122
READ_LEN = 50
PAIRS = 300
COMP = str.maketrans("ACGTN", "TGCAN")
OURS = ("s1_R1.fastq.gz", "s1_R2.fastq.gz", "s2_R1.fastq", "s2_R2.fastq", "s3.fastq.gz",
        "expected_quant.tsv", "expected_runinfo.tsv", "make_reads.py")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(text: str) -> dict:
    seqs, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def simulate(tx: dict, k: int, paired: bool) -> tuple:
    rng = random.Random(SEED + k)
    names = sorted(n for n, s in tx.items() if len(s) >= 150)
    weights = [math.exp(rng.gauss(0, 1.5)) for _ in names]
    r1, r2 = [], []
    for i in range(PAIRS):
        rid = f"s{k}_{i + 1:06d}"
        if rng.random() < 0.1:
            a = "".join(rng.choice("ACGT") for _ in range(READ_LEN))
            b = "".join(rng.choice("ACGT") for _ in range(READ_LEN))
        else:
            seq = tx[rng.choices(names, weights)[0]]
            flen = max(120, min(len(seq), round(rng.gauss(200, 20))))
            start = rng.randrange(0, len(seq) - flen + 1)
            frag = seq[start:start + flen]
            a, b = frag[:READ_LEN], revcomp(frag[-READ_LEN:])
            if rng.random() < 0.5:
                a, b = b, a
        q = "I" * READ_LEN
        r1.append(f"@{rid}/1\n{a}\n+\n{q}\n")
        r2.append(f"@{rid}/2\n{b}\n+\n{q}\n")
    return "".join(r1).encode(), ("".join(r2).encode() if paired else None)


def gz(data: bytes) -> bytes:
    return gzip.compress(data, compresslevel=9, mtime=0)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def cli_quant(qrnamap: str, tmp: Path, case: str, fastq1: Path, fastq2, tx_fa: Path, tx2gene: Path):
    work = tmp / case
    work.mkdir()
    sheet = work / "samples.tsv"
    sheet.write_text(f"sample\tfastq1\tfastq2\n{case}\t{fastq1}\t{fastq2 or ''}\n")
    cmd = [qrnamap, "quant", "--samples", str(sheet), "--transcripts", str(tx_fa), "--tx2gene", str(tx2gene),
           "--out", str(work / "out"), "--tmpdir", str(work / "tmp")]
    if fastq2 is None:
        cmd += ["--fragment-length", "200,20"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    counts = subprocess.run(["zstd", "-dc", str(work / "out" / f"{case}.genecounts.tsv.zst")],
                            check=True, capture_output=True, text=True).stdout
    info = json.loads((work / "out" / case / "run_info.json").read_text())
    return counts, info


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    out, qrnamap = Path(sys.argv[1]), sys.argv[2]
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        tx_fa, tx2gene = tmp / "tx.fa", tmp / "tx2gene.tsv"
        subprocess.run([qrnamap, "transcripts", "--reference", str(out / "mito.fna"), "--gff", str(out / "mito.gff3"),
                        "--out", str(tx_fa), "--tx2gene", str(tx2gene)], check=True)
        tx = read_fasta(tx_fa.read_text())

        s1 = simulate(tx, 1, True)
        s2 = simulate(tx, 2, True)
        s3 = simulate(tx, 3, False)
        files = {
            "s1_R1.fastq.gz": gz(s1[0]), "s1_R2.fastq.gz": gz(s1[1]),
            "s2_R1.fastq": s2[0], "s2_R2.fastq": s2[1],
            "s3.fastq.gz": gz(s3[0]),
        }
        for name, data in files.items():
            (out / name).write_bytes(data)

        first = lambda data, n: b"".join(data.splitlines(keepends=True)[:4 * n])
        plain = tmp / "plain"
        plain.mkdir()
        cases = {"s1": (s1[0], s1[1]), "s2": (s2[0], s2[1]), "s3": (s3[0], None),
                 "s1_first100": (first(s1[0], 100), first(s1[1], 100))}
        quant_rows, info_rows = ["case\tgene\tcount"], ["case\tpairs\tn_processed\tn_pseudoaligned"]
        for case, (m1, m2) in cases.items():
            p1 = plain / f"{case}_1.fastq"
            p1.write_bytes(m1)
            p2 = None
            if m2 is not None:
                p2 = plain / f"{case}_2.fastq"
                p2.write_bytes(m2)
            counts, info = cli_quant(qrnamap, tmp, case, p1, p2, tx_fa, tx2gene)
            for line in counts.splitlines()[1:]:  # the first line is the header `gene<TAB><sample>`
                if line.strip():
                    quant_rows.append(f"{case}\t{line}")
            pairs = m1.count(b"\n") // 4
            info_rows.append(f"{case}\t{pairs}\t{info['n_processed']}\t{info['n_pseudoaligned']}")
        (out / "expected_quant.tsv").write_text("\n".join(quant_rows) + "\n")
        (out / "expected_runinfo.tsv").write_text("\n".join(info_rows) + "\n")

    version = subprocess.run([qrnamap, "--version"], capture_output=True, text=True).stdout.strip()
    prov = out / "PROVENANCE.tsv"
    kept = [l for l in prov.read_text().splitlines() if l.split("\t")[0] not in OURS]
    for name in OURS:
        data = (out / name).read_bytes() if name != "make_reads.py" else Path(__file__).read_bytes()
        what = {"expected_quant.tsv": f"made by make_reads.py ({version} quant, kallisto 0.52.0)",
                "expected_runinfo.tsv": f"made by make_reads.py ({version} quant, kallisto 0.52.0)",
                "make_reads.py": "the script that made the read fixtures"}.get(name, f"made by make_reads.py (simulated, seed {SEED})")
        kept.append(f"{name}\t{len(data)}\t{sha256(data)}\t{what}")
    prov.write_text("\n".join(kept) + "\n")
    for name in OURS[:-1]:
        print(name, (out / name).stat().st_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
