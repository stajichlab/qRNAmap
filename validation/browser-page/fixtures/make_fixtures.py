#!/usr/bin/env python3
"""Make the small genome fixtures for the browser tests of the genome step.

Usage: python3 make_fixtures.py <yeast .fna.gz> <yeast .gff.gz> <output dir> [qrnamap binary]

Source: the S. cerevisiae S288C RefSeq files (GCF_000146045.2, R64) that `qrnamap fetch --source ncbi` downloaded
(copy in qRNAmap_data/spike0/fixtures/genomes/yeast/). The fixtures hold only the mitochondrial genome (NC_001224.1):
- mito.fna     the NC_001224.1 record of the FASTA, unchanged.
- mito.gff3    `##gff-version 3`, the `##sequence-region NC_001224.1` line, and every feature line of NC_001224.1, unchanged.
- mito.fna.gz, mito.gff3.gz   the same bytes, gzip (one member, mtime 0, level 9), so the output is the same on every run.
- PROVENANCE.tsv   SHA-256 and size of the sources and of each output.
- expected.tsv     genes and transcripts that `qrnamap transcripts` makes from mito.fna and mito.gff3 (if a binary is given).
Python standard library only. No network.
"""
import gzip
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

SEQ = "NC_001224.1"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    if len(sys.argv) not in (4, 5):
        print(__doc__, file=sys.stderr)
        return 2
    fna_src, gff_src, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    qrnamap = sys.argv[4] if len(sys.argv) == 5 else None
    out.mkdir(parents=True, exist_ok=True)
    fna_gz, gff_gz = fna_src.read_bytes(), gff_src.read_bytes()

    lines, keep = [], False
    for line in gzip.decompress(fna_gz).decode("ascii").splitlines(keepends=True):
        if line.startswith(">"):
            keep = line[1:].split()[0] == SEQ
        if keep:
            lines.append(line)
    fna = "".join(lines).encode("ascii")
    if not fna:
        raise SystemExit(f"{SEQ} not found in {fna_src}")

    gff_lines = ["##gff-version 3\n"]
    for line in gzip.decompress(gff_gz).decode("utf-8").splitlines(keepends=True):
        if line.startswith(f"##sequence-region {SEQ} ") or line.split("\t", 1)[0] == SEQ:
            gff_lines.append(line)
    gff = "".join(gff_lines).encode("utf-8")

    outputs = {
        "mito.fna": fna,
        "mito.gff3": gff,
        "mito.fna.gz": gzip.compress(fna, compresslevel=9, mtime=0),
        "mito.gff3.gz": gzip.compress(gff, compresslevel=9, mtime=0),
    }
    for name, data in outputs.items():
        (out / name).write_bytes(data)

    rows = ["file\tbytes\tsha256\twhat"]
    rows.append(f"{fna_src.name}\t{len(fna_gz)}\t{sha256(fna_gz)}\tsource FASTA (GCF_000146045.2, NCBI catalog URL)")
    rows.append(f"{gff_src.name}\t{len(gff_gz)}\t{sha256(gff_gz)}\tsource GFF3 (GCF_000146045.2, NCBI catalog URL)")
    for name, data in outputs.items():
        rows.append(f"{name}\t{len(data)}\t{sha256(data)}\tmade by make_fixtures.py ({SEQ} only)")
    (out / "PROVENANCE.tsv").write_text("\n".join(rows) + "\n")

    if qrnamap:
        with tempfile.TemporaryDirectory() as tmp:
            t = Path(tmp)
            subprocess.run([qrnamap, "transcripts", "--reference", str(out / "mito.fna"), "--gff", str(out / "mito.gff3"),
                            "--out", str(t / "tx.fa"), "--tx2gene", str(t / "tx2gene.tsv")], check=True)
            tx = (t / "tx.fa").read_bytes()
            n_tx = tx.count(b">")
            genes = {line.split("\t")[1] for line in (t / "tx2gene.tsv").read_text().splitlines() if line}
            version = subprocess.run([qrnamap, "--version"], capture_output=True, text=True).stdout.strip()
        (out / "expected.tsv").write_text(
            "fixture\tgenes\ttranscripts\ttranscripts_sha256\thow\n"
            f"mito\t{len(genes)}\t{n_tx}\t{sha256(tx)}\t{version} transcripts --reference mito.fna --gff mito.gff3 (genes: distinct gene names in tx2gene)\n")
    for name in sorted(p.name for p in out.iterdir()):
        print(name, (out / name).stat().st_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
