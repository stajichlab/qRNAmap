#!/usr/bin/env python3
"""Compare two transcript FASTA files by record id and sequence.

Usage:
  compare_gffread.py --ours A.fa[.gz|.zst] --gffread B.fa[.gz|.zst] [--gff G.gff3[.gz|.zst]] --out diff.tsv

The id is the first word of the header line. Sequences are compared after upper-casing.
Every record that is missing on one side or has a different sequence is written to --out
(one row per record). With --gff, each row also gets the GFF3 feature type, the gene_biotype
of the parent gene, and the number of exon and CDS children of that id, so that each
difference can be explained.

Input files can be plain, gzip (.gz) or zstd (.zst). zstd uses the `zstandard` module when
it is present, else the `zstd -dc` command.
Exit status: 0 when the comparison ran (even with differences), 1 on an input error.
"""
import argparse
import collections
import gzip
import io
import subprocess
import sys
import urllib.parse


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        try:
            import zstandard  # type: ignore

            fh = open(path, "rb")
            return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(fh))
        except ImportError:
            proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)
            return io.TextIOWrapper(proc.stdout)
    return open(path, "rt")


def read_fasta(path):
    seqs = {}
    dups = []
    name = None
    parts = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0] if len(line) > 1 else ""
                if name in seqs:
                    dups.append(name)
                parts = []
            elif line:
                parts.append(line.strip())
        if name is not None:
            seqs[name] = "".join(parts).upper()
    return seqs, dups


def read_gff(path):
    """Return id -> (type, parent id), gene id -> gene_biotype, id -> Counter of child types."""
    info = {}
    biotype = {}
    children = collections.defaultdict(collections.Counter)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9:
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = urllib.parse.unquote(v)
            fid = attrs.get("ID")
            parents = attrs.get("Parent", "").split(",") if attrs.get("Parent") else []
            if fid and fid not in info:
                info[fid] = (f[2], parents[0] if parents else "")
            if fid and "gene_biotype" in attrs:
                biotype[fid] = attrs["gene_biotype"]
            for p in parents:
                children[p][f[2]] += 1
    return info, biotype, children


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ours", required=True)
    ap.add_argument("--gffread", required=True)
    ap.add_argument("--gff")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ours, dup_o = read_fasta(a.ours)
    gff_r, dup_g = read_fasta(a.gffread)
    info, biotype, children = read_gff(a.gff) if a.gff else ({}, {}, {})

    def annot(i):
        t, parent = info.get(i, ("", ""))
        bt = biotype.get(parent, biotype.get(i, ""))
        ch = children.get(i, collections.Counter()) if children else collections.Counter()
        return [t, parent, bt, str(ch.get("exon", 0)), str(ch.get("CDS", 0))]

    rows = []
    for i in sorted(set(ours) | set(gff_r)):
        o, g = ours.get(i), gff_r.get(i)
        if o is None:
            status = "only_gffread"
        elif g is None:
            status = "only_ours"
        elif o != g:
            status = "seq_differs"
        else:
            continue
        rows.append([i, status, str(len(o)) if o is not None else "", str(len(g)) if g is not None else ""] + annot(i))

    with open(a.out, "w") as out:
        out.write("id\tstatus\tlen_ours\tlen_gffread\tgff_type\tparent\tgene_biotype\texons\tcds\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    same = sum(1 for i in ours if i in gff_r and ours[i] == gff_r[i])
    counts = collections.Counter(r[1] for r in rows)
    print(f"records ours={len(ours)} gffread={len(gff_r)} identical={same}")
    print(f"duplicate ids ours={len(dup_o)} gffread={len(dup_g)}")
    for k in ("only_ours", "only_gffread", "seq_differs"):
        print(f"{k}={counts.get(k, 0)}")
    by_type = collections.Counter((r[1], r[4], r[6]) for r in rows)
    for (st, t, bt), n in sorted(by_type.items()):
        print(f"  {st}\ttype={t or '-'}\tbiotype={bt or '-'}\t{n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
