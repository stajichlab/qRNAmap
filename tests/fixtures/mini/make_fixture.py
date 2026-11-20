#!/usr/bin/env python3
"""Write the mini fixture (genome.fa, genes.gff3) and print the expected transcripts.

Deterministic: a fixed linear congruential generator makes the bases.
The expected sequences are computed here, independently of the Rust code.
Run: python3 make_fixture.py <outdir>
"""
import sys


def lcg_bases(n, seed):
    x = seed
    out = []
    for _ in range(n):
        x = (1103515245 * x + 12345) % (1 << 31)
        out.append("ACGT"[(x >> 16) & 3])
    return out


def setb(s, pos1, text):
    """Write text at 1-based position pos1."""
    for i, c in enumerate(text):
        s[pos1 - 1 + i] = c


LEN = {"chrA": 1500, "chrB": 1500}
chrA = lcg_bases(LEN["chrA"], 7)
chrB = lcg_bases(LEN["chrB"], 11)

# Canonical splice motifs. Plus-strand intron: GT...AG. Minus-strand intron (plus-strand view): CT...AC.
setb(chrA, 161, "GT")
setb(chrA, 299, "AG")  # geneA intron 161-300
for a, b in [(641, 700), (741, 800), (741, 820)]:  # geneB introns
    setb(chrA, a, "CT")
    setb(chrA, b - 1, "AC")
setb(chrB, 431, "CT")
setb(chrB, 499, "AC")  # geneC intron 431-500

chrA = "".join(chrA)
chrB = "".join(chrB)
# Soft-masked stretch inside geneA exon 1 (R12): positions 141-180 lowercase.
chrA = chrA[:140] + chrA[140:180].lower() + chrA[180:]
genome = {"chrA": chrA, "chrB": chrB}

GFF = """##gff-version 3
##sequence-region chrA 1 1500
##sequence-region chrB 1 1500
chrA\tmini\tgene\t101\t350\t.\t+\t.\tID=geneA;Name=geneA
chrA\tmini\tmRNA\t101\t350\t.\t+\t.\tID=tA1;Parent=geneA
chrA\tmini\texon\t101\t160\t.\t+\t.\tParent=tA1
chrA\tmini\texon\t301\t350\t.\t+\t.\tParent=tA1
chrA\tmini\tCDS\t111\t160\t.\t+\t0\tParent=tA1
chrA\tmini\tCDS\t301\t340\t.\t+\t1\tParent=tA1
chrA\tmini\tgene\t601\t870\t.\t-\t.\tID=geneB;Name=geneB
chrA\tmini\tmRNA\t601\t840\t.\t-\t.\tID=tB1;Parent=geneB
chrA\tmini\texon\t601\t640\t.\t-\t.\tParent=tB1
chrA\tmini\texon\t701\t740\t.\t-\t.\tParent=tB1
chrA\tmini\texon\t801\t840\t.\t-\t.\tParent=tB1
chrA\tmini\tmRNA\t621\t870\t.\t-\t.\tID=tB2;Parent=geneB
chrA\tmini\texon\t621\t640\t.\t-\t.\tParent=tB2
chrA\tmini\texon\t701\t740\t.\t-\t.\tParent=tB2
chrA\tmini\texon\t821\t870\t.\t-\t.\tParent=tB2
chrB\tmini\tgene\t101\t180\t.\t+\t.\tID=geneR;Name=geneR
chrB\tmini\trRNA\t101\t180\t.\t+\t.\tID=tR1;Parent=geneR
chrB\tmini\texon\t101\t180\t.\t+\t.\tParent=tR1
chrB\tmini\tgene\t401\t540\t.\t-\t.\tID=geneC;Name=geneC
chrB\tmini\tmRNA\t401\t540\t.\t-\t.\tID=tC1;Parent=geneC
chrB\tmini\tCDS\t401\t430\t.\t-\t0\tParent=tC1
chrB\tmini\tCDS\t501\t540\t.\t-\t0\tParent=tC1
chrB\tmini\tgene\t701\t780\t.\t+\t.\tID=geneU;Name=geneU
chrB\tmini\tmRNA\t701\t780\t.\t+\t.\tID=tU1;Parent=geneU
chrB\tmini\tfive_prime_UTR\t701\t720\t.\t+\t.\tParent=tU1
chrB\tmini\tCDS\t721\t760\t.\t+\t0\tParent=tU1
chrB\tmini\tthree_prime_UTR\t761\t780\t.\t+\t.\tParent=tU1
"""

# transcript id -> (contig, strand, intervals), written by hand, not read from the GFF parser
TX = {
    "tA1": ("chrA", "+", [(101, 160), (301, 350)]),
    "tB1": ("chrA", "-", [(601, 640), (701, 740), (801, 840)]),
    "tB2": ("chrA", "-", [(621, 640), (701, 740), (821, 870)]),
    "tR1": ("chrB", "+", [(101, 180)]),
    "tC1": ("chrB", "-", [(401, 430), (501, 540)]),
    "tU1": ("chrB", "+", [(701, 780)]),
}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def tx_seq(contig, strand, ivs):
    s = "".join(genome[contig][a - 1:b] for a, b in ivs).upper()
    if strand == "-":
        s = "".join(COMP[c] for c in reversed(s))
    return s


def main():
    out = sys.argv[1]
    with open(f"{out}/genome.fa", "w") as f:
        for name, seq in genome.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    with open(f"{out}/genes.gff3", "w") as f:
        f.write(GFF)
    for tid, (c, st, ivs) in TX.items():
        print(f"{tid}\t{tx_seq(c, st, ivs)}")


main()
