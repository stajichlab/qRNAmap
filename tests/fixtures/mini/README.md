# Mini fixture: transcript builder, SAF and splice sites

Files:
- `genome.fa`: two contigs, `chrA` and `chrB`, 1500 bp each (3000 bp in total). Bases come from a fixed linear congruential generator. Positions 141-180 of `chrA` are lowercase (soft-masked) to test that transcripts are uppercase (ruling R12).
- `genes.gff3`: five genes (below).
- `make_fixture.py`: writes both files and prints the expected transcript sequences. It is an independent Python implementation (slice, join, reverse-complement). Run: `python3 make_fixture.py <dir>`.

Canonical splice motifs were written into the introns (plus strand GT...AG; minus strand CT...AC in plus-strand view). No tool in these tests checks them.

## Genes and expected transcripts

Coordinates are 1-based inclusive. Sequences are in transcript orientation, uppercase.

| Gene | Contig | Strand | Transcript | Kind | Intervals used | Rule |
|---|---|---|---|---|---|---|
| geneA | chrA | + | tA1 | mRNA | exons 101-160, 301-350 | exons (CDS 111-160, 301-340 ignored) |
| geneB | chrA | - | tB1 | mRNA | exons 601-640, 701-740, 801-840 | exons |
| geneB | chrA | - | tB2 | mRNA | exons 621-640, 701-740, 821-870 | exons |
| geneR | chrB | + | tR1 | rRNA | exon 101-180 | exons; rRNA kept |
| geneC | chrB | - | tC1 | mRNA | CDS 401-430, 501-540 | no exons, no UTR: CDS only |
| geneU | chrB | + | tU1 | mRNA | 5'UTR 701-720 + CDS 721-760 + 3'UTR 761-780 = 701-780 | no exons: CDS merged with UTRs |

Expected sequences (from `make_fixture.py`):

```
tA1  AAGTGCCTAGCATTACGGTCAAATCTGGTTATACTCGTTGCGTGTACTTCACGGATTATGCATGTGTCGCTTACGATCGATCGTTCATCGTGTGCTCGTATGCCAACGAT
tB1  TCCTCTACATGGCAATTGCCGTTTCAGCTTTCGCTGATCTTGCCGCATCTTCTTGAGCAATATAGCAGTGAAAAAACCGATCCCCAAAGAGTACAAGTGAACCATTATCAGAGAGATGCT
tB2  GAGACTGGAGCGTAGCCCACAGCGCAAGTATCCTCTACATGGCAATTGCCTGCCGCATCTTCTTGAGCAATATAGCAGTGAAAAAACCGATCCCCAAAGAGTACAAGTGA
tR1  TCGAAGCGCGTTTTAAATTCTGTTGGCGCGCCGACACGGTACCGAAACGTCCAAATTGAGCAATAATCCTCTGATCCAAG
tC1  ATAAAGTAACCACGGAGTAACCGTTCCGCCCCCGCGCGTGGTTCAACTCCGATGCGAAGGTTGGCGGCGT
tU1  GTAATAAACGCTTCAGCACTGTTGACAACCACCTCGATTGGCACATGGCGCAAATTAGGCGCGATGTCCCCCGGTGGGGC
```

Lengths: tA1 110, tB1 120, tB2 110, tR1 80, tC1 70, tU1 80.

## Expected SAF (union of exon sets per gene)

```
GeneID  Chr   Start  End  Strand
geneA   chrA  101    160  +
geneA   chrA  301    350  +
geneB   chrA  601    640  -
geneB   chrA  701    740  -
geneB   chrA  801    870  -    (tB1 801-840 and tB2 821-870 merge)
geneR   chrB  101    180  +
geneC   chrB  401    430  -
geneC   chrB  501    540  -
geneU   chrB  701    780  +
```

## Expected HISAT2 splice sites

`left` = 0-based position of the last base of the upstream exon; `right` = 0-based position of the first base of the downstream exon.

```
chrA  159  300  +    tA1
chrA  639  700  -    tB1 and tB2 (shared intron, one line)
chrA  739  800  -    tB1
chrA  739  820  -    tB2
chrB  429  500  -    tC1 (CDS-only gene)
```

## Comparison with gffread 0.12.9 and HISAT2 2.2.2 (`hisat2_extract_splice_sites.py`)

Run by `crates/qrnamap-core/tests/gffread_parity.rs`.

1. `gffread -g genome.fa -w out.fa genes.gff3`: all 6 transcripts match ours by id and sequence after uppercasing. No disagreement.
   - gffread keeps soft-masked bases lowercase; we write uppercase (R12). The test uppercases gffread's sequence before comparing.
   - gffread uses CDS as the exon set of the CDS-only gene tC1, and merges the UTR and CDS features of tU1 into one exon 701-780. Both agree with our rule (R4).
   - gffread keeps the rRNA transcript tR1. Agrees with our rule.
2. `gffread -T --force-exons` GTF, then `hisat2_extract_splice_sites.py`: the same 5 lines as ours. No disagreement.
3. Plain `gffread -T` GTF, then `hisat2_extract_splice_sites.py`: HISAT2 gives 4 lines; ours has one more, `chrB 429 500 -` (geneC).
   - Cause: plain `gffread -T` writes only `CDS` lines (no `exon` lines) for the CDS-only gene geneC, and the HISAT2 script reads `exon` lines only.
   - Decision: our rule is right, the difference is expected. R4 makes the CDS the exon set of a CDS-only transcript, and the transcript sequence (which gffread `-w` also builds from the CDS) has this intron spliced out. HISAT2 should know this junction. The test asserts that this line is the only difference.

## Rules copied from `hisat2_extract_splice_sites.py`

- Exons whose start minus the previous exon end is 5 or less are joined; no junction is written for that gap. `splice_sites` does the same.
- Output is sorted by (chrom, left, right, strand), chrom in byte order.

Known difference not covered by the fixture: the HISAT2 script drops single-base exons (`left >= right`); we keep them. This affects only GFF3 files with 1 bp exons.
