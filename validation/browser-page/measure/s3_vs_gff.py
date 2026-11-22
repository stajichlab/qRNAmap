#!/usr/bin/env python3
"""Compare the gene identifiers of Supplementary Table S3 (PMC11106425, sheet lkhl_basal_vs_H99_basal, dumped to TSV)
with the RefSeq H99 GFF3 (GCF_000149245.1) and, optionally, with a qrnamap gene table.

Usage: s3_vs_gff.py S3_BASAL_TSV GFF3_GZ [GENES_TSV]
Prints one TSV row per check. Stdlib only. Inputs are data; nothing in them is run.
"""
import gzip, sys, re

s3_path, gff_path = sys.argv[1], sys.argv[2]
genes_path = sys.argv[3] if len(sys.argv) > 3 else None

rows = [line.rstrip('\n').split('\t') for line in open(s3_path, encoding='utf-8')]
hdr_i = next(i for i, r in enumerate(rows) if r and r[0] == 'Gene_ID')
hdr = rows[hdr_i]
gi, lt, pi = hdr.index('Gene_ID'), hdr.index('locus_tag'), hdr.index('Protein_ID')
bt = hdr.index('gene_biotype')
fc = hdr.index('lkhl_basal/H99_basal.fc')
s3 = [r for r in rows[hdr_i + 1:] if len(r) > lt and r[gi].strip()]

gff_gene = {}  # locus_tag -> GeneID
gff_prot = set()
with gzip.open(gff_path, 'rt', encoding='utf-8') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = f[8]
        if f[2] == 'gene':
            m_lt = re.search(r'(?:^|;)locus_tag=([^;]+)', a)
            m_id = re.search(r'GeneID:(\d+)', a)
            if m_lt:
                gff_gene[m_lt.group(1)] = m_id.group(1) if m_id else ''
        elif f[2] == 'CDS':
            m = re.search(r'(?:^|;)protein_id=([^;]+)', a)
            if m:
                gff_prot.add(m.group(1))

lt_in = sum(1 for r in s3 if r[lt] in gff_gene)
pair_ok = sum(1 for r in s3 if gff_gene.get(r[lt]) == r[gi])
prot_ok = sum(1 for r in s3 if r[pi] in gff_prot)
up = sum(1 for r in s3 if float(r[fc]) > 0)
down = sum(1 for r in s3 if float(r[fc]) < 0)
print('check\tvalue')
print(f's3_rows\t{len(s3)}')
print(f's3_fc_positive\t{up}')
print(f's3_fc_negative\t{down}')
print(f'gff_genes_with_locus_tag\t{len(gff_gene)}')
print(f's3_locus_tag_in_gff\t{lt_in}')
print(f's3_locus_tag_and_GeneID_pair_equal_gff\t{pair_ok}')
print(f's3_protein_id_in_gff_cds\t{prot_ok}')
for biotype in sorted({r[bt] for r in s3}):
    print(f's3_gene_biotype_{biotype}\t{sum(1 for r in s3 if r[bt] == biotype)}')
print(f's3_protein_id_cells_with_semicolon\t{sum(1 for r in s3 if ";" in r[pi])}')
missing = [r[lt] for r in s3 if r[lt] not in gff_gene]
print(f's3_locus_tag_missing_examples\t{",".join(missing[:10])}')
if genes_path:
    ids = set()
    for line in open(genes_path, encoding='utf-8'):
        ids.add(line.split('\t', 1)[0].strip())
    print(f'genes_table_rows\t{len(ids)}')
    print(f's3_locus_tag_in_genes_table\t{sum(1 for r in s3 if r[lt] in ids)}')
    print(f's3_GeneID_in_genes_table\t{sum(1 for r in s3 if r[gi] in ids)}')
