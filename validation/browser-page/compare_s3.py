#!/usr/bin/env python3
"""Compare C. neoformans preview-depth CLI DESeq2 calls against paper Supplementary Table S3.

Criterion: adjusted p < 0.05 and |log2 fold change| >= 2.0.
"""
import csv
import sys

def main():
    s3_tsv = "/bigdata/stajichlab/jstajich/projects/qRNAmap_data/page-v1/dl_s3/tsv/lkhl_basal_vs_H99_basal.tsv"
    cli_csv = "validation/browser-page/expected/cneo-lkh1/deseq2.results.csv"

    # 1. Parse Table S3 genes
    s3_genes = set()
    s3_directions = {}
    with open(s3_tsv, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        for row in reader:
            if not row:
                continue
            if row[0] == "Gene_ID":
                header = row
                break
        assert header, "Header row starting with Gene_ID not found in S3 TSV"
        lt_idx = header.index("locus_tag")
        fc_idx = header.index("lkhl_basal/H99_basal.fc")

        for row in reader:
            if len(row) > lt_idx and row[lt_idx].strip():
                gene = row[lt_idx].strip()
                s3_genes.add(gene)
                try:
                    fc = float(row[fc_idx])
                    s3_directions[gene] = 1 if fc > 0 else -1
                except ValueError:
                    pass

    print(f"Paper Table S3 basal DE genes: {len(s3_genes)}")

    # 2. Parse CLI DESeq2 results
    cli_sig_genes = set()
    cli_directions = {}
    total_genes = 0
    with open(cli_csv, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_genes += 1
            gene = row["gene"].strip('"')
            padj_str = row["padj"]
            lfc_str = row["log2FoldChange"]
            if padj_str in ("NA", "") or lfc_str in ("NA", ""):
                continue
            try:
                padj = float(padj_str)
                lfc = float(lfc_str)
            except ValueError:
                continue

            if padj < 0.05 and abs(lfc) >= 2.0:
                cli_sig_genes.add(gene)
                cli_directions[gene] = 1 if lfc > 0 else -1

    print(f"Total genes tested in CLI: {total_genes}")
    print(f"CLI preview (N=75,000) DE genes (padj < 0.05, |log2FC| >= 2): {len(cli_sig_genes)}")

    # 3. Compare overlap
    overlap = cli_sig_genes & s3_genes
    cli_only = cli_sig_genes - s3_genes
    s3_only = s3_genes - cli_sig_genes

    concordant_dir = sum(1 for g in overlap if cli_directions.get(g) == s3_directions.get(g))

    print(f"Overlap (in both): {len(overlap)}")
    print(f"Direction concordance in overlap: {concordant_dir}/{len(overlap)}")
    print(f"CLI only: {len(cli_only)}")
    print(f"Table S3 only: {len(s3_only)}")
    jaccard = len(overlap) / len(cli_sig_genes | s3_genes) if (cli_sig_genes | s3_genes) else 0
    print(f"Jaccard similarity: {jaccard:.4f}")

    # Output report TSV
    out_tsv = "validation/browser-page/expected/cneo-lkh1/paper_s3_comparison.tsv"
    with open(out_tsv, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_genes\t{total_genes}\n")
        f.write(f"paper_s3_genes\t{len(s3_genes)}\n")
        f.write(f"cli_preview_de_genes\t{len(cli_sig_genes)}\n")
        f.write(f"overlap_genes\t{len(overlap)}\n")
        f.write(f"concordant_direction_genes\t{concordant_dir}\n")
        f.write(f"cli_only_genes\t{len(cli_only)}\n")
        f.write(f"s3_only_genes\t{len(s3_only)}\n")
        f.write(f"jaccard_index\t{jaccard:.4f}\n")
    print(f"Wrote summary to {out_tsv}")

if __name__ == "__main__":
    main()
