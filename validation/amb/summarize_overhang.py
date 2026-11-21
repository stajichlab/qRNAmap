#!/usr/bin/env python3
"""summarize_overhang.py: one table of the Salmon settings on AMB (default, A, B), from the small
TSVs that run_amb_salmon.sh and run_amb_salmon_overhang.sh wrote. It computes nothing new except
the pooled mapped fraction (sum num_mapped / sum num_processed) and min/max over samples.

Usage: python -I summarize_overhang.py RESULTS_AMB_DIR OUT.tsv
  RESULTS_AMB_DIR is validation/amb. Inputs: results_salmon/ (default),
  results_salmon_overhang/A/ and results_salmon_overhang/B/.
Output columns: setting, metric, value. Per-sample metrics have the sample in the metric name.
"""
import csv
import os
import sys

SETTINGS = [("default", "results_salmon"), ("A", "results_salmon_overhang/A"),
            ("B", "results_salmon_overhang/B")]


def rd(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    base, out = sys.argv[1], sys.argv[2]
    rows = []
    for name, sub in SETTINGS:
        d = os.path.join(base, sub)
        add = lambda m, v: rows.append({"setting": name, "metric": m, "value": v})  # noqa: E731
        ps = rd(os.path.join(d, "per_sample.tsv"))
        for r in ps:
            add(f"percent_mapped:{r['sample']}", f"{float(r['percent_mapped']):.2f}")
            add(f"num_mapped:{r['sample']}", r["num_mapped"])
        add("percent_mapped:pooled",
            f"{100 * sum(int(r['num_mapped']) for r in ps) / sum(int(r['num_processed']) for r in ps):.2f}")
        add("sample_task_seconds_sum", sum(int(r["seconds"]) for r in ps))
        for r in rd(os.path.join(d, "mode_agreement.tsv")):
            s = "pooled" if r["sample"] == "ALL_pooled" else r["sample"]
            add(f"spearman_vs_kallisto:{s}", r["spearman_all_genes"])
            add(f"zero_genes_salmon:{s}", r["zero_quant"])
            add(f"zero_genes_kallisto:{s}", r["zero_align"])
            add(f"total_salmon:{s}", r["total_quant"])
        de = {r["metric"]: r["value"] for r in rd(os.path.join(d, "de_overlap.tsv"))}
        for k in ("sig_quant", "sig_align", "sig_both", "sig_only_quant", "sig_only_align", "jaccard"):
            add(f"de_vs_kallisto:{k}", de[k])
        for r in rd(os.path.join(d, "zero_genes_by_length.tsv")):
            add(f"zero_salmon_nonzero_kallisto:{r['sample']}", r["salmon_zero_kallisto_nonzero"])
            add(f"of_these_all_tx_le_read_length:{r['sample']}", r["all_transcripts_le_read_length"])
            add(f"kallisto_reads_in_these_genes:{r['sample']}", r["kallisto_reads_in_these_genes"])
        for r in rd(os.path.join(d, "zero_genes_examples.tsv")):
            add(f"count_{r['name']}:{r['sample']}", r["salmon"])
        for r in rd(os.path.join(d, "geo", "counts_agreement.tsv")):
            if r["sample"] == "pooled" and r["comparison"] == "authors(x)_vs_quant(y)":
                add("authors_counts_spearman:pooled", r["spearman_raw"])
                add("authors_counts_pearson_log2p1:pooled", r["pearson_log2p1"])
                add("authors_counts_median_ratio:pooled", r["ratio_median"])
                add("authors_counts_total_ratio:pooled", r["total_ratio_y_over_x"])
            elif r["comparison"] == "authors(x)_vs_quant(y)":
                add(f"authors_counts_spearman:{r['sample']}", r["spearman_raw"])
        for r in rd(os.path.join(d, "geo", "de_agreement.tsv")):
            if r["comparison"] == "authors_published(x)_vs_deseq2_quant(y)":
                for k in ("sig_x_common", "sig_y_common", "sig_both", "jaccard", "spearman_lfc"):
                    if k in r:
                        add(f"authors_published_de:{k}", r[k])
            if r["comparison"] == "ii:deseq2_authors_counts(x)_vs_deseq2_quant(y)":
                add("deseq2_on_authors_counts_de:jaccard", r["jaccard"])
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["setting", "metric", "value"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
