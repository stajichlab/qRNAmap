"""prepare_geo.py: map the GSE80357 count columns to the AMB samples and build the inputs.

Run with the pixi python and `-I` (the GEO files are untrusted data; this script only parses them):

    python -I prepare_geo.py --geo-dir DIR --srx-srr FILE --sheet FILE \
        --quant-counts FILE --align-counts FILE --out-dir DIR --results-dir DIR

Inputs (all absolute paths):
  --geo-dir       directory with download/GSE80357_counts_of_mapped_reads.csv.gz and
                  gsm/GSM*.txt (GEO sample records, text form, view=full)
  --srx-srr       TSV run_accession, experiment_accession, sample_title (ENA filereport of the SRX)
  --sheet         the AMB sample sheet (sample, fastq1, ...); the run id is the fastq1 basename
                  before the first "_"
  --quant-counts, --align-counts   qrnamap counts.tsv.zst (gene x sample)
Outputs:
  OUT_DIR/authors_counts_amb.tsv    gene + the 6 AMB sample columns (our names), all authors' rows
  OUT_DIR/unmatched_ours.tsv, OUT_DIR/unmatched_theirs.tsv   full lists with read sums
  RESULTS_DIR/column_mapping.tsv    the mapping chain with evidence
  RESULTS_DIR/gene_id_overlap.tsv   set sizes and read totals per sample
  RESULTS_DIR/unmatched_summary.tsv per name class: genes, genes with reads, reads
  RESULTS_DIR/unmatched_top.tsv     the 15 unmatched genes with most reads on each side
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys


def read_gsm(path):
    """Return title, list of descriptions, SRX accession from a GEO sample text record."""
    title, descs, srx, gsm = None, [], None, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                gsm = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_title = "):
                title = line.split(" = ", 1)[1].strip()
            elif line.startswith("!Sample_description = "):
                descs.append(line.split(" = ", 1)[1].strip())
            elif line.startswith("!Sample_relation = SRA: "):
                m = re.search(r"term=(SRX\d+)", line)
                if m:
                    srx = m.group(1)
    return gsm, title, descs, srx


def read_matrix_zst(path):
    out = subprocess.run(["zstd", "-dc", path], check=True, capture_output=True).stdout.decode()
    rows = out.splitlines()
    header = rows[0].split("\t")
    data = {}
    for r in rows[1:]:
        f = r.split("\t")
        data[f[0]] = [int(x) for x in f[1:]]
    return header[1:], data


def name_class(g):
    if re.fullmatch(r"Q\d{4}", g):
        return "mitochondrial ORF (Q####)"
    if re.fullmatch(r"t[A-Z]\([ACGU]{3}\)[A-Z0-9]*", g):
        return "tRNA (t?(NNN)?)"
    if re.fullmatch(r"Y[A-P][LR]\d{3}[WC](-[A-Z])?", g):
        return "systematic ORF name (Y??###W/C)"
    if re.fullmatch(r"YNC[A-Q]\d{4}[WC]", g):
        return "SGD ncRNA systematic name (YNC?####W/C)"
    if re.fullmatch(r"(snR|SNR)\d+[a-z]?", g):
        return "snoRNA (snR)"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    for a in ("--geo-dir", "--srx-srr", "--sheet", "--quant-counts", "--align-counts",
              "--out-dir", "--results-dir"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)

    # Authors' counts file.
    cpath = os.path.join(args.geo_dir, "download", "GSE80357_counts_of_mapped_reads.csv.gz")
    with gzip.open(cpath, "rt", encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        theirs_all = {}
        for row in rd:
            if len(row) != len(header):
                sys.exit(f"bad row width in {cpath}: {row[:3]}")
            if row[0] in theirs_all:
                sys.exit(f"duplicate gene id {row[0]}")
            theirs_all[row[0]] = [int(x) for x in row[1:]]
    print(f"authors' counts: {len(theirs_all)} genes, columns {header}")

    # GEO sample records: description code -> GSM, title, SRX.
    srx_to_srr = {}
    with open(args.srx_srr, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].startswith("SRX"):
                srx_to_srr[f[1]] = f[0]
    run_to_sample = {}
    with open(args.sheet, encoding="utf-8") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            run = os.path.basename(r["fastq1"]).split("_")[0]
            run_to_sample[run] = r["sample"]
    gsm_dir = os.path.join(args.geo_dir, "gsm")
    code_rows = {}
    for fn in sorted(os.listdir(gsm_dir)):
        if not fn.endswith(".txt"):
            continue
        gsm, title, descs, srx = read_gsm(os.path.join(gsm_dir, fn))
        code = descs[0] if descs else None
        code_rows[code.upper()] = (gsm, title, code, srx)

    mapping = []
    for col in header[1:]:
        if col not in code_rows:
            sys.exit(f"column {col} has no GSM with that description code")
        gsm, title, code, srx = code_rows[col]
        srr = srx_to_srr.get(srx, "")
        sample = run_to_sample.get(srr, "")
        mapping.append((col, gsm, title, code, srx, srr, sample))
    with open(os.path.join(args.results_dir, "column_mapping.tsv"), "w", encoding="utf-8") as fh:
        fh.write("authors_column\tgsm\tgsm_title\tgsm_description_code\tsrx\tsrr\tour_sample\n")
        for m in mapping:
            fh.write("\t".join(m) + "\n")
    col_of = {m[6]: header.index(m[0]) - 1 for m in mapping if m[6]}

    # Our matrices.
    samples, quant = read_matrix_zst(args.quant_counts)
    samples_a, align = read_matrix_zst(args.align_counts)
    if samples != samples_a or set(quant) != set(align):
        sys.exit("quant and align matrices differ in samples or genes")
    if sorted(col_of) != sorted(samples):
        sys.exit(f"mapped samples {sorted(col_of)} != matrix samples {sorted(samples)}")
    theirs = {g: [v[col_of[s]] for s in samples] for g, v in theirs_all.items()}

    with open(os.path.join(args.out_dir, "authors_counts_amb.tsv"), "w", encoding="utf-8") as fh:
        fh.write("gene\t" + "\t".join(samples) + "\n")
        for g, v in theirs.items():
            fh.write(g + "\t" + "\t".join(str(x) for x in v) + "\n")

    ours_ids, their_ids = set(quant), set(theirs)
    inter = ours_ids & their_ids
    only_ours = ours_ids - their_ids
    only_theirs = their_ids - ours_ids

    def totals(ids, table):
        return [sum(table[g][i] for g in ids) for i in range(len(samples))]

    with open(os.path.join(args.results_dir, "gene_id_overlap.tsv"), "w", encoding="utf-8") as fh:
        fh.write("set\ttable\tgenes\t" + "\t".join(samples) + "\n")
        for name, ids, tabs in (("all", ours_ids, (("quant", quant), ("align", align))),
                                ("all", their_ids, (("authors", theirs),)),
                                ("intersection", inter,
                                 (("quant", quant), ("align", align), ("authors", theirs))),
                                ("only_ours", only_ours, (("quant", quant), ("align", align))),
                                ("only_theirs", only_theirs, (("authors", theirs),))):
            for tname, tab in tabs:
                fh.write(f"{name}\t{tname}\t{len(ids)}\t" +
                         "\t".join(str(x) for x in totals(ids, tab)) + "\n")

    # Unmatched lists and summaries.
    def write_list(path, ids, tabs):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("gene\tclass\t" + "\t".join(f"{t}_sum6" for t, _ in tabs) + "\n")
            for g in sorted(ids, key=lambda g: -sum(tabs[0][1][g])):
                fh.write(f"{g}\t{name_class(g)}\t" +
                         "\t".join(str(sum(tab[g])) for _, tab in tabs) + "\n")
    write_list(os.path.join(args.out_dir, "unmatched_ours.tsv"), only_ours,
               (("quant", quant), ("align", align)))
    write_list(os.path.join(args.out_dir, "unmatched_theirs.tsv"), only_theirs,
               (("authors", theirs),))

    with open(os.path.join(args.results_dir, "unmatched_summary.tsv"), "w",
              encoding="utf-8") as fh:
        fh.write("side\ttable\tclass\tgenes\tgenes_with_reads\treads_sum6\n")
        for side, ids, tabs in (("only_ours", only_ours, (("quant", quant), ("align", align))),
                                ("only_theirs", only_theirs, (("authors", theirs),))):
            classes = sorted({name_class(g) for g in ids})
            for tname, tab in tabs:
                for c in classes:
                    gs = [g for g in ids if name_class(g) == c]
                    fh.write(f"{side}\t{tname}\t{c}\t{len(gs)}\t"
                             f"{sum(1 for g in gs if sum(tab[g]) > 0)}\t"
                             f"{sum(sum(tab[g]) for g in gs)}\n")
    with open(os.path.join(args.results_dir, "unmatched_top.tsv"), "w", encoding="utf-8") as fh:
        fh.write("side\tgene\tclass\treads_sum6_primary\treads_sum6_other\n")
        for g in sorted(only_ours, key=lambda g: -sum(quant[g]))[:15]:
            fh.write(f"only_ours\t{g}\t{name_class(g)}\tquant={sum(quant[g])}\t"
                     f"align={sum(align[g])}\n")
        for g in sorted(only_theirs, key=lambda g: -sum(theirs[g]))[:15]:
            fh.write(f"only_theirs\t{g}\t{name_class(g)}\tauthors={sum(theirs[g])}\t-\n")
    print(f"ours {len(ours_ids)}, theirs {len(their_ids)}, intersection {len(inter)}, "
          f"only ours {len(only_ours)}, only theirs {len(only_theirs)}")


if __name__ == "__main__":
    main()
