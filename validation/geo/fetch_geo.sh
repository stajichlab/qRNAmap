#!/usr/bin/env bash
# fetch_geo.sh: download the GSE80357 series tables and the GEO/ENA records used for the mapping.
# Runs on the login node (HTTPS, under 1 MB in total). The downloaded files are untrusted data:
# they are only read by prepare_geo.py (python -I) and compare_geo.R, never executed.
#
# Usage: bash /ABS/validation/geo/fetch_geo.sh GEO_DIR
#   GEO_DIR  absolute path, for example /bigdata/stajichlab/jstajich/projects/qRNAmap_data/amb/geo
# Writes GEO_DIR/download/*.csv.gz, GEO_DIR/download/provenance.tsv, GEO_DIR/gsm/GSM*.txt,
# GEO_DIR/gse80357_series_full.txt and GEO_DIR/srx_to_srr.tsv.
# Note: on 2026-11-20 these commands were first run by hand in the shell with the same URLs;
# this script holds the same commands.
set -euo pipefail
GEO="${1:?GEO_DIR}"
mkdir -p "$GEO/download" "$GEO/gsm"
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80357/suppl
printf 'file\turl\tbytes\tmd5\tlast_modified\tdownloaded\n' > "$GEO/download/provenance.tsv"
for f in GSE80357_counts_of_mapped_reads.csv.gz GSE80357_log_fold_change_and_FDR_values.csv.gz \
         GSE80357_average_cpm_per_treatment_group.csv.gz; do
    lm=$(curl -sS -I "$BASE/$f" | tr -d '\r' | awk -F': ' 'tolower($1)=="last-modified"{print $2}')
    curl -sS -f -o "$GEO/download/$f" "$BASE/$f"
    gzip -t "$GEO/download/$f"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$BASE/$f" "$(stat -c %s "$GEO/download/$f")" \
        "$(md5sum < "$GEO/download/$f" | cut -d' ' -f1)" "$lm" "$(date -Is)" \
        >> "$GEO/download/provenance.tsv"
done
curl -sS -f -o "$GEO/gse80357_series_full.txt" \
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80357&targ=self&form=text&view=full"
for i in $(seq 2125072 2125089); do
    curl -sS -f -o "$GEO/gsm/GSM$i.txt" \
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM$i&targ=self&form=text&view=full"
    sleep 0.4
done
: > "$GEO/srx_to_srr.tsv"
for x in SRX1709853 SRX1709854 SRX1709855 SRX1709856 SRX1709857 SRX1709858; do
    curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$x&result=read_run&fields=run_accession,experiment_accession,sample_title&format=tsv" \
        | tail -n +2 >> "$GEO/srx_to_srr.tsv"
done
cat "$GEO/download/provenance.tsv"
