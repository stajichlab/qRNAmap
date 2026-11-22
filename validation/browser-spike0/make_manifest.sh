#!/usr/bin/env bash
# make_manifest.sh FIXTURE_DIR : write FIXTURE_DIR/manifest.tsv (file, bytes, sha256, how made) for every
# file under FIXTURE_DIR except manifest.tsv. The "how made" text comes from how_made.tsv (pattern<TAB>text,
# first matching glob wins) in the directory given by $SCRIPTS. No BASH_SOURCE.
set -euo pipefail
dir=${1:?FIXTURE_DIR}; : "${SCRIPTS:?}"
out="$dir/manifest.tsv"
printf 'file\tbytes\tsha256\thow_made\n' > "$out.tmp"
cd "$dir"
find . -type f ! -name manifest.tsv ! -name 'manifest.tsv.tmp' | sed 's#^\./##' | LC_ALL=C sort | while read -r f; do
    how=$(awk -F'\t' -v f="$f" 'BEGIN{h="UNSPECIFIED"} { p=$1; gsub(/\./,"\\.",p); gsub(/\*/,".*",p); if (f ~ "^" p "$") { h=$2; exit } } END{print h}' "$SCRIPTS/how_made.tsv")
    printf '%s\t%s\t%s\t%s\n' "$f" "$(stat -c %s "$f")" "$(sha256sum "$f" | cut -d' ' -f1)" "$how" >> "$out.tmp"
done
mv "$out.tmp" "$out"; wc -l "$out"
