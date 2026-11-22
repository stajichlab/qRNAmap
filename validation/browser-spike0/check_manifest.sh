#!/usr/bin/env bash
# check_manifest.sh FIXTURE_DIR [MANIFEST]
# Recomputes byte size and SHA-256 of every file in manifest.tsv (paths relative to FIXTURE_DIR).
# MANIFEST defaults to FIXTURE_DIR/manifest.tsv. Exits 1 on any missing file or mismatch. No BASH_SOURCE.
set -uo pipefail
dir=${1:?usage: check_manifest.sh ABSOLUTE_FIXTURE_DIR [MANIFEST]}
man=${2:-$dir/manifest.tsv}
case "$dir" in /*) ;; *) echo "error: fixture dir must be an absolute path" >&2; exit 2;; esac
[ -f "$man" ] || { echo "error: no manifest $man" >&2; exit 2; }
bad=0; n=0
while IFS=$'\t' read -r file bytes sha how; do
    [ "$file" = "file" ] && continue
    [ -z "$file" ] && continue
    n=$((n+1))
    if [ ! -f "$dir/$file" ]; then echo "MISSING $file"; bad=$((bad+1)); continue; fi
    b=$(stat -c %s "$dir/$file"); s=$(sha256sum "$dir/$file" | cut -d' ' -f1)
    if [ "$b" != "$bytes" ] || [ "$s" != "$sha" ]; then
        echo "MISMATCH $file (bytes $b vs $bytes; sha256 $s vs $sha)"; bad=$((bad+1))
    fi
done < "$man"
echo "checked $n files, $bad problem(s)"
[ "$bad" -eq 0 ]
