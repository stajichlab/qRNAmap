#!/usr/bin/env bash
# Check a vendored directory against the SHA-256 list in its PROVENANCE.md.
# Usage: check-vendor.sh <absolute directory>
# Exit 0: all files match. Exit 1: mismatch, missing file or extra file. Exit 2: usage error.
# No BASH_SOURCE is used.
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: check-vendor.sh <absolute directory>" >&2
  exit 2
fi
dir="$1"
case "$dir" in
  /*) ;;
  *) echo "error: the directory must be an absolute path: $dir" >&2; exit 2 ;;
esac
if [ ! -d "$dir" ]; then
  echo "error: not a directory: $dir" >&2
  exit 2
fi
prov="$dir/PROVENANCE.md"
if [ ! -f "$prov" ]; then
  echo "error: missing $prov" >&2
  exit 2
fi

# Expected: lines of the form "<64 hex>  <relative path>" in PROVENANCE.md.
expected="$(grep -E '^[0-9a-f]{64}  [^ ]+$' "$prov" | LC_ALL=C sort -k2)"
if [ -z "$expected" ]; then
  echo "error: no SHA-256 lines in $prov" >&2
  exit 2
fi

status=0
while read -r hash path; do
  if [ ! -f "$dir/$path" ]; then
    echo "MISSING: $path" >&2
    status=1
    continue
  fi
  actual="$(sha256sum "$dir/$path" | cut -d' ' -f1)"
  if [ "$actual" != "$hash" ]; then
    echo "MISMATCH: $path (expected $hash, got $actual)" >&2
    status=1
  fi
done <<< "$expected"

# Extra files: every file except PROVENANCE.md must be listed.
listed="$(printf '%s\n' "$expected" | cut -d' ' -f3- | sed 's/^ *//' | LC_ALL=C sort)"
found="$(cd "$dir" && find . -type f ! -path ./PROVENANCE.md | sed 's|^\./||' | LC_ALL=C sort)"
extra="$(LC_ALL=C comm -13 <(printf '%s\n' "$listed") <(printf '%s\n' "$found"))"
if [ -n "$extra" ]; then
  while read -r p; do echo "EXTRA: $p" >&2; done <<< "$extra"
  status=1
fi

if [ "$status" -eq 0 ]; then
  echo "ok: $(printf '%s\n' "$expected" | wc -l) files match $prov"
fi
exit "$status"
