#!/usr/bin/env bash
# Check a vendored directory against the SHA-256 list in its PROVENANCE.md.
# Usage: check-vendor.sh <absolute directory>
# Exit 0: all files match. Exit 1: mismatch, missing file, extra file, symlink,
# bad or duplicate row, or missing fence. Exit 2: usage error.
# The check trusts PROVENANCE.md itself. No BASH_SOURCE is used.
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

status=0

# The rows are inside the last fenced block (a line of three backticks opens and closes it).
fences="$(grep -c '^```$' "$prov" || true)"
if [ "$fences" -lt 2 ]; then
  echo "BAD: PROVENANCE.md needs an opening and a closing fence line" >&2
  exit 1
fi
block="$(awk '/^```$/ { n++; next } n % 2 == 1 { print }' "$prov")"

rows=""
while IFS= read -r line; do
  [ -z "$line" ] && continue
  if ! [[ "$line" =~ ^[0-9a-f]{64}\ \ [^\ ]+$ ]]; then
    echo "BAD ROW: $line" >&2
    status=1
    continue
  fi
  path="${line:66}"
  case "$path" in
    /*|./*|../*|*/../*|*/..|..)
      echo "BAD PATH: $path" >&2
      status=1
      continue ;;
  esac
  rows+="$line"$'\n'
done <<< "$block"

if [ -z "$rows" ]; then
  echo "BAD: no SHA-256 rows in $prov" >&2
  exit 1
fi
rows="$(printf '%s' "$rows" | LC_ALL=C sort -k2)"

dups="$(printf '%s\n' "$rows" | cut -c67- | LC_ALL=C sort | uniq -d)"
if [ -n "$dups" ]; then
  while IFS= read -r p; do echo "DUPLICATE: $p" >&2; done <<< "$dups"
  status=1
fi

while read -r hash path; do
  if [ -L "$dir/$path" ]; then
    echo "SYMLINK: $path" >&2
    status=1
  elif [ ! -f "$dir/$path" ]; then
    echo "MISSING: $path" >&2
    status=1
  else
    actual="$(sha256sum "$dir/$path" | cut -d' ' -f1)"
    if [ "$actual" != "$hash" ]; then
      echo "MISMATCH: $path (expected $hash, got $actual)" >&2
      status=1
    fi
  fi
done <<< "$rows"

# Extra entries: every file or symlink except PROVENANCE.md must be listed.
listed="$(printf '%s\n' "$rows" | cut -c67- | LC_ALL=C sort -u)"
found="$(cd "$dir" && find . ! -path ./PROVENANCE.md \( -type f -o -type l \) | sed 's|^\./||' | LC_ALL=C sort)"
extra="$(LC_ALL=C comm -13 <(printf '%s\n' "$listed") <(printf '%s\n' "$found"))"
if [ -n "$extra" ]; then
  while IFS= read -r p; do echo "EXTRA: $p" >&2; done <<< "$extra"
  status=1
fi

if [ "$status" -eq 0 ]; then
  echo "ok: $(printf '%s\n' "$rows" | wc -l) files match $prov"
fi
exit "$status"
