#!/usr/bin/env bash
# Fixture test of mirror-indexes-to-s3.sh with stub `gh` and `s3cmd` on PATH (no network, no bucket).
#   bash scripts/test-mirror-indexes.sh
# Cases: a good release uploads every .gz plus a complete checksums.sha256; a damaged gzip stops the run before any
# upload; a release checksum that disagrees with its file stops the run.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/bin" "$T/release"
SHA="sha256sum"; command -v sha256sum >/dev/null || SHA="shasum -a 256"

make_release() { # $1: mode (good | badgz | badsum)
  rm -rf "$T/release"/*; : > "$T/s3.log"
  for acc in GCF_1 GCF_2; do
    for kind in idx tx2gene.tsv gene_info.tsv genes.tsv; do echo "$acc $kind" | gzip -c > "$T/release/$acc.$kind.gz"; done
  done
  (cd "$T/release" && $SHA GCF_2.idx.gz GCF_2.genes.tsv.gz > checksums.sha256) # only one genome, as in indexes-v1
  case "$1" in
    badgz) head -c 10 "$T/release/GCF_1.idx.gz" > "$T/x" && mv "$T/x" "$T/release/GCF_1.idx.gz" ;;
    badsum) echo "different" | gzip -c > "$T/release/GCF_2.idx.gz" ;;
  esac
}

cat > "$T/bin/gh" <<EOF
#!/usr/bin/env bash
# gh release download TAG --repo R --dir DIR --clobber
dir=""; while [[ \$# -gt 0 ]]; do [[ "\$1" == --dir ]] && dir="\$2"; shift; done
cp "$T/release"/* "\$dir"/
EOF
cat > "$T/bin/s3cmd" <<EOF
#!/usr/bin/env bash
echo "\$*" >> "$T/s3.log"
if [[ "\$1" == ls ]]; then for f in $T/uploaded/*; do echo "2026-11-26 00:00 1 s3://b/\$(basename "\$f")"; done; exit 0; fi
if [[ "\$1" == put ]]; then mkdir -p "$T/uploaded"; for a in "\$@"; do [[ -f "\$a" ]] && cp "\$a" "$T/uploaded/"; done; fi
exit 0
EOF
if ! command -v sha256sum >/dev/null; then # macOS: a sha256sum that understands --check --quiet
  cat > "$T/bin/sha256sum" <<'EOF'
#!/usr/bin/env bash
args=(); for a in "$@"; do [[ "$a" == --quiet ]] || args+=("$a"); done
if [[ "${args[0]:-}" == --check ]]; then shasum -a 256 -c "${args[@]:1}" >/dev/null; else shasum -a 256 "${args[@]}"; fi
EOF
fi
chmod +x "$T/bin"/*
export PATH="$T/bin:$PATH" SCRATCH="$T/scratch"
mkdir -p "$SCRATCH"

fail() { echo "FAIL: $*" >&2; exit 1; }

make_release good; rm -rf "$T/uploaded"
bash "$HERE/mirror-indexes-to-s3.sh" indexes-test s3://b/ > "$T/out" 2>&1 || { cat "$T/out"; fail "good release"; }
[[ $(ls "$T/uploaded"/*.gz | wc -l) -eq 8 ]] || fail "8 .gz uploaded"
[[ $(wc -l < "$T/uploaded/checksums.sha256") -eq 8 ]] || fail "complete checksums.sha256 (8 lines)"
grep -q -- "--acl-public" "$T/s3.log" || fail "public-read"
grep -q -- "--mime-type=application/gzip" "$T/s3.log" || fail "gzip content type"
echo "ok: good release"

make_release badgz; rm -rf "$T/uploaded"
if bash "$HERE/mirror-indexes-to-s3.sh" indexes-test s3://b/ > "$T/out" 2>&1; then fail "damaged gzip accepted"; fi
grep -q put "$T/s3.log" && fail "upload after a damaged gzip"
echo "ok: damaged gzip stops before upload"

make_release badsum; rm -rf "$T/uploaded"
if bash "$HERE/mirror-indexes-to-s3.sh" indexes-test s3://b/ > "$T/out" 2>&1; then fail "bad checksum accepted"; fi
grep -q put "$T/s3.log" && fail "upload after a bad checksum"
echo "ok: release checksum mismatch stops before upload"

# A local set given as a relative SRC_DIR (as from the HPCC array): uploads every .gz and the joined .sha256 list.
make_release good; rm -rf "$T/uploaded"; mkdir -p "$T/local"; cp "$T/release"/*.gz "$T/local/"
(cd "$T/local" && $SHA GCF_1.idx.gz > GCF_1.sha256 && $SHA GCF_2.idx.gz > GCF_2.sha256)
(cd "$T" && SRC_DIR=local bash "$HERE/mirror-indexes-to-s3.sh" indexes-test s3://b/ > "$T/out" 2>&1) || { cat "$T/out"; fail "local relative SRC_DIR"; }
[[ $(ls "$T/uploaded"/*.gz | wc -l) -eq 8 ]] || fail "local set: 8 .gz uploaded"
echo "ok: local set with a relative SRC_DIR"
echo "all mirror tests passed"
