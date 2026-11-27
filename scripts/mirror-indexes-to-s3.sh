#!/usr/bin/env bash
# Mirror the pre-built kallisto indexes of a GitHub release (made by .github/workflows/build-indexes.yml) to the lab S3
# bucket, where the web page and the desktop app can fetch them. GitHub release downloads cannot be fetched from a web
# page: github.com answers with a redirect that has no Access-Control-Allow-Origin header. The S3 bucket has a CORS rule.
#
# Usage: [SRC_DIR=dir] mirror-indexes-to-s3.sh [TAG] [S3_PREFIX]
#   SRC_DIR    take the files from this folder instead of the GitHub release (its <acc>.sha256 files are checked)
#   TAG        release tag (default indexes-v1); a draft release works (gh reads it with your login)
#   S3_PREFIX  destination (default s3://stajichlab/qRNAmap/<TAG>/)
# Needs: gh (logged in), s3cmd (~/.s3cfg), sha256sum. Run it on the HPCC as a job: the files are staged in $SCRATCH.
#   sbatch -p short -c 2 --mem 4G -t 1:00:00 -o bench/<id>_s3-indexes/mirror.%j.log scripts/mirror-indexes-to-s3.sh
#
# Objects are public-read. They stay gzip files as built (Content-Type application/gzip, no Content-Encoding): the page
# gunzips them itself, so a browser must not decompress them on the way.
# Check without uploading: DRY_RUN=1 (downloads and checks, prints the s3cmd commands).
set -euo pipefail

TAG="${1:-indexes-v1}"
DEST="${2:-s3://stajichlab/qRNAmap/${TAG}/}"
REPO="${REPO:-stajichlab/qRNAmap}"
[[ "$DEST" == */ ]] || DEST="$DEST/"

# SRC_DIR as an absolute path: the script works inside $WORK.
if [[ -n "${SRC_DIR:-}" ]]; then SRC_DIR="$(cd "$SRC_DIR" && pwd)"; fi
WORK="${SCRATCH:-${TMPDIR:-/tmp}}/qrnamap_mirror_${TAG}_$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"

if [[ -n "${SRC_DIR:-}" ]]; then
  # A local set (for example made by an HPCC array of build-one-index.sh) instead of a GitHub release.
  echo "== copy local set $SRC_DIR into $WORK"
  cp "$SRC_DIR"/*.gz "$WORK"/
  cat "$SRC_DIR"/*.sha256 2>/dev/null | sort -u -k2 > "$WORK/checksums.sha256" || true
  [[ -s "$WORK/checksums.sha256" ]] || rm -f "$WORK/checksums.sha256"
else
  echo "== download release $TAG of $REPO into $WORK"
  start=$(date +%s)
  gh release download "$TAG" --repo "$REPO" --dir "$WORK" --clobber
  echo "downloaded $(ls | wc -l) files, $(du -sh . | cut -f1), in $(( $(date +%s) - start )) s"
fi

echo "== check the files"
# checksums.sha256 of indexes-v1 lists only one genome (the build artifacts overwrote each other; fixed in
# build-indexes.yml for later tags). So: check the files it lists, test every gzip (CRC and length of each member),
# then write a complete list from the downloaded files and upload that.
if [[ -f checksums.sha256 ]]; then
  sha256sum --check --quiet checksums.sha256
  echo "$(wc -l < checksums.sha256) files listed in checksums.sha256 match"
  mv checksums.sha256 checksums.release.sha256
fi
n=0
for f in ./*.gz; do gzip -t "$f"; n=$((n + 1)); done
echo "gzip -t passed for all $n files"
sha256sum ./*.gz | sed 's#  \./#  #' > checksums.sha256
if [[ -f checksums.release.sha256 ]]; then
  # every line of the release list must be in the new list
  sort checksums.release.sha256 | comm -23 - <(sort checksums.sha256) | { ! grep . ; } || { echo "release checksums disagree with the files" >&2; exit 1; }
fi

run() { if [[ "${DRY_RUN:-0}" == 1 ]]; then echo "DRY_RUN: $*"; else "$@"; fi; }

echo "== upload to $DEST"
start=$(date +%s)
# Cache for a day: the files of one tag do not change (a rebuild gets a new tag).
run s3cmd put --acl-public --no-preserve --mime-type=application/gzip \
  --add-header="Cache-Control:public, max-age=86400" \
  ./*.gz "$DEST"
run s3cmd put --acl-public --no-preserve --mime-type=text/plain \
  --add-header="Cache-Control:public, max-age=86400" \
  checksums.sha256 "$DEST"
echo "uploaded in $(( $(date +%s) - start )) s"

echo "== verify"
if [[ "${DRY_RUN:-0}" != 1 ]]; then
  remote=$(s3cmd ls "$DEST" | grep -c '\.gz$' || true)
  local_n=$(ls ./*.gz | wc -l)
  echo "remote .gz objects: $remote, local: $local_n"
  [[ "$remote" -ge "$local_n" ]] || { echo "upload incomplete" >&2; exit 1; }
fi
echo "done"
