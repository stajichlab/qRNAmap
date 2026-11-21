#!/usr/bin/env bash
# Browser spike 0: copy the spike workflow and scripts into the genomics-wasm fork and push branch qrnamap-spike0.
# No BASH_SOURCE: every path comes from the arguments (absolute paths).
#
# Usage: sync.sh SRC_DIR FORK_DIR "commit message" [--push]
#   SRC_DIR   absolute path of validation/browser-spike0 in the qRNAmap repository
#   FORK_DIR  absolute path of a clone of https://github.com/stajichlab/genomics-wasm
#   --push    push to origin (must be stajichlab/genomics-wasm). Without it, the script only commits.
#
# Copies:  SRC_DIR/workflow/qrnamap-spike0.yml -> FORK_DIR/.github/workflows/qrnamap-spike0.yml
#          SRC_DIR/{*.mjs,*.py,deseq2_ref.R,manifest.tsv} -> FORK_DIR/spike0/
# The script refuses to run when the fork branch has changes outside these two paths relative to the pin.
set -euo pipefail
PIN=0b015f63e8048eb57b72ffd2706b79b803862b7c
BRANCH=qrnamap-spike0
SRC=${1:?SRC_DIR}; FORK=${2:?FORK_DIR}; MSG=${3:?commit message}; PUSH=${4:-}
case "$SRC$FORK" in *..*) echo "no .. in paths" >&2; exit 2;; esac
[[ "$SRC" = /* && "$FORK" = /* ]] || { echo "SRC_DIR and FORK_DIR must be absolute" >&2; exit 2; }
[ -f "$SRC/workflow/qrnamap-spike0.yml" ] || { echo "missing $SRC/workflow/qrnamap-spike0.yml" >&2; exit 2; }
url=$(git -C "$FORK" remote get-url origin)
case "$url" in
  https://github.com/stajichlab/genomics-wasm.git|https://github.com/stajichlab/genomics-wasm|git@github.com:stajichlab/genomics-wasm.git) ;;
  *) echo "origin of $FORK is $url, not stajichlab/genomics-wasm" >&2; exit 2;;
esac
git -C "$FORK" fetch -q origin "$BRANCH"
cur=$(git -C "$FORK" rev-parse --abbrev-ref HEAD)
if [ "$cur" != "$BRANCH" ]; then git -C "$FORK" checkout -q -B "$BRANCH" "origin/$BRANCH"; fi
git -C "$FORK" merge-base --is-ancestor "$PIN" HEAD || { echo "pin $PIN is not an ancestor of $BRANCH" >&2; exit 2; }

mkdir -p "$FORK/.github/workflows" "$FORK/spike0"
cp "$SRC/workflow/qrnamap-spike0.yml" "$FORK/.github/workflows/qrnamap-spike0.yml"
shopt -s nullglob
for f in "$SRC"/*.mjs "$SRC"/*.py "$SRC/deseq2_ref.R" "$SRC/manifest.tsv"; do cp "$f" "$FORK/spike0/"; done
git -C "$FORK" add .github/workflows/qrnamap-spike0.yml spike0
if git -C "$FORK" diff --cached --quiet; then
  echo "no changes to commit"
else
  git -C "$FORK" commit -q -m "$MSG"
fi
# Only the spike workflow and spike0/ may differ from the pinned upstream tree.
other=$(git -C "$FORK" diff --name-only "$PIN" HEAD | grep -v -e '^\.github/workflows/qrnamap-spike0\.yml$' -e '^spike0/' || true)
[ -z "$other" ] || { echo "unexpected changes outside the spike paths:" >&2; echo "$other" >&2; exit 3; }
git -C "$FORK" diff --stat "$PIN" HEAD
if [ "$PUSH" = "--push" ]; then
  git -C "$FORK" push origin "HEAD:refs/heads/$BRANCH"
fi
echo "fork head: $(git -C "$FORK" rev-parse HEAD)"
