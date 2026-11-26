#!/usr/bin/env bash
# Build the pre-built index set on a Slurm cluster: one array task per genome of
# crates/genome-catalog/data/indexes-manifest-targets.json, each running scripts/build-one-index.sh (kallisto.wasm index
# plus its self-check). The same work as .github/workflows/build-indexes.yml, for the HPCC.
#
#   OUT=bench/<id>_indexes-v2/out sbatch --array=0-49%10 -p short -c 2 --mem 8G -t 1:00:00 \
#     -o bench/<id>_indexes-v2/logs/%a.log scripts/build-indexes-slurm.sh
# Then upload the set:  SRC_DIR=$OUT scripts/mirror-indexes-to-s3.sh indexes-v2
# Needs: target/release/qrnamap (cargo build --release -p qrnamap-cli), node, curl, gzip, sha256sum.
# Work files go to $SCRATCH (build-one-index.sh); OUT keeps only the four .gz files and the .sha256 of each genome.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
: "${OUT:?set OUT to the output folder}"
i="${SLURM_ARRAY_TASK_ID:?run as a Slurm array (or set SLURM_ARRAY_TASK_ID)}"
export PATH="$PWD/target/release:$PATH"
module load nodejs 2>/dev/null || true
command -v qrnamap >/dev/null || { echo "qrnamap not built: cargo build --release -p qrnamap-cli" >&2; exit 1; }

read -r acc fna gff < <(python3 - "$i" <<'EOF'
import json, sys
t = json.load(open('crates/genome-catalog/data/indexes-manifest-targets.json'))
i = int(sys.argv[1])
if i < len(t):
    g = t[i]
    print(g['accession'], g['fasta']['url'], g['gff']['url'])
EOF
)
[[ -n "${acc:-}" ]] || { echo "no genome at index $i"; exit 0; }
mkdir -p "$OUT"
start=$(date +%s)
bash scripts/build-one-index.sh "$acc" "$fna" "$gff" "$OUT"
echo "$acc done in $(( $(date +%s) - start )) s"
