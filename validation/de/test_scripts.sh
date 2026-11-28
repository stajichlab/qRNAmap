#!/usr/bin/env bash
# Fixture test of page_de.mjs and compare_de.py on a small synthetic matrix (no R needed).
#   bash validation/de/test_scripts.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
node -e '
let s = 3; const r = () => ((s = (s * 1103515245 + 12345) % 2147483648) / 2147483648);
const n = ["c1","c2","c3","t1","t2","t3"]; const L = ["gene\t" + n.join("\t")];
for (let g = 0; g < 300; g++) { const mu = 50 + 2000 * r();
  L.push("g" + g + "\t" + n.map((x) => Math.round(mu * (x[0] === "t" && g < 15 ? 4 : 1) * Math.exp(0.15 * (r() + r() + r() - 1.5)))).join("\t")); }
require("fs").writeFileSync(process.argv[1] + "/counts.tsv", L.join("\n") + "\n");
require("fs").writeFileSync(process.argv[1] + "/coldata.tsv", "sample\tcondition\n" + n.map((x) => x + "\t" + (x[0] === "t" ? "trt" : "ctl")).join("\n") + "\n");
' "$T"
node "$HERE/page_de.mjs" "$T/counts.tsv" "$T/coldata.tsv" condition trt ctl "$T/page.tsv" | grep -q "300 genes" || { echo "FAIL: page_de"; exit 1; }
[[ $(head -1 "$T/page.tsv") == $'gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj' ]] || { echo "FAIL: header"; exit 1; }
python3 "$HERE/compare_de.py" "$T/page.tsv" "$T/page.tsv" > "$T/cmp"
grep -q "reference only 0, other only 0" "$T/cmp" && grep -q "log2FoldChange: Pearson 1.0000, Spearman 1.0000" "$T/cmp" || { cat "$T/cmp"; echo "FAIL: self-compare"; exit 1; }
grep -Eq "padj < 0.05: both (1[5-9]|2[0-9])," "$T/cmp" || { cat "$T/cmp"; echo "FAIL: the 15 changed genes"; exit 1; }
echo "ok: DE comparison scripts"
