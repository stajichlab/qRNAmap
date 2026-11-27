// The page's DE (site/analyze/de-stats.js) on a count matrix, written like deseq2_reference.R.
//   node page_de.mjs COUNTS.tsv COLDATA.tsv FACTOR LEVEL REFERENCE OUT.tsv
import { readFileSync, writeFileSync } from 'node:fs';
import { computeDifferentialExpression } from '../../site/analyze/de-stats.js';

const [counts, coldata, factor, level, reference, out] = process.argv.slice(2);
if (!out) { console.error('usage: node page_de.mjs COUNTS COLDATA FACTOR LEVEL REFERENCE OUT'); process.exit(2); }
const [head, ...rows] = readFileSync(coldata, 'utf8').trim().split('\n').map((l) => l.split('\t'));
const samples = rows.map((r) => ({ name: r[head.indexOf('sample')], covariates: Object.fromEntries(head.map((h, i) => [h, r[i]])) }));
const res = computeDifferentialExpression(readFileSync(counts, 'utf8'), samples, { factor, level, reference });
const fmt = (x) => (Number.isFinite(x) ? String(x) : 'NA');
writeFileSync(out, ['gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj',
  ...res.map((r) => [r.gene, r.baseMean, r.log2FoldChange, r.lfcSE, r.stat, r.pvalue, r.padj].map((v, i) => (i ? fmt(v) : v)).join('\t'))].join('\n') + '\n');
console.log(`page DE: ${res.length} genes, padj < 0.05: ${res.filter((r) => r.padj < 0.05).length}`);
