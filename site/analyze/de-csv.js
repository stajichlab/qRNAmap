// The results CSV of DESeq2 (`qrnamap de`, and the desktop app's "Run DESeq2 (R)"), read for the page's volcano and DE table.

/**
 * Read the results CSV that DESeq2 writes (`write.csv`: a quoted header, quoted gene names, NA for missing) into the rows
 * the volcano and the DE table use: {gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj}, with NaN for NA (a NaN
 * padj never passes a cutoff, as DESeq2's NA means "not tested").
 */
export function parseDeseq2Csv(text) {
  const splitCsv = (line) => {
    const out = [];
    let cur = '';
    let quoted = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (quoted) {
        if (c === '"' && line[i + 1] === '"') { cur += '"'; i++; } else if (c === '"') quoted = false; else cur += c;
      } else if (c === '"') quoted = true;
      else if (c === ',') { out.push(cur); cur = ''; } else cur += c;
    }
    out.push(cur);
    return out;
  };
  const lines = String(text || '').replace(/\r/g, '').split('\n').filter((l) => l.length);
  if (!lines.length) return [];
  const head = splitCsv(lines[0]);
  const at = (name) => head.indexOf(name);
  const need = ['gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'];
  const missing = need.filter((n) => at(n) < 0);
  if (missing.length) throw new Error(`The DESeq2 results have no column ${missing.join(', ')}.`);
  const num = (v) => (v === 'NA' || v === '' || v === undefined ? NaN : Number(v));
  return lines.slice(1).map((line) => {
    const c = splitCsv(line);
    return {
      gene: c[at('gene')],
      baseMean: num(c[at('baseMean')]),
      log2FoldChange: num(c[at('log2FoldChange')]),
      lfcSE: num(c[at('lfcSE')]),
      stat: num(c[at('stat')]),
      pvalue: num(c[at('pvalue')]),
      padj: num(c[at('padj')]),
    };
  });
}
