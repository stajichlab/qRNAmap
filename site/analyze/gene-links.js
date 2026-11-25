// External database links for genes on the results page. Built from the Dbxref entries of the
// GFF gene feature (the `dbxref` column of gene_info.tsv) and from the gene label itself.
// Pure functions: no DOM, so units.test.mjs can test them in Node.

/** FungiDB gene permalink; the record id is the locus tag (for example AFUA_1G00100, AFUB_102960). */
const fungidbUrl = (id) => `https://fungidb.org/gene/${id}`;

/**
 * Links for one gene, most specific database first.
 *
 * `gene` is the gene label (locus tag, the `gene` column of gene_info.tsv); `info` is its parsed
 * row { symbol, geneId, dbxref }. Returns [] for genes with nothing to link.
 */
export function geneExternalLinks(gene, info) {
  if (!info) return [];
  const links = [];
  const seen = new Set();
  const add = (name, url) => {
    if (url && !seen.has(url)) {
      seen.add(url);
      links.push({ name, url });
    }
  };
  const dbxref = info.dbxref || [];
  const sgd = dbxref.find((x) => x.startsWith('SGD:'));
  if (sgd && sgd.length > 4) {
    add('SGD', `https://www.yeastgenome.org/locus/${sgd.slice(4)}`);
  }
  const fungidb = dbxref.find((x) => x.startsWith('FungiDB:'));
  if (fungidb && fungidb.length > 8) {
    add('FungiDB', fungidbUrl(fungidb.slice(8)));
  }
  // Aspergillus fumigatus locus tags (Af293 AFUA_1G00100, A1163 AFUB_102960) are FungiDB record
  // ids even when the GFF carries no Dbxref; the A1163 GenBank annotation has none at all.
  if (/^AFU[A-Z]_\d/.test(gene || '')) {
    add('FungiDB', fungidbUrl(gene));
  }
  if (/^\d+$/.test(info.geneId || '')) {
    add('NCBI Gene', `https://www.ncbi.nlm.nih.gov/gene/${info.geneId}`);
  }
  return links;
}
