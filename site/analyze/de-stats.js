// Quick two-group differential expression for the page (a screen, not DESeq2).
//
// Counts are scaled by DESeq2-style median-of-ratios size factors, then each sample's value is
// log2(normalized + 0.5). Each gene gets a two-sample t-test on that log scale with a moderated
// (empirical-Bayes-shrunk) pooled variance, as in limma: s2_post = (d0*s0^2 + d*s2) / (d0 + d),
// with the prior s0^2 the median gene variance among expressed genes and a fixed prior weight
// d0 = PRIOR_DF. P-values come from Student's t with d + d0 degrees of freedom, then
// Benjamini-Hochberg. The fold change and its standard error are on the same (log2) scale, so the
// test does not lose power as read depth grows. For publication use `qrnamap de` (DESeq2).

/** Prior degrees of freedom of the moderated variance (limma's d0 is typically 2-10 for RNA-seq). */
export const PRIOR_DF = 4;
/** Genes with a mean normalized count below this do not set the prior variance. */
const PRIOR_MIN_MEAN = 10;

/** ln Γ(x) (Lanczos, g = 7). */
function lnGamma(x) {
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lnGamma(1 - x);
  x -= 1;
  let a = c[0];
  const t = x + 7.5;
  for (let i = 1; i < 9; i++) a += c[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

/** Continued fraction of the regularized incomplete beta (Numerical Recipes betacf). */
function betaCf(a, b, x) {
  const FPMIN = 1e-300;
  let c = 1;
  let d = 1 - ((a + b) * x) / (a + 1);
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 300; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((a + m2 - 1) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = (-(a + m) * (a + b + m) * x) / ((a + m2) * (a + m2 + 1));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-14) break;
  }
  return h;
}

/** Regularized incomplete beta I_x(a, b). */
export function incompleteBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(lnGamma(a + b) - lnGamma(a) - lnGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? (bt * betaCf(a, b, x)) / a : 1 - (bt * betaCf(b, a, 1 - x)) / b;
}

/** Two-sided p-value of Student's t with df degrees of freedom. */
export function tTwoSided(t, df) {
  if (!Number.isFinite(t)) return Number.isNaN(t) ? 1 : 0;
  return Math.min(1, incompleteBeta(df / (df + t * t), df / 2, 0.5));
}

/** Benjamini-Hochberg adjusted p-values, in input order. */
export function bhAdjust(pvalues) {
  const m = pvalues.length;
  const order = pvalues.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]);
  const out = new Array(m);
  let cumMin = 1;
  for (let k = m - 1; k >= 0; k--) {
    cumMin = Math.min(cumMin, (order[k][0] * m) / (k + 1), 1);
    out[order[k][1]] = cumMin;
  }
  return out;
}

/** DESeq2 median-of-ratios size factors (column totals when fewer than 10 genes are non-zero everywhere). */
export function sizeFactors(rows, nSamples) {
  const logGeo = rows.map((row) => (row.every((v) => v > 0) ? row.reduce((s, v) => s + Math.log(v), 0) / nSamples : null));
  const usable = logGeo.filter((g) => g !== null).length;
  const sf = new Array(nSamples).fill(1);
  if (usable >= 10) {
    for (let j = 0; j < nSamples; j++) {
      const ratios = [];
      rows.forEach((row, i) => { if (logGeo[i] !== null) ratios.push(row[j] / Math.exp(logGeo[i])); });
      ratios.sort((a, b) => a - b);
      const mid = Math.floor(ratios.length / 2);
      sf[j] = ratios.length % 2 ? ratios[mid] : (ratios[mid - 1] + ratios[mid]) / 2;
    }
  } else {
    const totals = new Array(nSamples).fill(0);
    for (const row of rows) row.forEach((v, j) => { totals[j] += v; });
    const mean = totals.reduce((a, b) => a + b, 0) / nSamples || 1;
    for (let j = 0; j < nSamples; j++) sf[j] = (totals[j] || 1) / mean;
  }
  return sf;
}

function median(xs) {
  if (!xs.length) return NaN;
  const s = [...xs].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
}

/**
 * Two-group differential expression of a gene-count matrix (TSV, first column gene, one column per sample).
 * `samples`: [{name, covariates}], `contrast`: {factor, level, reference}. Returns one row per gene:
 * {gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj}. With no replicate in either group there is
 * no variance to test against: pvalue and padj are 1.
 */
export function computeDifferentialExpression(matrixTsv, samples, contrast) {
  const lines = matrixTsv.trim().split('\n');
  if (lines.length <= 1) return [];
  const sampleCols = lines[0].split('\t').slice(1);
  const factor = contrast?.factor || 'condition';
  const groupOf = new Map(samples.map((s) => [s.name, String(s.covariates?.[factor] ?? '').trim()]));
  const trt = [];
  const ref = [];
  sampleCols.forEach((name, j) => {
    const g = groupOf.get(name);
    if (g === contrast?.level) trt.push(j);
    else if (g === contrast?.reference) ref.push(j);
  });
  if (!trt.length || !ref.length) return [];

  const genes = [];
  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split('\t');
    genes.push(parts[0]);
    rows.push(parts.slice(1).map(Number));
  }
  const sf = sizeFactors(rows, sampleCols.length);
  const nT = trt.length;
  const nR = ref.length;
  const d = nT + nR - 2;

  const genesStats = rows.map((row) => {
    const norm = row.map((v, j) => v / (sf[j] || 1));
    const lt = trt.map((j) => Math.log2(norm[j] + 0.5));
    const lr = ref.map((j) => Math.log2(norm[j] + 0.5));
    const mt = lt.reduce((a, b) => a + b, 0) / nT;
    const mr = lr.reduce((a, b) => a + b, 0) / nR;
    const ss = lt.reduce((a, v) => a + (v - mt) ** 2, 0) + lr.reduce((a, v) => a + (v - mr) ** 2, 0);
    const baseMean = [...trt, ...ref].reduce((a, j) => a + norm[j], 0) / (nT + nR);
    return { baseMean, lfc: mt - mr, s2: d > 0 ? ss / d : NaN };
  });

  let s0sq = median(genesStats.filter((g) => g.baseMean >= PRIOR_MIN_MEAN && g.s2 > 0).map((g) => g.s2));
  if (!(s0sq > 0)) s0sq = median(genesStats.filter((g) => g.s2 > 0).map((g) => g.s2));
  const d0 = s0sq > 0 ? PRIOR_DF : 0;
  const dfTotal = d + d0;

  const results = genesStats.map((g, i) => {
    if (d === 0 || dfTotal === 0 || g.baseMean === 0) {
      return { gene: genes[i], baseMean: g.baseMean, log2FoldChange: g.lfc, lfcSE: NaN, stat: 0, pvalue: 1, padj: 1 };
    }
    const s2post = (d0 * (s0sq || 0) + d * g.s2) / dfTotal;
    const se = Math.sqrt(s2post * (1 / nT + 1 / nR));
    const stat = se > 0 ? g.lfc / se : 0;
    return { gene: genes[i], baseMean: g.baseMean, log2FoldChange: g.lfc, lfcSE: se, stat, pvalue: tTwoSided(stat, dfTotal), padj: 1 };
  });
  const padj = bhAdjust(results.map((r) => r.pvalue));
  results.forEach((r, i) => { r.padj = padj[i]; });
  return results;
}
