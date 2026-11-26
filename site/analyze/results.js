// Step 3: results view, downloads and plots.
import { assemble_matrix, qc, de_command, validate_sheet } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { readStore, writeStore } from '../persist.js?v=__BUILD__';
import { geneExternalLinks } from './gene-links.js?v=__BUILD__';

const STORE_KEY = 'results';

/** Escape text for use in SVG markup (text content and double-quoted attributes). */
const escSvg = (s) =>
  String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

/** Parse gene_info.tsv into a map of gene -> { tokens, symbol, geneId, dbxref }. */
export function parseAliasIndex(geneInfoTsv) {
  const map = new Map();
  if (!geneInfoTsv) return map;
  const lines = geneInfoTsv.trim().split('\n');
  if (lines.length <= 1) return map;
  // gene_info.tsv may come from an older build without the dbxref column.
  const dbxrefCol = lines[0].split('\t').indexOf('dbxref');
  for (const line of lines.slice(1)) {
    const cols = line.split('\t');
    if (!cols.length) continue;
    const gene = cols[0];
    // gene, symbol, gene_id, old_locus_tag, transcripts, proteins
    const tokens = cols.slice(0, 6).filter(Boolean).join(' ').toLowerCase();
    const dbxref = dbxrefCol >= 0 && cols[dbxrefCol] ? cols[dbxrefCol].split(',').filter(Boolean) : [];
    map.set(gene, { tokens, symbol: cols[1] || '', geneId: cols[2] || '', dbxref });
  }
  return map;
}

/** Format a number with commas. */
const fmt = (n) => Number(n).toLocaleString('en-US');

/** Helper to trigger file download in browser. */
function downloadFile(name, text) {
  downloadBlob(text, name, 'text/tab-separated-values;charset=utf-8');
}

/** Download `text` as `name` with the given MIME type. */
function downloadBlob(text, name, mime) {
  const blob = new Blob([text], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = name;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** Render SVG correlation heatmap. */
/** RGB color interpolation between multi-stop color ramps. */
function interpolateRgb(stops, t) {
  const clamped = Math.max(0, Math.min(1, t));
  for (let i = 0; i < stops.length - 1; i++) {
    const [p0, c0] = stops[i];
    const [p1, c1] = stops[i + 1];
    if (clamped >= p0 && clamped <= p1) {
      const f = p1 === p0 ? 0 : (clamped - p0) / (p1 - p0);
      return [
        Math.round(c0[0] + f * (c1[0] - c0[0])),
        Math.round(c0[1] + f * (c1[1] - c0[1])),
        Math.round(c0[2] + f * (c1[2] - c0[2])),
      ];
    }
  }
  return stops[stops.length - 1][1];
}

/** Perceived relative luminance of an [R, G, B] tuple (0..1). */
function rgbLuminance(rgb) {
  const a = rgb.map((v) => {
    const s = v / 255;
    return s <= 0.03928 ? s / 12.92 : Math.pow((s + 0.055) / 1.055, 2.4);
  });
  return a[0] * 0.2126 + a[1] * 0.7152 + a[2] * 0.0722;
}

/** Sequential YlGnBu palette for non-negative correlations. */
const HEATMAP_YL_GN_BU = [
  [0.00, [255, 255, 217]], // Pale cream (#ffffd9)
  [0.25, [199, 233, 180]], // Soft green (#c7e9b4)
  [0.50, [65, 182, 196]],  // Teal (#41b6c4)
  [0.75, [34, 94, 168]],   // Blue (#225ea8)
  [1.00, [8, 29, 88]],     // Dark navy (#081d58)
];

/** Diverging RdBu palette when correlations include negative values. */
const HEATMAP_RD_BU = [
  [0.00, [178, 24, 43]],   // Negative: crimson red (#b2182b)
  [0.50, [247, 247, 247]], // 0.0: neutral off-white (#f7f7f7)
  [1.00, [33, 102, 172]],  // Positive: deep blue (#2166ac)
];

export function renderHeatmapSvg(samples, correlation, previewLabel) {
  const n = samples.length;
  const padLeft = 85;
  const padTop = 60;
  const cellSize = Math.max(35, Math.min(60, Math.floor(300 / n)));
  const matrixW = n * cellSize;
  const matrixH = n * cellSize;

  // Extract all off-diagonal correlations to compute dynamic data range
  const offDiag = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i !== j) {
        const v = correlation?.[i]?.[j];
        if (typeof v === 'number' && !Number.isNaN(v)) offDiag.push(v);
      }
    }
  }

  const hasOff = offDiag.length > 0;
  const minOff = hasOff ? Math.min(...offDiag) : 0;
  const maxOff = hasOff ? Math.max(...offDiag) : 1;
  const isDiverging = minOff < 0;

  // Palette selection and value mapper
  const stops = isDiverging ? HEATMAP_RD_BU : HEATMAP_YL_GN_BU;
  const maxAbs = Math.max(Math.abs(minOff), Math.abs(maxOff), 0.01);

  const getCellColor = (r, isDiag) => {
    if (isDiag) {
      // Diagonal is self-correlation (1.000). Use distinctive deep slate navy.
      return { fill: '#162238', text: '#ffffff' };
    }
    let t = 0.5;
    if (isDiverging) {
      // Diverging around 0: [-maxAbs, +maxAbs] -> [0, 1]
      t = (r + maxAbs) / (2 * maxAbs);
    } else {
      // Sequential: map observed off-diagonal range [minOff, maxOff] -> [0, 1]
      const span = maxOff - minOff;
      t = span > 0.001 ? (r - minOff) / span : 0.5;
    }
    const rgb = interpolateRgb(stops, t);
    const lum = rgbLuminance(rgb);
    const fill = `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
    const text = lum < 0.42 ? '#ffffff' : '#14171a';
    return { fill, text };
  };

  // Colorbar legend layout
  const legendW = 14;
  const legendX = padLeft + matrixW + 28;
  const legendY = padTop;
  const legendH = Math.min(matrixH, 200);
  const width = legendX + legendW + 65;
  const height = padTop + matrixH + 45;

  let rects = '';
  let labels = '';

  for (let i = 0; i < n; i++) {
    // Row label
    labels += `<text x="${padLeft - 8}" y="${padTop + i * cellSize + cellSize / 2 + 4}" text-anchor="end" font-size="12" fill="currentColor">${samples[i]}</text>`;
    // Col label
    labels += `<text x="${padLeft + i * cellSize + cellSize / 2}" y="${padTop - 8}" text-anchor="middle" font-size="12" fill="currentColor">${samples[i]}</text>`;

    for (let j = 0; j < n; j++) {
      const isDiag = i === j;
      const val = correlation?.[i]?.[j];
      const r = typeof val === 'number' && !Number.isNaN(val) ? val : (isDiag ? 1 : 0);
      const { fill, text } = getCellColor(r, isDiag);
      const x = padLeft + j * cellSize;
      const y = padTop + i * cellSize;

      rects += `<rect x="${x}" y="${y}" width="${cellSize}" height="${cellSize}" fill="${fill}" stroke="var(--line, #ccc)">` +
        `<title>${samples[i]} vs ${samples[j]}: ${r.toFixed(4)}${isDiag ? ' (self)' : ''}</title></rect>`;
      rects += `<text x="${x + cellSize / 2}" y="${y + cellSize / 2 + 4}" text-anchor="middle" font-size="11" font-weight="${isDiag ? '600' : '500'}" fill="${text}">${r.toFixed(3)}</text>`;
    }
  }

  // Build SVG gradient for the legend bar
  const gradStops = stops.map(([pos, rgb]) => {
    // Top of legend corresponds to t=1, bottom to t=0
    const pct = ((1 - pos) * 100).toFixed(0);
    return `<stop offset="${pct}%" stop-color="rgb(${rgb[0]},${rgb[1]},${rgb[2]})" />`;
  }).join('');

  const topLabel = isDiverging ? `+${maxAbs.toFixed(2)}` : maxOff.toFixed(3);
  const midLabel = isDiverging ? '0.00' : ((minOff + maxOff) / 2).toFixed(3);
  const btmLabel = isDiverging ? `-${maxAbs.toFixed(2)}` : minOff.toFixed(3);

  const legendSvg =
    `<defs>` +
      `<linearGradient id="heat-leg-grad" x1="0" y1="0" x2="0" y2="1">${gradStops}</linearGradient>` +
    `</defs>` +
    `<g class="heatmap-legend" transform="translate(${legendX}, ${legendY})">` +
      `<text x="0" y="-8" font-size="10" font-weight="600" fill="currentColor">Pearson r</text>` +
      `<rect x="0" y="0" width="${legendW}" height="${legendH}" fill="url(#heat-leg-grad)" stroke="var(--line, #ccc)" />` +
      `<text x="${legendW + 6}" y="10" font-size="10" fill="currentColor">${topLabel}</text>` +
      (legendH > 80 ? `<text x="${legendW + 6}" y="${legendH / 2 + 3}" font-size="9" fill="var(--muted, #666)">${midLabel}</text>` : '') +
      `<text x="${legendW + 6}" y="${legendH}" font-size="10" fill="currentColor">${btmLabel}</text>` +
      (hasOff && !isDiverging && maxOff < 0.95 ? `<text x="0" y="${legendH + 16}" font-size="9" fill="var(--muted, #666)">(off-diag)</text>` : '') +
    `</g>`;

  return `<svg class="heatmap" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="Sample correlation heatmap">` +
    `<text x="10" y="20" font-weight="600" font-size="14" fill="currentColor">Sample Correlation Matrix</text>` +
    `<text x="10" y="38" font-size="11" fill="var(--muted, #666)">${previewLabel}</text>` +
    labels + rects + legendSvg +
    `</svg>`;
}

/** Categorical palette for PCA and condition groups. */
const PALETTE = [
  '#0b5cad', // Blue
  '#d95f02', // Orange/Vermilion
  '#1b9e77', // Teal/Green
  '#7570b3', // Purple
  '#e7298a', // Magenta/Pink
  '#e6ab02', // Gold
  '#a6761d', // Brown
  '#66a61e', // Olive
];

/** Render SVG PCA scatter plot with condition group colors and legend. */
export function renderPcaSvg(pcaData, previewLabel, sampleGroups = new Map()) {
  const scores = pcaData.scores;
  const varExp = pcaData.variance_explained || [0, 0];
  const uniqueGroups = [...new Set(scores.map((s) => sampleGroups.get(s.sample) || ''))].filter(Boolean);
  const colorMap = new Map();
  uniqueGroups.forEach((g, i) => colorMap.set(g, PALETTE[i % PALETTE.length]));

  const hasGroups = uniqueGroups.length >= 2;
  const width = hasGroups ? 520 : 450;
  const height = 340;
  const pad = { top: 60, right: hasGroups ? 110 : 30, bottom: 50, left: 60 };
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;

  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const s of scores) {
    if (s.pc1 < minX) minX = s.pc1;
    if (s.pc1 > maxX) maxX = s.pc1;
    if (s.pc2 < minY) minY = s.pc2;
    if (s.pc2 > maxY) maxY = s.pc2;
  }
  // Add margin
  const spanX = Math.max(1, maxX - minX);
  const spanY = Math.max(1, maxY - minY);
  minX -= spanX * 0.15;
  maxX += spanX * 0.15;
  minY -= spanY * 0.15;
  maxY += spanY * 0.15;

  const scaleX = (x) => pad.left + ((x - minX) / (maxX - minX)) * plotW;
  const scaleY = (y) => pad.top + plotH - ((y - minY) / (maxY - minY)) * plotH;

  const pc1Label = `PC1 (${(varExp[0] * 100).toFixed(1)}%)`;
  const pc2Label = `PC2 (${(varExp[1] * 100).toFixed(1)}%)`;

  let points = '';
  for (const s of scores) {
    const cx = scaleX(s.pc1);
    const cy = scaleY(s.pc2);
    const grp = sampleGroups.get(s.sample) || '';
    const col = (grp && colorMap.get(grp)) || 'var(--accent, #0b5cad)';
    const tip = grp ? `${s.sample} (${grp}): PC1=${s.pc1.toFixed(2)}, PC2=${s.pc2.toFixed(2)}` : `${s.sample}: PC1=${s.pc1.toFixed(2)}, PC2=${s.pc2.toFixed(2)}`;
    points += `<circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="6" fill="${col}" stroke="#ffffff" stroke-width="1.5"><title>${tip}</title></circle>`;
    points += `<text x="${(cx + 8).toFixed(1)}" y="${(cy + 4).toFixed(1)}" font-size="11" fill="currentColor">${s.sample}</text>`;
  }

  // Zero axes if within range
  let zeroAxes = '';
  if (minX <= 0 && maxX >= 0) {
    const zx = scaleX(0);
    zeroAxes += `<line x1="${zx}" y1="${pad.top}" x2="${zx}" y2="${pad.top + plotH}" stroke="var(--line, #ccc)" stroke-dasharray="3,3" />`;
  }
  if (minY <= 0 && maxY >= 0) {
    const zy = scaleY(0);
    zeroAxes += `<line x1="${pad.left}" y1="${zy}" x2="${pad.left + plotW}" y2="${zy}" stroke="var(--line, #ccc)" stroke-dasharray="3,3" />`;
  }

  let legend = '';
  if (hasGroups) {
    const lx = pad.left + plotW + 15;
    let ly = pad.top + 15;
    legend += `<g class="pca-legend">`;
    legend += `<text x="${lx}" y="${ly - 5}" font-size="11" font-weight="600" fill="currentColor">Group</text>`;
    for (const g of uniqueGroups) {
      const col = colorMap.get(g);
      legend += `<circle cx="${lx + 6}" cy="${ly + 10}" r="5" fill="${col}" stroke="#ffffff" stroke-width="1" />`;
      legend += `<text x="${lx + 16}" y="${ly + 14}" font-size="11" fill="currentColor">${g}</text>`;
      ly += 22;
    }
    legend += `</g>`;
  }

  return `<svg class="pca-plot" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="Principal Component Analysis plot">` +
    `<text x="10" y="20" font-weight="600" font-size="13" fill="currentColor">PCA of log2(CPM+1), top 500 variable genes (not DESeq2-normalised)</text>` +
    `<text x="10" y="38" font-size="11" fill="var(--muted, #666)">${previewLabel}</text>` +
    `<rect x="${pad.left}" y="${pad.top}" width="${plotW}" height="${plotH}" fill="none" stroke="var(--line, #ccc)" />` +
    zeroAxes +
    `<text x="${pad.left + plotW / 2}" y="${height - 12}" text-anchor="middle" font-size="12" fill="currentColor">${pc1Label}</text>` +
    `<text x="16" y="${pad.top + plotH / 2}" text-anchor="middle" font-size="12" fill="currentColor" transform="rotate(-90 16 ${pad.top + plotH / 2})">${pc2Label}</text>` +
    points +
    legend +
    `</svg>`;
}

function normalCdfComplement(z) {
  // Abramowitz and Stegun formula 7.1.26 approximation
  if (z >= 8) return 0;
  const t = 1.0 / (1.0 + 0.2316419 * z);
  const poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + 1.330274429 * t))));
  const pdf = (1.0 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * z * z);
  return pdf * poly;
}

/** Compute two-group differential expression statistics for count matrix. */
export function computeDifferentialExpression(matrixTsv, samples, contrast) {
  const lines = matrixTsv.trim().split('\n');
  if (lines.length <= 1) return [];
  const header = lines[0].split('\t');
  const sampleCols = header.slice(1);

  const factor = contrast?.factor || 'condition';
  const trtLevel = contrast?.level;
  const refLevel = contrast?.reference;

  const sampleMap = new Map();
  for (const s of samples) {
    sampleMap.set(s.name, String(s.covariates?.[factor] ?? '').trim());
  }

  const trtIdxs = [];
  const refIdxs = [];
  for (let j = 0; j < sampleCols.length; j++) {
    const sName = sampleCols[j];
    const grp = sampleMap.get(sName);
    if (grp === trtLevel) trtIdxs.push(j);
    else if (grp === refLevel) refIdxs.push(j);
  }

  if (!trtIdxs.length || !refIdxs.length) return [];

  const genes = [];
  const rawMatrix = [];
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split('\t');
    genes.push(parts[0]);
    rawMatrix.push(parts.slice(1).map(Number));
  }
  const nGenes = genes.length;
  const nSamples = sampleCols.length;

  // Size factors (DESeq2 median of ratios)
  const sizeFactors = new Array(nSamples).fill(1);
  const logGeoMeans = new Float64Array(nGenes);
  let validGeoGenes = 0;

  for (let i = 0; i < nGenes; i++) {
    let sumLog = 0;
    let allPos = true;
    const row = rawMatrix[i];
    for (let j = 0; j < nSamples; j++) {
      if (row[j] <= 0) { allPos = false; break; }
      sumLog += Math.log(row[j]);
    }
    if (allPos) {
      logGeoMeans[i] = sumLog / nSamples;
      validGeoGenes++;
    } else {
      logGeoMeans[i] = -Infinity;
    }
  }

  if (validGeoGenes >= 10) {
    for (let j = 0; j < nSamples; j++) {
      const ratios = [];
      for (let i = 0; i < nGenes; i++) {
        if (logGeoMeans[i] !== -Infinity && rawMatrix[i][j] > 0) {
          ratios.push(rawMatrix[i][j] / Math.exp(logGeoMeans[i]));
        }
      }
      ratios.sort((a, b) => a - b);
      if (ratios.length) {
        const mid = Math.floor(ratios.length / 2);
        sizeFactors[j] = ratios.length % 2 === 1 ? ratios[mid] : (ratios[mid - 1] + ratios[mid]) / 2;
      }
    }
  } else {
    const totals = new Float64Array(nSamples);
    for (let i = 0; i < nGenes; i++) {
      for (let j = 0; j < nSamples; j++) totals[j] += rawMatrix[i][j];
    }
    const meanTot = totals.reduce((a, b) => a + b, 0) / nSamples;
    for (let j = 0; j < nSamples; j++) sizeFactors[j] = (totals[j] || 1) / (meanTot || 1);
  }

  const results = [];
  const nTrt = trtIdxs.length;
  const nRef = refIdxs.length;

  for (let i = 0; i < nGenes; i++) {
    const row = rawMatrix[i];
    let sumTrt = 0;
    let sumRef = 0;
    const trtVals = new Float64Array(nTrt);
    const refVals = new Float64Array(nRef);

    for (let k = 0; k < nTrt; k++) {
      const j = trtIdxs[k];
      const norm = row[j] / (sizeFactors[j] || 1);
      trtVals[k] = norm;
      sumTrt += norm;
    }
    for (let k = 0; k < nRef; k++) {
      const j = refIdxs[k];
      const norm = row[j] / (sizeFactors[j] || 1);
      refVals[k] = norm;
      sumRef += norm;
    }

    const meanTrt = sumTrt / nTrt;
    const meanRef = sumRef / nRef;
    const baseMean = (meanTrt + meanRef) / 2;

    const lfc = Math.log2(meanTrt + 0.5) - Math.log2(meanRef + 0.5);

    let varTrt = 0;
    for (let k = 0; k < nTrt; k++) varTrt += (trtVals[k] - meanTrt) ** 2;
    varTrt = nTrt > 1 ? varTrt / (nTrt - 1) : 0;

    let varRef = 0;
    for (let k = 0; k < nRef; k++) varRef += (refVals[k] - meanRef) ** 2;
    varRef = nRef > 1 ? varRef / (nRef - 1) : 0;

    const se = Math.sqrt(varTrt / nTrt + varRef / nRef + 0.05);
    const stat = lfc / (se || 1);
    const z = Math.abs(stat);
    const pval = baseMean === 0 ? 1.0 : Math.min(1.0, 2.0 * normalCdfComplement(z));

    results.push({
      gene: genes[i],
      baseMean,
      log2FoldChange: lfc,
      stat,
      pvalue: pval,
      padj: 1.0,
    });
  }

  // Benjamini-Hochberg adjustment
  const indexed = results.map((r, i) => ({ idx: i, pvalue: r.pvalue }));
  indexed.sort((a, b) => a.pvalue - b.pvalue);
  const m = indexed.length;
  let cumMin = 1.0;

  for (let k = m - 1; k >= 0; k--) {
    const rawP = indexed[k].pvalue;
    const rank = k + 1;
    const adj = Math.min(1.0, (rawP * m) / rank);
    cumMin = Math.min(cumMin, adj);
    results[indexed[k].idx].padj = cumMin;
  }

  return results;
}

/** Render SVG Volcano plot. */
export function renderVolcanoSvg(deRows, aliasIndex = new Map(), { fcCutoff = 1.0, padjCutoff = 0.05 } = {}) {
  const width = 560;
  const height = 400;
  const pad = { top: 45, right: 30, bottom: 50, left: 60 };
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;

  let maxAbsFc = 2.5;
  let maxNegLogP = 2.0;

  for (const r of deRows) {
    if (Number.isFinite(r.log2FoldChange)) {
      const a = Math.abs(r.log2FoldChange);
      if (a > maxAbsFc) maxAbsFc = a;
    }
    const p = r.padj ?? r.pvalue;
    if (Number.isFinite(p) && p > 0) {
      const nlp = -Math.log10(p);
      if (nlp > maxNegLogP) maxNegLogP = nlp;
    }
  }

  maxNegLogP = Math.min(15, Math.ceil(maxNegLogP * 1.15));
  const spanX = Math.ceil(maxAbsFc * 1.1 * 2) / 2;
  const minX = -spanX;
  const maxX = spanX;
  const minY = 0;
  const maxY = Math.max(2, maxNegLogP);

  const scaleX = (x) => pad.left + ((x - minX) / (maxX - minX)) * plotW;
  const scaleY = (y) => pad.top + plotH - ((y - minY) / (maxY - minY)) * plotH;

  const pLineY = scaleY(-Math.log10(padjCutoff));
  const fcLeftX = scaleX(-fcCutoff);
  const fcRightX = scaleX(fcCutoff);

  let grid = '';
  if (pLineY >= pad.top && pLineY <= pad.top + plotH) {
    grid += `<line x1="${pad.left}" y1="${pLineY.toFixed(1)}" x2="${pad.left + plotW}" y2="${pLineY.toFixed(1)}" stroke="#e41a1c" stroke-dasharray="3,3" stroke-width="1"><title>padj = ${padjCutoff}</title></line>`;
  }
  if (fcLeftX >= pad.left && fcLeftX <= pad.left + plotW) {
    grid += `<line x1="${fcLeftX.toFixed(1)}" y1="${pad.top}" x2="${fcLeftX.toFixed(1)}" y2="${pad.top + plotH}" stroke="var(--line, #999)" stroke-dasharray="3,3" stroke-width="1"><title>log2FC = -${fcCutoff}</title></line>`;
  }
  if (fcRightX >= pad.left && fcRightX <= pad.left + plotW) {
    grid += `<line x1="${fcRightX.toFixed(1)}" y1="${pad.top}" x2="${fcRightX.toFixed(1)}" y2="${pad.top + plotH}" stroke="var(--line, #999)" stroke-dasharray="3,3" stroke-width="1"><title>log2FC = ${fcCutoff}</title></line>`;
  }

  let points = '';
  let nUp = 0;
  let nDown = 0;
  let nNs = 0;

  const topUp = [];
  const topDown = [];

  for (const r of deRows) {
    const lfc = r.log2FoldChange;
    const p = r.padj ?? r.pvalue;
    if (!Number.isFinite(lfc) || !Number.isFinite(p) || p <= 0) continue;

    const nlp = Math.min(maxY, -Math.log10(p));
    const cx = scaleX(lfc);
    const cy = scaleY(nlp);

    const isSig = p < padjCutoff;
    const isUp = isSig && lfc >= fcCutoff;
    const isDown = isSig && lfc <= -fcCutoff;

    let col = '#999999';
    let opacity = '0.35';
    let radius = '3.5';

    if (isUp) {
      col = '#e41a1c';
      opacity = '0.85';
      radius = '4.5';
      nUp++;
      topUp.push({ ...r, nlp, cx, cy });
    } else if (isDown) {
      col = '#0b5cad';
      opacity = '0.85';
      radius = '4.5';
      nDown++;
      topDown.push({ ...r, nlp, cx, cy });
    } else {
      nNs++;
    }

    const alias = aliasIndex.get(r.gene);
    const sym = alias?.symbol ? ` (${alias.symbol})` : '';
    const tip = `${r.gene}${sym}: log2FC=${lfc.toFixed(2)}, padj=${p.toExponential(2)}`;
    points += `<circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${radius}" fill="${col}" fill-opacity="${opacity}"><title>${escSvg(tip)}</title></circle>`;
  }

  topUp.sort((a, b) => (b.nlp * b.log2FoldChange) - (a.nlp * a.log2FoldChange));
  topDown.sort((a, b) => (b.nlp * Math.abs(b.log2FoldChange)) - (a.nlp * Math.abs(a.log2FoldChange)));

  function layoutSideLabels(items, isRight) {
    if (!items.length) return [];
    const list = items.map((item) => {
      const alias = aliasIndex.get(item.gene);
      const txt = alias?.symbol || item.gene;
      const links = geneExternalLinks(item.gene, alias);
      const p = item.padj ?? item.pvalue;
      let tip = `${item.gene}${alias?.symbol ? ` (${alias.symbol})` : ''}: log2FC=${item.log2FoldChange.toFixed(2)}, padj=${p.toExponential(2)}`;
      for (const l of links) tip += `\nOpen ${l.name}: ${l.url}`;
      return {
        item,
        txt,
        alias,
        links,
        tip,
        isRight,
        origX: isRight ? item.cx + 8 : item.cx - 8,
        origY: item.cy + 3,
        x: isRight ? item.cx + 8 : item.cx - 8,
        y: item.cy + 3,
        anchor: isRight ? 'start' : 'end',
      };
    });

    // Sort by y coordinate (top to bottom)
    list.sort((a, b) => a.y - b.y);

    const minGapY = 14;
    const minY = pad.top + 10;
    const maxY = pad.top + plotH - 5;

    // Iterative repulsion passes to eliminate vertical overlap
    for (let iter = 0; iter < 12; iter++) {
      for (let i = 0; i < list.length - 1; i++) {
        const a = list[i];
        const b = list[i + 1];
        const overlap = minGapY - (b.y - a.y);
        if (overlap > 0) {
          a.y -= overlap * 0.55;
          b.y += overlap * 0.55;
        }
      }
      for (let i = list.length - 1; i > 0; i--) {
        const a = list[i - 1];
        const b = list[i];
        const overlap = minGapY - (b.y - a.y);
        if (overlap > 0) {
          a.y -= overlap * 0.55;
          b.y += overlap * 0.55;
        }
      }
      for (const l of list) {
        if (l.y < minY) l.y = minY;
        if (l.y > maxY) l.y = maxY;
      }
    }

    // Ensure strictly no overlap after bounding
    for (let i = 1; i < list.length; i++) {
      if (list[i].y - list[i - 1].y < minGapY) {
        list[i].y = Math.min(maxY, list[i - 1].y + minGapY);
      }
    }

    // Offset / jitter horizontally if labels were displaced
    for (let i = 0; i < list.length; i++) {
      const l = list[i];
      const dy = Math.abs(l.y - l.origY);
      if (dy > 3) {
        const stagger = (i % 2 === 1) ? 6 : 0;
        l.x = isRight
          ? Math.min(width - pad.right - 5, l.origX + 4 + stagger)
          : Math.max(pad.left + 5, l.origX - 4 - stagger);
      }
    }

    return list;
  }

  let lines = '';
  let labels = '';
  let hitAreas = '';
  const layouted = [
    ...layoutSideLabels(topUp.slice(0, 5), true),
    ...layoutSideLabels(topDown.slice(0, 5), false),
  ];

  for (const l of layouted) {
    const item = l.item;
    const isRight = l.isRight;
    const dy = Math.abs(l.y - l.origY);
    const dx = Math.abs(l.x - l.origX);
    const isDisplaced = dy > 3 || dx > 3;

    if (isDisplaced) {
      const startX = item.cx + (isRight ? 4.5 : -4.5);
      const startY = item.cy;
      const endX = isRight ? l.x - 3 : l.x + 3;
      const endY = l.y - 3;
      lines += `<line x1="${startX.toFixed(1)}" y1="${startY.toFixed(1)}" x2="${endX.toFixed(1)}" y2="${endY.toFixed(1)}" stroke="var(--line, #777)" stroke-width="0.8" stroke-dasharray="2,2" stroke-opacity="0.85"><title>${escSvg(l.tip)}</title></line>`;
    }

    const textAttrs = `x="${l.x.toFixed(1)}" y="${l.y.toFixed(1)}" font-size="10" font-weight="600" text-anchor="${l.anchor}" fill="currentColor" paint-order="stroke fill" stroke="var(--bg, #fff)" stroke-width="2.5" stroke-linejoin="round"`;
    if (l.links.length) {
      const primary = l.links[0];
      // The label and the point itself link out; the mouseover names every linked database.
      labels += `<a class="volcano-gene-link" href="${escSvg(primary.url)}" target="_blank" rel="noopener noreferrer"><title>${escSvg(l.tip)}</title><text ${textAttrs}>${escSvg(l.txt)}</text></a>`;
      hitAreas += `<a class="volcano-gene-link" href="${escSvg(primary.url)}" target="_blank" rel="noopener noreferrer"><title>${escSvg(l.tip)}</title><circle cx="${item.cx.toFixed(1)}" cy="${item.cy.toFixed(1)}" r="8" fill="#000" fill-opacity="0" stroke="none"/></a>`;
    } else {
      labels += `<text ${textAttrs}><title>${escSvg(l.tip)}</title>${escSvg(l.txt)}</text>`;
    }
  }

  const legend = `<g class="volcano-legend">` +
    `<circle cx="${pad.left + 15}" cy="22" r="4.5" fill="#e41a1c" />` +
    `<text x="${pad.left + 25}" y="26" font-size="11" fill="currentColor">Up: ${nUp.toLocaleString()}</text>` +
    `<circle cx="${pad.left + 115}" cy="22" r="4.5" fill="#0b5cad" />` +
    `<text x="${pad.left + 125}" y="26" font-size="11" fill="currentColor">Down: ${nDown.toLocaleString()}</text>` +
    `<circle cx="${pad.left + 225}" cy="22" r="3.5" fill="#999999" />` +
    `<text x="${pad.left + 235}" y="26" font-size="11" fill="currentColor">Not sig: ${nNs.toLocaleString()}</text>` +
    `</g>`;

  return `<svg class="volcano-plot" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="Volcano plot">` +
    legend +
    `<rect x="${pad.left}" y="${pad.top}" width="${plotW}" height="${plotH}" fill="none" stroke="var(--line, #ccc)" />` +
    grid +
    lines +
    points +
    labels +
    hitAreas +
    `<text x="${pad.left + plotW / 2}" y="${height - 12}" text-anchor="middle" font-size="12" fill="currentColor">log2(Fold Change)</text>` +
    `<text x="16" y="${pad.top + plotH / 2}" text-anchor="middle" font-size="12" fill="currentColor" transform="rotate(-90 16 ${pad.top + plotH / 2})">-log10(p-adj)</text>` +
    `</svg>`;
}

export function createResultsStep({ $, app, shell }) {
  let currentResults = null;
  let aliasIndex = new Map();
  let matrixRows = [];
  let sortCol = 0;
  let sortDir = 'asc';

  function previewLabelText() {
    return app.state?.preview?.label || 'preview of the first N read pairs';
  }

  function getPreviewN() {
    return app.state?.preview?.n || app.state?.samples?.find((s) => s.status === 'done')?.previewN || 0;
  }

  function countsFileName() {
    const n = getPreviewN();
    return n ? `counts.preview-${n}.tsv` : 'counts.preview.tsv';
  }

  function getFolder() {
    const val = $('reads-folder')?.value.trim() || '';
    return val.replace(/\/+$/, '');
  }

  function getSamplesTsv() {
    if (!currentResults) return '';
    const folder = getFolder();
    const doneSamples = currentResults.samples.filter((s) => s.status === 'done');
    if (!doneSamples.length) return '';

    // Covariates columns
    const covKeys = Object.keys(doneSamples[0].covariates || {});
    const head = ['sample', 'fastq1', 'fastq2', ...covKeys];
    const lines = [head.join('\t')];

    for (const s of doneSamples) {
      const f1 = s.fastq1 ? (folder ? `${folder}/${s.fastq1}` : s.fastq1) : '';
      const f2 = s.fastq2 ? (folder ? `${folder}/${s.fastq2}` : s.fastq2) : '';
      const row = [s.name, f1, f2, ...covKeys.map((k) => s.covariates[k] ?? '')];
      lines.push(row.join('\t'));
    }
    return lines.join('\n') + '\n';
  }

  function getDeCommand() {
    if (!currentResults) return '';
    const doneSamples = currentResults.samples.filter((s) => s.status === 'done');
    const covKeys = Object.keys(doneSamples[0]?.covariates || {});
    const colName = covKeys.includes('condition') ? 'condition' : covKeys[0] || 'condition';
    const values = [...new Set(doneSamples.map((s) => String(s.covariates?.[colName] ?? '').trim()).filter(Boolean))];
    const design = `~ ${colName}`;
    const contrast = values.length >= 2 ? `${colName},${values[1]},${values[0]}` : `${colName},trt,ctl`;
    return de_command('samples.tsv', countsFileName(), design, contrast) + '\n';
  }

  function getCountsTsv() {
    return currentResults?.matrixTsv || '';
  }

  function getGeneInfoTsv() {
    return app.state?.geneInfoTsv || currentResults?.geneInfoTsv || '';
  }

  // Expose on window for programmatic access in tests
  window.qrnamapResults = {
    getCountsTsv,
    getSamplesTsv,
    getDeCommand,
    getGeneInfoTsv,
    getDeTsv,
  };

  function updateReadsFolderNote() {
    const note = $('reads-folder-note');
    if (!note) return;
    const folder = getFolder();
    if (folder) {
      note.textContent = `Using folder "${folder}" in samples.tsv and de-command.txt.`;
    } else {
      note.textContent = 'Using bare file names in samples.tsv and de-command.txt.';
    }
  }

  $('reads-folder')?.addEventListener('input', updateReadsFolderNote);

  // Download buttons
  $('dl-counts')?.addEventListener('click', () => {
    downloadFile(countsFileName(), getCountsTsv());
  });

  $('dl-sheet')?.addEventListener('click', () => {
    downloadFile('samples.tsv', getSamplesTsv());
  });

  $('dl-de')?.addEventListener('click', () => {
    downloadFile('de-command.txt', getDeCommand());
  });

  $('dl-gene-info')?.addEventListener('click', () => {
    downloadFile('gene_info.tsv', getGeneInfoTsv());
  });

  // Render mapping summary
  function renderMappingSummary(samples, preview) {
    const table = $('results-mapping-table');
    const labelEl = $('mapping-summary-label');
    if (labelEl) labelEl.textContent = preview;

    const done = samples.filter((s) => s.status === 'done');
    let html = '<thead><tr><th>Sample</th><th>Reads processed</th><th>Pseudoaligned</th><th>Pseudoaligned %</th><th>Note</th></tr></thead><tbody>';
    for (const s of done) {
      const note = s.allPairs ? `all pairs (${fmt(s.nPairs)})` : preview;
      html += `<tr>` +
        `<td>${s.name}</td>` +
        `<td>${fmt(s.processed)}</td>` +
        `<td>${fmt(s.pseudoaligned)}</td>` +
        `<td>${(s.fraction * 100).toFixed(1)}%</td>` +
        `<td>${note}</td>` +
        `</tr>`;
    }
    html += '</tbody>';
    table.innerHTML = html;
  }

  // Render count table
  function renderCountTable(sampleNames, preview) {
    const table = $('results-count-table');
    const labelEl = $('count-table-label');
    if (labelEl) labelEl.textContent = preview;

    const query = $('count-search')?.value.trim().toLowerCase() || '';

    let visibleRows = matrixRows;
    if (query) {
      visibleRows = matrixRows.filter((r) => {
        const info = aliasIndex.get(r.gene);
        const tokens = info?.tokens || r.gene.toLowerCase();
        return tokens.includes(query);
      });
    }

    // Sort
    visibleRows.sort((a, b) => {
      let va = sortCol === 0 ? a.gene : a.counts[sortCol - 1] ?? 0;
      let vb = sortCol === 0 ? b.gene : b.counts[sortCol - 1] ?? 0;
      if (sortCol > 0) {
        va = Number(va) || 0;
        vb = Number(vb) || 0;
      }
      if (va < vb) return sortDir === 'asc' ? -1 : 1;
      if (va > vb) return sortDir === 'asc' ? 1 : -1;
      return 0;
    });

    let head = '<thead><tr>';
    head += `<th data-col="0" data-sort="${sortCol === 0 ? sortDir : ''}">Gene</th>`;
    sampleNames.forEach((s, idx) => {
      const c = idx + 1;
      head += `<th data-col="${c}" data-sort="${sortCol === c ? sortDir : ''}">${s}</th>`;
    });
    head += '</tr></thead>';

    // Limit initial display rows to 300 for UI responsiveness
    const limit = Math.min(visibleRows.length, 300);
    let body = '<tbody>';
    for (let i = 0; i < limit; i++) {
      const r = visibleRows[i];
      body += `<tr><td>${r.gene}</td>`;
      for (const val of r.counts) {
        body += `<td>${fmt(val)}</td>`;
      }
      body += '</tr>';
    }
    if (visibleRows.length === 0) {
      body += `<tr><td colspan="${sampleNames.length + 1}" class="small">No matching genes found.</td></tr>`;
    } else if (visibleRows.length > limit) {
      body += `<tr><td colspan="${sampleNames.length + 1}" class="small">Showing first ${limit} of ${visibleRows.length} genes. Use search above to find specific genes.</td></tr>`;
    }
    body += '</tbody>';

    table.innerHTML = head + body;

    table.querySelectorAll('th').forEach((th) => {
      th.addEventListener('click', () => {
        const col = Number(th.dataset.col);
        if (sortCol === col) {
          sortDir = sortDir === 'asc' ? 'desc' : 'asc';
        } else {
          sortCol = col;
          sortDir = 'asc';
        }
        renderCountTable(sampleNames, preview);
      });
    });
  }

  $('count-search')?.addEventListener('input', () => {
    if (currentResults) {
      renderCountTable(currentResults.sampleNames, previewLabelText());
    }
  });

  // Render QC plots
  function renderQc(doneSamples, qcData, preview) {
    const heatContainer = $('heatmap-container');
    const pcaContainer = $('pca-container');

    // Heatmap
    if (doneSamples.length < 2 || !qcData?.correlation || qcData.correlation_error) {
      heatContainer.innerHTML = `<p class="qc-message">${qcData?.correlation_error || 'Correlation heatmap needs at least two samples.'}</p>`;
    } else {
      heatContainer.innerHTML = renderHeatmapSvg(doneSamples.map((s) => s.name), qcData.correlation, preview);
    }

    // PCA
    if (doneSamples.length < 3 || !qcData?.pca || qcData.pca_error) {
      pcaContainer.innerHTML = `<p class="qc-message">${qcData?.pca_error || 'PCA needs at least three samples.'}</p>`;
    } else {
      const sampleGroups = new Map();
      const covKeys = Object.keys(doneSamples[0]?.covariates || {});
      const groupKey = covKeys.includes('condition') ? 'condition' : (covKeys.includes('genotype') ? 'genotype' : covKeys[0]);
      if (groupKey) {
        for (const s of doneSamples) {
          sampleGroups.set(s.name, String(s.covariates?.[groupKey] ?? '').trim() || 'default');
        }
      }
      pcaContainer.innerHTML = renderPcaSvg(qcData.pca, preview, sampleGroups);
    }
  }

  let currentDeResults = [];
  let currentContrast = null;

  function updateVolcanoPlot() {
    const container = $('volcano-container');
    const summaryEl = $('volcano-summary');
    if (!container || !currentDeResults.length) return;

    const fcCutoff = Math.max(0.1, parseFloat($('volcano-fc-cutoff')?.value) || 1.0);
    const padjCutoff = Math.max(0.0001, parseFloat($('volcano-padj-cutoff')?.value) || 0.05);

    container.innerHTML = renderVolcanoSvg(currentDeResults, aliasIndex, { fcCutoff, padjCutoff });

    const nUp = currentDeResults.filter((r) => r.padj < padjCutoff && r.log2FoldChange >= fcCutoff).length;
    const nDown = currentDeResults.filter((r) => r.padj < padjCutoff && r.log2FoldChange <= -fcCutoff).length;
    const nTested = currentDeResults.length;

    if (summaryEl) {
      summaryEl.textContent = `Significant (|log2 FC| \u2265 ${fcCutoff.toFixed(1)}, padj < ${padjCutoff}): ${nUp} upregulated, ${nDown} downregulated, ${nTested - nUp - nDown} unchanged out of ${nTested.toLocaleString()} genes.`;
    }
  }

  function renderVolcano(doneSamples, matrixTsv) {
    const container = $('volcano-container');
    const labelEl = $('de-comparison-label');
    const summaryEl = $('volcano-summary');
    if (!container) return;

    const covKeys = Object.keys(doneSamples[0]?.covariates || {});
    const colName = covKeys.includes('condition') ? 'condition' : (covKeys.includes('genotype') ? 'genotype' : covKeys[0]);

    if (!colName) {
      if (labelEl) labelEl.textContent = 'No covariate column found for differential expression.';
      container.innerHTML = '<p class="qc-message">Add covariate groups (e.g. condition or genotype) to compare samples.</p>';
      if (summaryEl) summaryEl.textContent = '';
      currentDeResults = [];
      return;
    }

    const values = [...new Set(doneSamples.map((s) => String(s.covariates?.[colName] ?? '').trim()).filter(Boolean))];
    if (values.length < 2) {
      if (labelEl) labelEl.textContent = `Only one group (${values[0] || 'unknown'}) found in column "${colName}".`;
      container.innerHTML = '<p class="qc-message">Differential expression comparison requires samples from at least two groups.</p>';
      if (summaryEl) summaryEl.textContent = '';
      currentDeResults = [];
      return;
    }

    currentContrast = { factor: colName, level: values[1], reference: values[0] };
    if (labelEl) {
      labelEl.textContent = `Comparison: "${currentContrast.level}" vs "${currentContrast.reference}" (~ ${currentContrast.factor})`;
    }

    currentDeResults = computeDifferentialExpression(matrixTsv, doneSamples, currentContrast);
    updateVolcanoPlot();
  }

  $('volcano-fc-cutoff')?.addEventListener('input', () => updateVolcanoPlot());
  $('volcano-padj-cutoff')?.addEventListener('input', () => updateVolcanoPlot());

  function getDeTsv() {
    if (!currentDeResults.length) return '';
    const fcCutoff = Math.max(0.1, parseFloat($('volcano-fc-cutoff')?.value) || 1.0);
    const padjCutoff = Math.max(0.0001, parseFloat($('volcano-padj-cutoff')?.value) || 0.05);

    const head = ['gene', 'symbol', 'baseMean', 'log2FoldChange', 'stat', 'pvalue', 'padj', 'significant'];
    const lines = [head.join('\t')];

    for (const r of currentDeResults) {
      const alias = aliasIndex.get(r.gene);
      const sym = alias?.symbol || '';
      let sig = 'no';
      if (r.padj < padjCutoff && r.log2FoldChange >= fcCutoff) sig = 'up';
      else if (r.padj < padjCutoff && r.log2FoldChange <= -fcCutoff) sig = 'down';

      lines.push([
        r.gene,
        sym,
        r.baseMean.toFixed(2),
        r.log2FoldChange.toFixed(4),
        r.stat.toFixed(4),
        r.pvalue.toExponential(4),
        r.padj.toExponential(4),
        sig,
      ].join('\t'));
    }
    return lines.join('\n') + '\n';
  }

  $('dl-de-tsv')?.addEventListener('click', () => {
    const tsv = getDeTsv();
    if (tsv) {
      downloadBlob(tsv, 'differential_expression.tsv', 'text/tab-separated-values');
    }
  });

  // Display results from completed samples
  async function displayResults() {
    const doneSamples = (app.state.samples || []).filter((s) => s.status === 'done');
    if (!doneSamples.length) return;

    const preview = previewLabelText();
    const prevEl = $('results-preview-label');
    if (prevEl) prevEl.textContent = preview;

    // Parse alias index from geneInfoTsv
    aliasIndex = parseAliasIndex(app.state.geneInfoTsv);

    // Assemble matrix
    const matrixPayload = JSON.stringify({
      tx2gene: app.state.tx2gene,
      samples: doneSamples.map((s) => ({ name: s.name, counts: s.geneCounts })),
    });
    const assembled = JSON.parse(assemble_matrix(matrixPayload));
    const matrixTsv = assembled.tsv;
    const sampleNames = doneSamples.map((s) => s.name);

    // Parse rows for count table
    const lines = matrixTsv.trim().split('\n');
    matrixRows = lines.slice(1).map((line) => {
      const parts = line.split('\t');
      return { gene: parts[0], counts: parts.slice(1).map(Number) };
    });

    // Run QC
    let qcData = null;
    try {
      qcData = JSON.parse(qc(matrixTsv));
    } catch (e) {
      qcData = { correlation_error: e.message || String(e), pca_error: e.message || String(e) };
    }

    currentResults = {
      genomeKey: app.state.genome?.accession || app.state.genome?.label,
      previewN: getPreviewN(),
      matrixTsv,
      samples: app.state.samples,
      sampleNames,
      qcData,
      geneInfoTsv: app.state.geneInfoTsv,
    };

    updateReadsFolderNote();
    renderMappingSummary(app.state.samples, preview);
    renderCountTable(sampleNames, preview);
    renderQc(doneSamples, qcData, preview);
    renderVolcano(doneSamples, matrixTsv);

    // Save to persist
    const storeEntry = {
      genomeKey: currentResults.genomeKey,
      previewN: currentResults.previewN,
      matrixTsv: currentResults.matrixTsv,
      samples: currentResults.samples.map((s) => ({
        name: s.name,
        covariates: s.covariates,
        fastq1: s.fastq1,
        fastq2: s.fastq2,
        status: s.status,
        nPairs: s.nPairs,
        processed: s.processed,
        pseudoaligned: s.pseudoaligned,
        fraction: s.fraction,
        allPairs: s.allPairs,
        previewN: s.previewN,
        readLength: s.readLength,
        geneCounts: s.geneCounts,
      })),
      preview: app.state.preview,
      geneInfoTsv: app.state.geneInfoTsv,
    };
    writeStore({ [STORE_KEY]: storeEntry }).catch(() => {});

    shell.completeStep('samples');
  }

  document.addEventListener('qrnamap:samples-done', () => {
    displayResults();
  });

  // Reset when worker replaced or genome changed
  window.addEventListener('qrnamap:worker-replaced', () => {
    currentResults = null;
    matrixRows = [];
    $('results-mapping-table').innerHTML = '';
    $('results-count-table').innerHTML = '';
    $('heatmap-container').innerHTML = '';
    $('pca-container').innerHTML = '';
    $('volcano-container').innerHTML = '';
    $('volcano-summary').textContent = '';
    currentDeResults = [];
    writeStore({ [STORE_KEY]: null }).catch(() => {});
  });

  // Restore saved results on load if matching genome & N
  async function tryRestore() {
    try {
      const stored = (await readStore([STORE_KEY]))[STORE_KEY];
      if (!stored) return;
      const currentGenomeKey = app.state.genome?.accession || app.state.genome?.label;
      if (!currentGenomeKey || stored.genomeKey !== currentGenomeKey) return;
      const currentN = getPreviewN();
      if (currentN && stored.previewN && stored.previewN !== currentN) return;

      // Restore
      app.state.samples = stored.samples;
      app.state.preview = stored.preview || { n: stored.previewN, label: `preview of the first ${fmt(stored.previewN)} pairs` };
      app.state.geneInfoTsv = stored.geneInfoTsv || app.state.geneInfoTsv;
      displayResults();
    } catch (e) {
      // ignore
    }
  }

  document.addEventListener('qrnamap:genome', () => {
    if (app.state.indexReady && !currentResults) {
      tryRestore();
    }
  });

  if (app.state.indexReady) {
    tryRestore();
  }
}
