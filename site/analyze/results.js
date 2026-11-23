// Step 3: results view, downloads and plots.
import { assemble_matrix, qc, de_command, validate_sheet } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { readStore, writeStore } from '../persist.js?v=__BUILD__';

const STORE_KEY = 'results';

/** Parse gene_info.tsv into a map of gene -> lowercased searchable tokens. */
export function parseAliasIndex(geneInfoTsv) {
  const map = new Map();
  if (!geneInfoTsv) return map;
  const lines = geneInfoTsv.trim().split('\n');
  if (lines.length <= 1) return map;
  for (const line of lines.slice(1)) {
    const cols = line.split('\t');
    if (!cols.length) continue;
    const gene = cols[0];
    // gene, symbol, gene_id, old_locus_tag, transcripts, proteins
    const tokens = cols.slice(0, 6).filter(Boolean).join(' ').toLowerCase();
    map.set(gene, tokens);
  }
  return map;
}

/** Format a number with commas. */
const fmt = (n) => Number(n).toLocaleString('en-US');

/** Helper to trigger file download in browser. */
function downloadFile(name, text) {
  const blob = new Blob([text], { type: 'text/tab-separated-values;charset=utf-8' });
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
export function renderHeatmapSvg(samples, correlation, previewLabel) {
  const n = samples.length;
  const padLeft = 80;
  const padTop = 60;
  const cellSize = Math.max(35, Math.min(60, Math.floor(300 / n)));
  const width = padLeft + n * cellSize + 20;
  const height = padTop + n * cellSize + 40;

  let rects = '';
  let labels = '';

  for (let i = 0; i < n; i++) {
    // Row label
    labels += `<text x="${padLeft - 8}" y="${padTop + i * cellSize + cellSize / 2 + 4}" text-anchor="end" font-size="12" fill="currentColor">${samples[i]}</text>`;
    // Col label
    labels += `<text x="${padLeft + i * cellSize + cellSize / 2}" y="${padTop - 8}" text-anchor="middle" font-size="12" fill="currentColor">${samples[i]}</text>`;

    for (let j = 0; j < n; j++) {
      const val = correlation[i][j];
      const r = typeof val === 'number' && !Number.isNaN(val) ? val : 0;
      // Interpolate color from light (0.5 or lower) to solid blue (1.0)
      const alpha = Math.max(0.1, Math.min(1.0, (r - 0.5) / 0.5));
      const fill = `rgba(11, 92, 173, ${alpha.toFixed(2)})`;
      const textColor = alpha > 0.65 ? '#ffffff' : '#1d2024';
      const x = padLeft + j * cellSize;
      const y = padTop + i * cellSize;

      rects += `<rect x="${x}" y="${y}" width="${cellSize}" height="${cellSize}" fill="${fill}" stroke="var(--line, #ccc)">` +
        `<title>${samples[i]} vs ${samples[j]}: ${r.toFixed(4)}</title></rect>`;
      rects += `<text x="${x + cellSize / 2}" y="${y + cellSize / 2 + 4}" text-anchor="middle" font-size="11" fill="${textColor}">${r.toFixed(3)}</text>`;
    }
  }

  return `<svg class="heatmap" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="Sample correlation heatmap">` +
    `<text x="10" y="20" font-weight="600" font-size="14" fill="currentColor">Sample Correlation Matrix</text>` +
    `<text x="10" y="38" font-size="11" fill="var(--muted, #666)">${previewLabel}</text>` +
    labels + rects +
    `</svg>`;
}

/** Render SVG PCA scatter plot. */
export function renderPcaSvg(pcaData, previewLabel) {
  const scores = pcaData.scores;
  const varExp = pcaData.variance_explained || [0, 0];
  const width = 450;
  const height = 340;
  const pad = { top: 60, right: 30, bottom: 50, left: 60 };
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
    points += `<circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="6" fill="var(--accent, #0b5cad)" stroke="#ffffff" stroke-width="1.5"><title>${s.sample}: PC1=${s.pc1.toFixed(2)}, PC2=${s.pc2.toFixed(2)}</title></circle>`;
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

  return `<svg class="pca-plot" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" role="img" aria-label="Principal Component Analysis plot">` +
    `<text x="10" y="20" font-weight="600" font-size="13" fill="currentColor">PCA of log2(CPM+1), top 500 variable genes (not DESeq2-normalised)</text>` +
    `<text x="10" y="38" font-size="11" fill="var(--muted, #666)">${previewLabel}</text>` +
    `<rect x="${pad.left}" y="${pad.top}" width="${plotW}" height="${plotH}" fill="none" stroke="var(--line, #ccc)" />` +
    zeroAxes +
    `<text x="${pad.left + plotW / 2}" y="${height - 12}" text-anchor="middle" font-size="12" fill="currentColor">${pc1Label}</text>` +
    `<text x="16" y="${pad.top + plotH / 2}" text-anchor="middle" font-size="12" fill="currentColor" transform="rotate(-90 16 ${pad.top + plotH / 2})">${pc2Label}</text>` +
    points +
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
        const tokens = aliasIndex.get(r.gene) || r.gene.toLowerCase();
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
      pcaContainer.innerHTML = renderPcaSvg(qcData.pca, preview);
    }
  }

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
