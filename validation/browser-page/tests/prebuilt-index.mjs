// Measurement: time to "index ready" with a pre-built index from the S3 mirror, against the index built in the page.
//   DIST=<abs dist-site> ENGINE=chromium RUNS=3 GENOMES=GCF_000146045.2,GCF_000002655.1 \
//     node validation/browser-page/tests/prebuilt-index.mjs > results.tsv
// The page is served on http://localhost:8000: an origin the bucket's CORS rule allows (so this needs no CORS change).
// For each genome and run: a fresh browser context, `loadCatalogAccession` (the same code as the "Download and build" button),
// timed from the call to "index ready", once with the pre-built index (default) and once with `?prebuilt=0` (NCBI download
// and kallisto index in the page). Also records the bytes fetched per host, and whether the pre-built path really
// was used (the page falls back to building silently when the download fails).
// Compatibility check (CHECK_COUNTS=1, S. cerevisiae only): the fixture reads s1 (simulated from the yeast mitochondrial
// genome, ../fixtures) are counted with both indexes; the gene counts must be equal.
// Output: TSV on stdout, one row per run. Progress on stderr.
import http from 'node:http';
import { readFileSync, existsSync, statSync } from 'node:fs';
import { join, resolve, extname, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(process.env.DIST || resolve(here, '../../../dist-site'));
const FIX = resolve(here, '../fixtures');
const ENGINE = process.env.ENGINE || 'chromium';
const RUNS = Number(process.env.RUNS || 3);
const GENOMES = (process.env.GENOMES || 'GCF_000146045.2,GCF_000149245.1,GCF_000002655.1').split(',');
const MODES = (process.env.MODES || 'prebuilt,page').split(',');
const PORT = 8000;
const TIMEOUT = 15 * 60_000;

const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
const server = http.createServer((req, res) => {
  const path = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  let file = join(root, path);
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, 'index.html');
  if (!file.startsWith(root) || !existsSync(file)) { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { 'Content-Type': TYPES[extname(file)] || 'application/octet-stream' });
  res.end(readFileSync(file));
});
await new Promise((r) => server.listen(PORT, 'localhost', r));
const base = `http://localhost:${PORT}`;
const browser = await { chromium, firefox, webkit }[ENGINE].launch(
  process.env.BROWSER_PATH ? { executablePath: process.env.BROWSER_PATH } : {});

async function oneRun(accession, mode, checkCounts) {
  const context = await browser.newContext();
  const page = await context.newPage();
  const bytes = {};
  page.on('response', async (r) => {
    try {
      const host = new URL(r.url()).host;
      if (host.startsWith('localhost')) return;
      const len = Number(r.headers()['content-length'] || 0) || (await r.body()).length;
      bytes[host] = (bytes[host] || 0) + len;
    } catch { /* aborted or redirect */ }
  });
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  page.on('console', (m) => { if (m.type() === 'warning' && /Pre-built index/.test(m.text())) errors.push(m.text()); });
  await page.goto(`${base}/analyze/${mode === 'page' ? '?prebuilt=0' : ''}`);
  await page.waitForSelector('#step-genome[data-ready="menu"]', { timeout: 60_000 });
  const t0 = Date.now();
  await page.evaluate((acc) => { window.qrnamapApp.loadCatalogAccession(acc); }, accession);
  await page.waitForSelector('#step-genome[data-phase="ready"], #genome-error:not([hidden])', { timeout: TIMEOUT });
  const seconds = (Date.now() - t0) / 1000;
  const info = await page.evaluate(() => ({
    error: document.getElementById('genome-error').hidden ? '' : document.getElementById('genome-error').textContent,
    readyTime: document.getElementById('ready-time').textContent,
    nTranscripts: window.qrnamapApp.state.nTranscripts,
    nGenes: window.qrnamapApp.state.nGenes,
  }));
  let counts = null;
  if (checkCounts && !info.error) {
    await page.setInputFiles('#sample-table td[data-mate="1"] input', join(FIX, 's1_R1.fastq.gz'));
    await page.setInputFiles('#sample-table td[data-mate="2"] input', join(FIX, 's1_R2.fastq.gz'));
    await page.fill('#sample-table .sample-name', 'S1');
    await page.click('#run-samples');
    await page.waitForFunction(() => ['done', 'cancelled', 'lost'].includes(window.qrnamapApp.state.samplesPhase), null, { timeout: TIMEOUT });
    const s = await page.evaluate(() => window.qrnamapApp.state.samples[0]);
    counts = { status: s.status, error: s.error, processed: s.processed, pseudoaligned: s.pseudoaligned, geneCounts: s.geneCounts };
  }
  await context.close();
  return { seconds, bytes, errors, ...info, counts };
}

const cols = ['engine', 'accession', 'mode', 'run', 'seconds', 'ready_time', 'prebuilt_used', 'n_transcripts', 'n_genes', 's3_bytes', 'ncbi_bytes', 'other_bytes', 'counts_equal', 'error'];
console.log(cols.join('\t'));
const reference = {};
for (const accession of GENOMES) {
  for (let run = 1; run <= RUNS; run++) {
    for (const mode of MODES) {
      const check = process.env.CHECK_COUNTS === '1' && accession === 'GCF_000146045.2' && run === 1;
      process.stderr.write(`${ENGINE} ${accession} ${mode} run ${run}...\n`);
      let r;
      try { r = await oneRun(accession, mode, check); } catch (e) { r = { seconds: NaN, bytes: {}, errors: [], error: e.message.split('\n')[0] }; }
      const s3 = Object.entries(r.bytes).filter(([h]) => /nrp-nautilus/.test(h)).reduce((a, [, n]) => a + n, 0);
      const ncbi = Object.entries(r.bytes).filter(([h]) => /ncbi/.test(h)).reduce((a, [, n]) => a + n, 0);
      const other = Object.values(r.bytes).reduce((a, n) => a + n, 0) - s3 - ncbi;
      let equal = '';
      if (r.counts) {
        reference[mode] = r.counts;
        if (reference.prebuilt && reference.page) {
          equal = String(JSON.stringify(reference.prebuilt.geneCounts) === JSON.stringify(reference.page.geneCounts)
            && reference.prebuilt.pseudoaligned === reference.page.pseudoaligned);
          process.stderr.write(`counts: prebuilt ${reference.prebuilt.pseudoaligned}/${reference.prebuilt.processed} pseudoaligned, page ${reference.page.pseudoaligned}/${reference.page.processed}; equal=${equal}\n`);
        }
        if (r.counts.status !== 'done') r.error = `quant: ${r.counts.error}`;
      }
      const prebuiltUsed = /pre-built/.test(r.readyTime || '');
      console.log([ENGINE, accession, mode, run, r.seconds?.toFixed?.(2), r.readyTime || '', prebuiltUsed, r.nTranscripts ?? '', r.nGenes ?? '',
        s3, ncbi, other, equal, [r.error, ...(r.errors || [])].filter(Boolean).join(' | ').replace(/\s+/g, ' ')].join('\t'));
    }
  }
}
await browser.close();
server.close();
