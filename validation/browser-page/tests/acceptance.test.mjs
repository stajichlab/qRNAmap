// Acceptance test: compare browser testset count tables against native CLI expected references.
// Run with:
//   DIST=<dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/acceptance.test.mjs
import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import http from 'node:http';
import { readFileSync, statSync, existsSync } from 'node:fs';
import { join, resolve, extname, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(process.env.DIST || resolve(here, '../../../dist-site'));
const EXPECTED_DIR = resolve(here, '../expected');
const DATA_ROOT = '/bigdata/stajichlab/jstajich/projects/qRNAmap_data';
const PREVIEW_DIR = join(DATA_ROOT, 'page-v1-preview-reads');
const ENGINE = process.env.ENGINE || 'chromium';
const WAIT_MS = 600_000; // 10 minutes for 6 samples

const TYPES = {
  '.html': 'text/html',
  '.js': 'text/javascript',
  '.mjs': 'text/javascript',
  '.css': 'text/css',
  '.wasm': 'application/wasm',
  '.json': 'application/json',
};

const server = http.createServer((req, res) => {
  const url = new URL(req.url, 'http://x');
  const path = decodeURIComponent(url.pathname);
  let file = resolve(join(root, path));
  if (!file.startsWith(root)) { res.writeHead(403).end(); return; }
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, 'index.html');
  if (!existsSync(file)) { res.writeHead(404).end('not found'); return; }
  res.writeHead(200, {
    'content-type': TYPES[extname(file)] ?? 'application/octet-stream',
    'cache-control': 'no-store',
  }).end(readFileSync(file));
});

let browser;
let base;
let launchSkipped = false;

before(async () => {
  assert.ok(existsSync(join(root, 'analyze/testsets.js')), `dist-site missing at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;

  try {
    browser = await { chromium, firefox, webkit }[ENGINE].launch();
  } catch (err) {
    if (ENGINE === 'webkit') {
      console.warn(`WebKit launch skipped on this host: ${err.message}`);
      launchSkipped = true;
    } else {
      throw err;
    }
  }
});

after(async () => {
  await browser?.close();
  server.close();
});

async function open(t) {
  const context = await browser.newContext({ acceptDownloads: true });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  page.on('console', (msg) => {
    console.log(`[browser ${msg.type()}] ${msg.text()}`);
  });
  t.after(async () => {
    await context.close();
  });
  await page.goto(`${base}/analyze/`);
  await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 30_000 });
  return { page, errors };
}

function diffTables(expectedTsv, actualTsv) {
  const expLines = expectedTsv.trim().split('\n');
  const actLines = actualTsv.trim().split('\n');
  if (expLines[0] !== actLines[0]) {
    return `Header mismatch:\nexpected: "${expLines[0]}"\nactual:   "${actLines[0]}"`;
  }
  const expMap = new Map(expLines.slice(1).map((l) => [l.split('\t')[0], l]));
  const actMap = new Map(actLines.slice(1).map((l) => [l.split('\t')[0], l]));

  const diffs = [];
  for (const [gene, expRow] of expMap) {
    const actRow = actMap.get(gene);
    if (!actRow) diffs.push(`Gene ${gene} missing from actual counts`);
    else if (expRow !== actRow) {
      diffs.push(`Mismatch at gene ${gene}:\nexpected: ${expRow}\nactual:   ${actRow}`);
    }
  }
  for (const [gene] of actMap) {
    if (!expMap.has(gene)) diffs.push(`Extra gene ${gene} in actual counts`);
  }
  if (diffs.length) return diffs.join('\n\n');
  if (expLines.length !== actLines.length) {
    return `Line count mismatch: expected ${expLines.length} lines, got ${actLines.length}`;
  }
  return null;
  if (expLines.length !== actLines.length) {
    return `Line count mismatch: expected ${expLines.length} lines, got ${actLines.length}`;
  }
  return null;
}

test('Acceptance: yeast-ampb counts table matches native CLI reference', async (t) => {
  if (launchSkipped) {
    t.skip(`Skipped on ${ENGINE}: engine could not launch`);
    return;
  }
  const expectedFile = join(EXPECTED_DIR, 'yeast-ampb/counts.tsv');
  assert.ok(existsSync(expectedFile), `Reference counts missing at ${expectedFile}`);
  const expectedCounts = readFileSync(expectedFile, 'utf8');

  const { page, errors } = await open(t);

  // Intercept NCBI catalog genome download with local Yeast genome files
  const yeastFna = join(DATA_ROOT, 'spike0/fixtures/genomes/yeast/GCF_000146045.2_R64_genomic.fna.gz');
  const yeastGff = join(DATA_ROOT, 'spike0/fixtures/genomes/yeast/GCF_000146045.2_R64_genomic.gff.gz');
  await page.route('**/ftp.ncbi.nlm.nih.gov/**', (route) => {
    const url = route.request().url();
    if (url.includes('.fna.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: readFileSync(yeastFna) });
    } else if (url.includes('.gff.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: readFileSync(yeastGff) });
    } else {
      route.continue();
    }
  });

  // Intercept ENA read downloads with local preview reads
  const yeastReadsDir = join(PREVIEW_DIR, 'yeast-ampb');
  await page.route('**/ftp.sra.ebi.ac.uk/**', (route) => {
    const url = route.request().url();
    const filename = url.split('/').pop();
    const localFile = join(yeastReadsDir, filename);
    if (existsSync(localFile)) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: readFileSync(localFile) });
    } else {
      route.continue();
    }
  });

  // Click Load test set: Yeast AmpB
  await page.click('#load-testset-yeast');

  // Verify sample table gets populated with 6 rows
  await page.waitForFunction(() => document.querySelectorAll('#sample-table tbody tr').length === 6, null, { timeout: WAIT_MS });

  // Wait for run to end in results view
  await page.waitForFunction(() => document.getElementById('step-results')?.dataset.state === 'active', null, { timeout: WAIT_MS });

  // Retrieve counts table from page
  const actualCounts = await page.evaluate(() => window.qrnamapResults?.getCountsTsv());
  assert.ok(actualCounts && actualCounts.length > 0, 'actual counts table is non-empty');

  const diff = diffTables(expectedCounts, actualCounts);
  if (diff) {
    assert.fail(`Counts table differs from CLI reference:\n${diff}`);
  }

  // Exact byte comparison
  assert.equal(actualCounts.trim(), expectedCounts.trim(), 'Counts TSV text is identical to CLI reference');
  assert.deepEqual(errors, []);
});

test('Acceptance: cneo-lkh1 counts table matches native CLI reference', async (t) => {
  if (launchSkipped) {
    t.skip(`Skipped on ${ENGINE}: engine could not launch`);
    return;
  }
  const expectedFile = join(EXPECTED_DIR, 'cneo-lkh1/counts.tsv');
  assert.ok(existsSync(expectedFile), `Reference counts missing at ${expectedFile}`);
  const expectedCounts = readFileSync(expectedFile, 'utf8');

  const { page, errors } = await open(t);

  // Intercept NCBI catalog genome download with local H99 genome files
  const h99Fna = join(DATA_ROOT, 'page-v1/dl_h99/GCF_000149245.1_CNA3_genomic.fna.gz');
  const h99Gff = join(DATA_ROOT, 'page-v1/dl_h99/GCF_000149245.1_CNA3_genomic.gff.gz');
  await page.route('**/ftp.ncbi.nlm.nih.gov/**', (route) => {
    const url = route.request().url();
    if (url.includes('.fna.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: readFileSync(h99Fna) });
    } else if (url.includes('.gff.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: readFileSync(h99Gff) });
    } else {
      route.continue();
    }
  });

  // Intercept ENA read downloads with local preview reads
  const cneoReadsDir = join(PREVIEW_DIR, 'cneo-lkh1');
  await page.route('**/ftp.sra.ebi.ac.uk/**', (route) => {
    const url = route.request().url();
    const filename = url.split('/').pop();
    const localFile = join(cneoReadsDir, filename);
    if (existsSync(localFile)) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: readFileSync(localFile) });
    } else {
      route.continue();
    }
  });

  // Click Load test set: C. neoformans H99
  await page.click('#load-testset-cneo');

  // Verify sample table gets populated with 6 rows
  await page.waitForFunction(() => document.querySelectorAll('#sample-table tbody tr').length === 6, null, { timeout: WAIT_MS });

  // Wait for run to end in results view
  await page.waitForFunction(() => document.getElementById('step-results')?.dataset.state === 'active', null, { timeout: WAIT_MS });

  // Retrieve counts table from page
  const actualCounts = await page.evaluate(() => window.qrnamapResults?.getCountsTsv());
  assert.ok(actualCounts && actualCounts.length > 0, 'actual counts table is non-empty');

  const diff = diffTables(expectedCounts, actualCounts);
  if (diff) {
    assert.fail(`Counts table differs from CLI reference:\n${diff}`);
  }

  // Exact byte comparison
  assert.equal(actualCounts.trim(), expectedCounts.trim(), 'Counts TSV text is identical to CLI reference');
  assert.deepEqual(errors, []);
});
