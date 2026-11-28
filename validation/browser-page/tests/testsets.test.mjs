// Browser tests of Task 8: test sets (site/analyze/testsets.js and site/testsets.json).
//   DIST=<absolute path of dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/testsets.test.mjs
import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import http from 'node:http';
import { readFileSync, statSync, existsSync } from 'node:fs';
import { join, resolve, extname, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const here = dirname(fileURLToPath(import.meta.url));
const FIX = resolve(here, '../fixtures');
const root = resolve(process.env.DIST || resolve(here, '../../../dist-site'));
const ENGINE = process.env.ENGINE || 'chromium';
const WAIT_MS = 120_000;

const fixture = (name) => readFileSync(join(FIX, name));

// ------------------------------------------------------------------ static server
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
const server = http.createServer((req, res) => {
  const url = new URL(req.url, 'http://x');
  const path = decodeURIComponent(url.pathname);
  const headers = (type) => ({ 'content-type': type, 'cache-control': 'no-store' });
  let file = resolve(join(root, path));
  if (!file.startsWith(root)) { res.writeHead(403).end(); return; }
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, 'index.html');
  if (!existsSync(file)) { res.writeHead(404).end('not found'); return; }
  res.writeHead(200, headers(TYPES[extname(file)] ?? 'application/octet-stream')).end(readFileSync(file));
});

let browser;
let base;
before(async () => {
  assert.ok(existsSync(join(root, 'analyze/testsets.js')), `no dist-site with testsets.js at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;
  browser = await { chromium, firefox, webkit }[ENGINE].launch();
});
after(async () => {
  await browser?.close();
  server.close();
});

/** The /analyze/ address with `query`, always with `prebuilt=0`: these tests check the in-page build from the NCBI fixtures, and
 * the catalog test genomes have pre-built indexes on S3, which the page would otherwise fetch over the network. */
const analyzeUrl = (query = '') => {
  const params = new URLSearchParams(query);
  if (!params.has('prebuilt')) params.set('prebuilt', '0');
  return `${base}/analyze/?${params}`;
};
async function open(t, query = '') {
  const context = await browser.newContext({ acceptDownloads: true });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  t.after(async () => {
    await context.close();
  });
  await page.goto(analyzeUrl(query));
  await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 30_000 });
  return { page, errors };
}

// ------------------------------------------------------------------ tests


// The test-set reads: the S3 mirror of release testsets-v1 (site/analyze/testsets.js), and ENA when that fails.
const isTestsetReadUrl = (url) => url.hostname === 'ftp.sra.ebi.ac.uk' || url.pathname.includes('/qRNAmap/testsets-v1/');
test('testsets.json exists and defines yeast-ampb and cneo-lkh1 with 3+3 samples and ENA metadata', () => {
  const data = JSON.parse(readFileSync(join(root, 'testsets.json'), 'utf8'));
  assert.equal(data.length, 2);

  const yeast = data.find((d) => d.id === 'yeast-ampb');
  assert.ok(yeast);
  assert.equal(yeast.catalog_accession, 'GCF_000146045.2');
  assert.equal(yeast.samples.length, 6);
  assert.equal(yeast.contrast, 'condition,AmpB,control');
  assert.ok(yeast.selection_reason.length > 20);
  assert.ok(yeast.samples.every((s) => s.sample && s.run && s.fastq1_url && s.fastq2_url && s.read_count > 0));

  const cneo = data.find((d) => d.id === 'cneo-lkh1');
  assert.ok(cneo);
  assert.equal(cneo.catalog_accession, 'GCF_000149245.1');
  assert.equal(cneo.samples.length, 6);
  assert.equal(cneo.contrast, 'genotype,lkh1_delta,wild_type');
  assert.ok(cneo.selection_reason.length > 20);
  assert.ok(cneo.samples.every((s) => s.sample && s.run && s.fastq1_url && s.fastq2_url && s.read_count > 0));
});

test('Load test set button fills genome, fills sample table with 6 rows, fetches reads, and ends in results view', async (t) => {
  const { page, errors } = await open(t);

  // Intercept NCBI catalog genome download with small fixture genome
  await page.route('**/ftp.ncbi.nlm.nih.gov/**', (route) => {
    const url = route.request().url();
    if (url.includes('.fna.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: fixture('mito.fna.gz') });
    } else if (url.includes('.gff.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: fixture('mito.gff3.gz') });
    } else {
      route.continue();
    }
  });

  // Intercept the S3 mirror of the test-set reads and ENA with fixture reads
  await page.route(isTestsetReadUrl, (route) => {
    const url = route.request().url();
    if (url.includes('_1.fastq.gz') || url.includes('_1.fastq')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: fixture('s1_R1.fastq.gz') });
    } else if (url.includes('_2.fastq.gz') || url.includes('_2.fastq')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', headers: { 'access-control-allow-origin': '*' }, body: fixture('s1_R2.fastq.gz') });
    } else {
      route.continue();
    }
  });

  // Click Load test set: Yeast AmpB
  await page.click('#load-testset-yeast');

  // Verify sample table gets populated with 6 rows
  await page.waitForFunction(() => document.querySelectorAll('#sample-table tbody tr').length === 6, null, { timeout: 30_000 });

  // Wait for run to end in results view
  await page.waitForFunction(() => document.getElementById('step-results')?.dataset.state === 'active', null, { timeout: WAIT_MS });

  // Results view checks
  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'active');
  const mappingCount = await page.locator('#results-mapping-table tbody tr').count();
  assert.equal(mappingCount, 6, 'mapping summary lists all 6 test set samples');

  assert.deepEqual(errors, []);
});

test('ENA failure (HTTP 404 or network error) gives clear message in testset-errors', async (t) => {
  const { page, errors } = await open(t);

  // Route NCBI ok
  await page.route('**/ftp.ncbi.nlm.nih.gov/**', (route) => {
    const url = route.request().url();
    if (url.includes('.fna.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', body: fixture('mito.fna.gz') });
    } else if (url.includes('.gff.gz')) {
      route.fulfill({ status: 200, contentType: 'application/gzip', body: fixture('mito.gff3.gz') });
    } else {
      route.continue();
    }
  });

  // Route ENA and the S3 test-set reads with 404
  await page.route(isTestsetReadUrl, (route) => {
    route.fulfill({ status: 404, statusText: 'Not Found', body: 'Not found' });
  });

  // Click Load test set
  await page.click('#load-testset-yeast');

  // Wait for error to appear
  await page.waitForSelector('#testset-errors:not([hidden])', { timeout: 30_000 });
  const errorText = await page.textContent('#testset-errors');
  assert.match(errorText, /ENA|HTTP 404|failed/i);

  assert.deepEqual(errors, []);
});
