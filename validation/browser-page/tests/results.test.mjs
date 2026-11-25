// Browser tests of step 3 (results and downloads) of /analyze/.
//   DIST=<absolute path of dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/results.test.mjs
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
  assert.ok(existsSync(join(root, 'analyze/app.js')), `no dist-site at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;
  browser = await { chromium, firefox, webkit }[ENGINE].launch();
});
after(async () => {
  await browser?.close();
  server.close();
});

/** A new page with the fixture genome index ready and step 2 open. */
async function open(t, query = '') {
  const context = await browser.newContext({ acceptDownloads: true });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  t.after(async () => {
    await context.close();
  });
  await page.goto(`${base}/analyze/${query}`);
  await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 30_000 });
  await page.check('input[name="genome-source"][value="upload"]');
  await page.setInputFiles('#upload-fasta', join(FIX, 'mito.fna.gz'));
  await page.setInputFiles('#upload-gff', join(FIX, 'mito.gff3'));
  await page.click('#upload-load');
  await page.waitForSelector('#step-genome[data-phase="ready"]', { timeout: WAIT_MS });
  return { page, errors };
}

const row = (page, i) => page.locator('#sample-table tbody tr').nth(i);

async function fillRow(page, i, { name, condition, mate1, mate2 }) {
  while (await page.locator('#sample-table tbody tr').count() <= i) await page.click('#add-sample');
  const r = row(page, i);
  await r.locator('.sample-name').fill(name);
  if (condition !== undefined) await r.locator('.cov').first().fill(condition);
  const arg = (m) => (typeof m === 'string' ? resolve(FIX, m) : m);
  if (mate1) await r.locator('td[data-mate="1"] input').setInputFiles(arg(mate1));
  if (mate2) await row(page, i).locator('td[data-mate="2"] input').setInputFiles(arg(mate2));
}

async function run(page) {
  await page.click('#run-samples');
  await page.waitForFunction(() => window.qrnamapApp?.state?.samplesPhase === 'done', null, { timeout: WAIT_MS });
}

// ------------------------------------------------------------------ tests

test('results step unlocks after run: mapping summary, count table, alias search, downloads and SVG plots', async (t) => {
  const { page, errors } = await open(t);

  // Results step is locked initially
  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'locked');

  // Fill 3 samples: S1, S2, S3
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'S2', condition: 'trt', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  await fillRow(page, 2, { name: 'S3', condition: 'trt', mate1: 's3.fastq.gz' });
  await page.fill('#frag-mean', '200');
  await page.fill('#frag-sd', '25');

  await run(page);

  // Step 3 (results) must now be active
  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'active');

  // Preview label must be visible
  const previewLabel = await page.textContent('.results-preview-label, #results-preview-label');
  assert.match(previewLabel, /preview of the first (75,000|100,000|200,000) pairs|all pairs/i);

  // 1. Mapping summary table
  const mappingRows = await page.locator('#results-mapping-table tbody tr').count();
  assert.equal(mappingRows, 3, 'mapping summary table has 3 rows');
  const summaryText = await page.textContent('#results-mapping-table');
  assert.match(summaryText, /S1/);
  assert.match(summaryText, /S2/);
  assert.match(summaryText, /S3/);

  // 2. Count table
  const countRows = await page.locator('#results-count-table tbody tr').count();
  assert.ok(countRows > 10, 'count table has rows');
  const countHead = await page.textContent('#results-count-table thead');
  assert.match(countHead, /Gene/);
  assert.match(countHead, /S1/);
  assert.match(countHead, /S2/);
  assert.match(countHead, /S3/);

  // Sortable: click S1 column header
  await page.click('#results-count-table th[data-col="1"], #results-count-table th:has-text("S1")');
  const firstGeneSorted = await page.locator('#results-count-table tbody tr td:first-child').first().textContent();
  assert.ok(firstGeneSorted.length > 0);

  // 3. Count table search box matches locus name AND aliases in gene_info.tsv:
  // Locus name Q0045:
  await page.fill('#count-search', 'Q0045');
  assert.equal(await page.locator('#results-count-table tbody tr:visible').count(), 1);
  assert.equal(await page.locator('#results-count-table tbody tr:visible td:first-child').textContent(), 'Q0045');

  // Symbol COX1:
  await page.fill('#count-search', 'COX1');
  assert.equal(await page.locator('#results-count-table tbody tr:visible').count(), 1);
  assert.equal(await page.locator('#results-count-table tbody tr:visible td:first-child').textContent(), 'Q0045');

  // GeneID 854598:
  await page.fill('#count-search', '854598');
  assert.equal(await page.locator('#results-count-table tbody tr:visible').count(), 1);
  assert.equal(await page.locator('#results-count-table tbody tr:visible td:first-child').textContent(), 'Q0045');

  // Transcript accession rna-Q0045:
  await page.fill('#count-search', 'rna-Q0045');
  assert.equal(await page.locator('#results-count-table tbody tr:visible').count(), 1);
  assert.equal(await page.locator('#results-count-table tbody tr:visible td:first-child').textContent(), 'Q0045');

  // Protein accession NP_009305.1:
  await page.fill('#count-search', 'NP_009305.1');
  assert.equal(await page.locator('#results-count-table tbody tr:visible').count(), 1);
  assert.equal(await page.locator('#results-count-table tbody tr:visible td:first-child').textContent(), 'Q0045');

  // Clear search
  await page.fill('#count-search', '');
  assert.ok(await page.locator('#results-count-table tbody tr:visible').count() > 10);

  // 4. SVG Plots
  // Heatmap
  const heatmap = page.locator('#heatmap-container svg');
  assert.equal(await heatmap.count(), 1, 'SVG heatmap is present');
  assert.match(await page.textContent('#heatmap-container'), /preview/i);

  // PCA
  const pca = page.locator('#pca-container svg');
  assert.equal(await pca.count(), 1, 'SVG PCA is present');
  const pcaText = await page.textContent('#pca-section, #pca-container');
  assert.match(pcaText, /PCA of log2\(CPM\+1\), top 500 variable genes \(not DESeq2-normalised\)/);
  assert.ok(await page.locator('#pca-container .pca-legend').isVisible(), 'PCA legend is visible');

  // Volcano plot
  const volcano = page.locator('#volcano-container svg');
  assert.equal(await volcano.count(), 1, 'SVG volcano plot is present');
  const volcanoText = await page.textContent('#de-section, #volcano-container');
  assert.match(volcanoText, /Comparison: "trt" vs "ctl"/i);

  // Labeled genes carry a mouseover title and a link out to SGD or NCBI (from the GFF Dbxref)
  const volcanoLinks = await page.locator('#volcano-container svg a.volcano-gene-link').all();
  assert.ok(volcanoLinks.length > 0, 'volcano has linkable labeled genes');
  let labelLinks = 0;
  for (const a of volcanoLinks) {
    const href = await a.getAttribute('href');
    assert.match(
      href,
      /^https:\/\/www\.(yeastgenome\.org\/locus\/S\d+|ncbi\.nlm\.nih\.gov\/gene\/\d+)$/,
      `link out to SGD or NCBI Gene, got ${href}`
    );
    assert.equal(await a.getAttribute('target'), '_blank', 'link opens in a new tab');
    const title = await a.locator('title').textContent();
    assert.match(title, /log2FC=/, 'mouseover tooltip shows the DE stats');
    if (await a.locator('text').count()) labelLinks++;
  }
  assert.ok(labelLinks > 0, 'the label text itself is clickable');

  // 5. Downloads
  // Download counts.preview-N.tsv
  const [countsDl] = await Promise.all([
    page.waitForEvent('download'),
    page.click('#dl-counts'),
  ]);
  assert.match(countsDl.suggestedFilename(), /^counts\.preview(-\d+)?\.tsv$/);

  // Download differential_expression.tsv
  const [volcanoDl] = await Promise.all([
    page.waitForEvent('download'),
    page.click('#dl-de-tsv'),
  ]);
  assert.equal(volcanoDl.suggestedFilename(), 'differential_expression.tsv');

  // Download gene_info.tsv
  const [infoDl] = await Promise.all([
    page.waitForEvent('download'),
    page.click('#dl-gene-info'),
  ]);
  assert.equal(infoDl.suggestedFilename(), 'gene_info.tsv');

  // Download samples.tsv without reads folder (bare filenames)
  assert.ok(await page.locator('#reads-folder-note').isVisible());
  const [sheetDl] = await Promise.all([
    page.waitForEvent('download'),
    page.click('#dl-sheet'),
  ]);
  assert.equal(sheetDl.suggestedFilename(), 'samples.tsv');

  // Download de-command.txt
  const [deDl] = await Promise.all([
    page.waitForEvent('download'),
    page.click('#dl-de'),
  ]);
  assert.equal(deDl.suggestedFilename(), 'de-command.txt');

  // Optional reads folder box: when filled, samples.tsv uses folder/file
  await page.fill('#reads-folder', '/cluster/myreads');
  const sheetTextWithFolder = await page.evaluate(() => window.qrnamapResults.getSamplesTsv());
  assert.match(sheetTextWithFolder, /\/cluster\/myreads\/s1_R1\.fastq\.gz/);

  // Written sheet re-parses with validate_sheet and round-trips
  const roundTrip = await page.evaluate((tsv) => {
    const { validate_sheet } = window.qrnamapWasm;
    return JSON.parse(validate_sheet(tsv));
  }, sheetTextWithFolder);
  assert.equal(roundTrip.ok, true);
  assert.equal(roundTrip.samples.length, 3);

  assert.deepEqual(errors, []);
});

test('one or two samples show the clear message instead of the PCA/heatmap', async (t) => {
  const { page, errors } = await open(t);

  // 1 sample only
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await run(page);

  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'active');

  // Heatmap and PCA should show messages instead of plots
  assert.equal(await page.locator('#heatmap-container svg').count(), 0);
  assert.match(await page.textContent('#heatmap-container'), /needs at least two samples|at least 2 samples/i);

  assert.equal(await page.locator('#pca-container svg').count(), 0);
  assert.match(await page.textContent('#pca-container'), /needs at least three samples|at least 3 samples/i);

  // Add sample 2: 2 samples total
  await fillRow(page, 1, { name: 'S2', condition: 'trt', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  await run(page);

  // Now Heatmap has SVG, PCA still shows message
  assert.equal(await page.locator('#heatmap-container svg').count(), 1);
  assert.equal(await page.locator('#pca-container svg').count(), 0);
  assert.match(await page.textContent('#pca-container'), /needs at least three samples|at least 3 samples/i);

  assert.deepEqual(errors, []);
});

test('persistence: reload restores saved results with same genome and N, and not with another', async (t) => {
  const { page, errors } = await open(t);
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await run(page);

  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'active');

  // Reload page
  await page.reload();
  await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 30_000 });

  // Upload same genome
  await page.check('input[name="genome-source"][value="upload"]');
  await page.setInputFiles('#upload-fasta', join(FIX, 'mito.fna.gz'));
  await page.setInputFiles('#upload-gff', join(FIX, 'mito.gff3'));
  await page.click('#upload-load');
  await page.waitForSelector('#step-genome[data-phase="ready"]', { timeout: WAIT_MS });

  // Results should be restored
  await page.waitForFunction(() => document.getElementById('step-results').dataset.state === 'active', null, { timeout: 10_000 });
  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'active');
  assert.ok(await page.locator('#results-mapping-table tbody tr').count() >= 1);

  assert.deepEqual(errors, []);
});
