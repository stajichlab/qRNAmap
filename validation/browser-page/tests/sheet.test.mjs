// Browser tests of the sample sheet builder (site/sheet/).
//   DIST=<absolute path of dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/sheet.test.mjs
// SCREENSHOT_DIR=<dir> also saves a screenshot of the filled page.
import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import http from 'node:http';
import { readFileSync, statSync, existsSync } from 'node:fs';
import { join, resolve, extname, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const here = dirname(fileURLToPath(import.meta.url));
const root = resolve(process.env.DIST || resolve(here, '../../../dist-site'));
const ENGINE = process.env.ENGINE || 'chromium';

const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
const server = http.createServer((req, res) => {
  const path = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  let file = resolve(join(root, path));
  if (!file.startsWith(root)) { res.writeHead(403).end(); return; }
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, 'index.html');
  if (!existsSync(file)) { res.writeHead(404).end('not found'); return; }
  res.writeHead(200, { 'content-type': TYPES[extname(file)] ?? 'application/octet-stream', 'cache-control': 'no-store' }).end(readFileSync(file));
});

let browser;
let base;
before(async () => {
  assert.ok(existsSync(join(root, 'sheet/index.html')), `no dist-site/sheet at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;
  browser = await { chromium, firefox, webkit }[ENGINE].launch();
});
after(async () => {
  await browser?.close();
  server.close();
});

async function open(t, context) {
  const ctx = context || await browser.newContext({ acceptDownloads: true });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  if (!context) t.after(() => ctx.close());
  await page.goto(`${base}/sheet/`);
  await page.waitForSelector('body[data-ready="true"]', { timeout: 30_000 });
  return { page, errors, ctx };
}

const NAMES = [
  'WT_rep1_S1_L001_R1_001.fastq.gz', 'WT_rep1_S1_L001_R2_001.fastq.gz',
  'WT_rep2_S2_L001_R1_001.fastq.gz', 'WT_rep2_S2_L001_R2_001.fastq.gz',
  'dcat_rep1_S3_L001_R1_001.fastq.gz', 'dcat_rep1_S3_L001_R2_001.fastq.gz',
  'dcat_rep1_S3_L002_R1_001.fastq.gz', 'dcat_rep1_S3_L002_R2_001.fastq.gz',
  'dcat_rep2_S4_L001_I1_001.fastq.gz', 'md5sum.txt',
  'lone_a.fq.gz', 'lone_b.fq.gz',
];

const rowNames = (page) => page.$$eval('#builder-table .sample-name', (xs) => xs.map((x) => x.value));
/** The table row of a sample (by the name input's current value). */
const rowOf = async (page, name) => {
  const i = (await rowNames(page)).indexOf(name);
  assert.ok(i >= 0, `no row "${name}"`);
  return page.locator('#builder-table tbody tr').nth(i);
};
const select = async (page, names) => {
  for (const n of names) await (await rowOf(page, n)).locator('.row-select').check();
};

test('pasted names become paired samples; skips listed; treatments set; the sheet validates and downloads', async (t) => {
  const { page, errors } = await open(t);
  await page.click('#paste-box summary');
  await page.fill('#paste-names', NAMES.join('\n'));
  await page.click('#btn-use-pasted');
  assert.deepEqual(await rowNames(page), ['dcat_rep1', 'lone_a', 'lone_b', 'WT_rep1', 'WT_rep2']);
  assert.match(await page.textContent('#group-summary'), /5 sample\(s\): 3 paired-end, 2 single-end\./);
  assert.match(await page.textContent('#skipped-summary'), /2 file\(s\) left out/);
  assert.match(await page.textContent('#skipped-list'), /index read/);
  assert.match(await page.textContent('#page-notes'), /dcat_rep1/); // two lanes: the Analyze page needs them joined

  // The default "condition" column is empty: the check names the empty cells.
  assert.equal(await page.isVisible('#sheet-check'), true);
  assert.match(await page.textContent('#sheet-check'), /condition/);

  // Pair the two single-end rows by hand.
  await select(page, ['lone_a', 'lone_b']);
  await page.click('#btn-pair');
  assert.deepEqual(await rowNames(page), ['dcat_rep1', 'lone_a', 'WT_rep1', 'WT_rep2']);
  const lone = (await rowOf(page, 'lone_a'));
  assert.match(await lone.locator('td.files').nth(1).textContent(), /lone_b\.fq\.gz/);

  // Columns from the name parts, then rename and bulk-set.
  await page.click('#btn-name-parts');
  const heads = await page.$$eval('#builder-table .col-name', (xs) => xs.map((x) => x.value));
  assert.deepEqual(heads, ['condition', 'part1', 'part2']);
  await select(page, ['WT_rep1', 'WT_rep2']);
  await page.selectOption('#bulk-col', { label: 'condition' });
  await page.fill('#bulk-value', 'control');
  await page.click('#btn-bulk-set');
  await page.locator('#builder-table .row-select').first().check(); // dcat_rep1
  await (await rowOf(page, 'WT_rep1')).locator('.row-select').uncheck();
  await (await rowOf(page, 'WT_rep2')).locator('.row-select').uncheck();
  await select(page, ['lone_a']);
  await page.fill('#bulk-value', 'treated');
  await page.click('#btn-bulk-set');
  // The part columns are not wanted: remove part1, rename part2 to replicate.
  await page.click('#builder-table button.remove-col[aria-label="Remove column part1"]');
  await page.fill('#builder-table .col-name >> nth=1', 'replicate');
  await (await rowOf(page, 'lone_a')).locator('.cov').nth(1).fill('rep1');
  assert.equal(await page.isVisible('#sheet-check'), false, await page.textContent('#sheet-check'));
  assert.match(await page.textContent('#sheet-ok'), /The sheet is valid: 4 samples, 2 covariate column/);

  await page.fill('#path-prefix', '/bigdata/lab/run1');
  if (process.env.SCREENSHOT_DIR) await page.screenshot({ path: join(process.env.SCREENSHOT_DIR, `sheet-${ENGINE}.png`), fullPage: true });
  const [download] = await Promise.all([page.waitForEvent('download'), page.click('#btn-download')]);
  assert.equal(download.suggestedFilename(), 'samples.tsv');
  const text = readFileSync(await download.path(), 'utf8');
  const lines = text.trimEnd().split('\n');
  assert.equal(lines[0], 'sample\tfastq1\tfastq2\tcondition\treplicate');
  assert.equal(lines[1], 'dcat_rep1\t/bigdata/lab/run1/dcat_rep1_S3_L001_R1_001.fastq.gz,/bigdata/lab/run1/dcat_rep1_S3_L002_R1_001.fastq.gz\t' +
    '/bigdata/lab/run1/dcat_rep1_S3_L001_R2_001.fastq.gz,/bigdata/lab/run1/dcat_rep1_S3_L002_R2_001.fastq.gz\ttreated\trep1');
  assert.equal(lines[2], 'lone_a\t/bigdata/lab/run1/lone_a.fq.gz\t/bigdata/lab/run1/lone_b.fq.gz\ttreated\trep1');
  assert.equal(lines[3], 'WT_rep1\t/bigdata/lab/run1/WT_rep1_S1_L001_R1_001.fastq.gz\t/bigdata/lab/run1/WT_rep1_S1_L001_R2_001.fastq.gz\tcontrol\trep1');
  assert.deepEqual(errors, []);
});

test('picked files give names only; _1/_2 mates can be turned off; split, swap; the draft survives a reload', async (t) => {
  const context = await browser.newContext({ acceptDownloads: true });
  t.after(() => context.close());
  const { page, errors } = await open(t, context);
  const files = ['SRR1_1.fastq.gz', 'SRR1_2.fastq.gz', 'SRR2_1.fastq.gz', 'SRR2_2.fastq.gz']
    .map((name) => ({ name, mimeType: 'application/gzip', buffer: Buffer.from('not read') }));
  await page.setInputFiles('#pick-files', files);
  assert.deepEqual(await rowNames(page), ['SRR1', 'SRR2']);
  assert.match(await page.textContent('#names-status'), /4 name\(s\) added/);

  await page.uncheck('#numeric-mates');
  assert.deepEqual(await rowNames(page), ['SRR1_1', 'SRR1_2', 'SRR2_1', 'SRR2_2']);
  await page.check('#numeric-mates');
  assert.deepEqual(await rowNames(page), ['SRR1', 'SRR2']);

  await select(page, ['SRR2']);
  await page.click('#btn-swap');
  const srr2 = (await rowOf(page, 'SRR2'));
  assert.match(await srr2.locator('td.files').first().textContent(), /SRR2_2/);
  await select(page, ['SRR2']);
  await page.click('#btn-split');
  assert.deepEqual(await rowNames(page), ['SRR1', 'SRR2-1', 'SRR2-2']);
  await page.click('#btn-swap');
  assert.match(await page.textContent('#action-error'), /Select one or more paired rows/);

  await page.locator('#builder-table .cov').first().fill('ctl');
  await page.reload();
  await page.waitForSelector('body[data-ready="true"]');
  assert.deepEqual(await rowNames(page), ['SRR1', 'SRR2-1', 'SRR2-2']);
  assert.equal(await page.locator('#builder-table .cov').first().inputValue(), 'ctl');
  assert.match(await page.textContent('#names-status'), /from your last visit/);

  await page.click('#btn-start-over');
  assert.equal(await page.isVisible('#sheet-table-step'), false);
  assert.deepEqual(errors, []);
});

test('the landing page and the Analyze page link to the builder', async (t) => {
  const ctx = await browser.newContext();
  t.after(() => ctx.close());
  const page = await ctx.newPage();
  await page.goto(`${base}/`);
  assert.equal(await page.getAttribute('#sheet-link', 'href'), 'sheet/');
  await page.goto(`${base}/analyze/`);
  assert.equal(await page.locator('#sheet-box a[href="../sheet/"]').count(), 1);
});
