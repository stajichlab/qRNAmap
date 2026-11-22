// Browser tests of step 1 (genome) of /analyze/.
//   DIST=<absolute path of dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/genome.test.mjs
// A local static server serves dist-site. NCBI is never contacted: page.route answers the catalog URLs of
// S. cerevisiae S288C (GCF_000146045.2) with the small fixtures in ../fixtures (the yeast mitochondrial genome).
// The expected gene and transcript counts come from ../fixtures/expected.tsv (made with `qrnamap transcripts`).
// REAL_NCBI=1 adds one test that downloads S288C from NCBI through the page (not for CI).
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
const YEAST = 'GCF_000146045.2';
const NCBI_DIR = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/';
const FNA_URL = `${NCBI_DIR}GCF_000146045.2_R64_genomic.fna.gz`;
const GFF_URL = `${NCBI_DIR}GCF_000146045.2_R64_genomic.gff.gz`;
const READY_MS = 120_000;

const fixture = (name) => readFileSync(join(FIX, name));
const [, expectedRow] = readFileSync(join(FIX, 'expected.tsv'), 'utf8').trim().split('\n');
const [, expGenes, expTranscripts] = expectedRow.split('\t');
const EXPECTED = { genes: Number(expGenes), transcripts: Number(expTranscripts) };

// ------------------------------------------------------------------ static server, with two test switches for kallisto.mjs
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
const hooks = { kallisto: 'serve', release: null, heldAborted: false };
const server = http.createServer(async (req, res) => {
  const path = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  const headers = (type) => ({ 'content-type': type, 'cache-control': 'no-store' });
  if (path === '/vendor/kallisto/kallisto.mjs' && hooks.kallisto === 'oom') {
    // A kallisto module that fails as an out-of-memory WebAssembly build fails.
    res.writeHead(200, headers('text/javascript'))
      .end("export default async function () { throw new RangeError('WebAssembly.Memory(): could not allocate memory'); }\n");
    return;
  }
  if (path === '/vendor/kallisto/kallisto.mjs' && hooks.kallisto === 'hold') {
    hooks.heldAborted = false;
    // A diagnostic: the browser closes this connection when the worker that asked for the module goes away.
    res.on('close', () => { if (!res.writableEnded) hooks.heldAborted = true; });
    await new Promise((r) => { hooks.release = r; });
  }
  let file = resolve(join(root, path));
  if (!file.startsWith(root)) { res.writeHead(403).end(); return; }
  if (existsSync(file) && statSync(file).isDirectory()) file = join(file, 'index.html');
  if (!existsSync(file)) { res.writeHead(404).end('not found'); return; }
  res.writeHead(200, headers(TYPES[extname(file)] ?? 'application/octet-stream')).end(readFileSync(file));
});

let browser;
let base;
before(async () => {
  assert.ok(existsSync(join(root, 'analyze/index.html')), `no dist-site at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;
  browser = await { chromium, firefox, webkit }[ENGINE].launch();
});
after(async () => {
  await browser?.close();
  server.close();
});

/** A new context and page. `ncbi(route, url)` answers NCBI requests (default: the fixtures for the S288C URLs, 404 otherwise). */
async function open(t, { query = '', ncbi, context } = {}) {
  context = context || await browser.newContext();
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  const ncbiRequests = [];
  await page.route('https://ftp.ncbi.nlm.nih.gov/**', async (route) => {
    const url = route.request().url();
    ncbiRequests.push(url);
    if (ncbi) return ncbi(route, url);
    if (url === FNA_URL) return route.fulfill({ status: 200, body: fixture('mito.fna.gz'), headers: { 'access-control-allow-origin': '*' } });
    if (url === GFF_URL) return route.fulfill({ status: 200, body: fixture('mito.gff3.gz'), headers: { 'access-control-allow-origin': '*' } });
    return route.fulfill({ status: 404, body: 'not found', headers: { 'access-control-allow-origin': '*' } });
  });
  await page.goto(`${base}/analyze/${query}`);
  await page.waitForSelector('#step-genome[data-ready="menu"]', { timeout: 30_000 });
  t.after(async () => {
    hooks.kallisto = 'serve';
    hooks.release?.();
    hooks.release = null;
    await context.close();
  });
  return { page, errors, ncbiRequests, context };
}

async function chooseYeast(page) {
  await page.fill('#genome-filter', YEAST);
  await page.selectOption('#genome-menu', YEAST);
  await page.click('#genome-load');
}

async function upload(page, fasta, gff) {
  await page.check('input[name="genome-source"][value="upload"]');
  await page.setInputFiles('#upload-fasta', fasta);
  await page.setInputFiles('#upload-gff', gff);
  await page.click('#upload-load');
}

const phase = (page) => page.getAttribute('#step-genome', 'data-phase');
const waitPhase = (page, p, timeout = READY_MS) => page.waitForSelector(`#step-genome[data-phase="${p}"]`, { timeout });
const appState = (page) => page.evaluate(() => {
  const { genome, nTranscripts, nGenes, indexReady, phase } = window.qrnamapApp.state;
  return { genome, nTranscripts, nGenes, indexReady, phase };
});

/** Wait until the step is back at "choose genome" with an error, and return the error text. */
async function errorText(page) {
  await page.waitForSelector('#genome-error:not([hidden])', { timeout: READY_MS });
  assert.equal(await phase(page), 'choose');
  const s = await appState(page);
  assert.equal(s.indexReady, false);
  assert.equal(s.genome, null);
  assert.equal(await page.getAttribute('#step-samples', 'data-state'), 'locked');
  return page.textContent('#genome-error');
}

const fileArg = (name, buffer) => ({ name, mimeType: 'application/octet-stream', buffer });

// ------------------------------------------------------------------ tests
test('the menu lists the catalog, RefSeq first, and the search filters it', async (t) => {
  const { page, errors } = await open(t);
  const groups = await page.$$eval('#genome-menu optgroup', (gs) => gs.map((g) => `${g.label}:${g.children.length}`));
  assert.deepEqual(groups, ['NCBI RefSeq:153', 'NCBI GenBank:1']);
  await page.fill('#genome-filter', 'cerevisiae');
  const options = await page.$$eval('#genome-menu option', (os) => os.map((o) => o.textContent));
  assert.deepEqual(options, ['Saccharomyces cerevisiae S288C | GCF_000146045.2 | 12.1 Mb | 6,477 genes']);
  assert.match(await page.textContent('#construction-notice'), /Under construction: the genome step and the quantification work; the results view is not built yet\. Tested: the genome step in headless Chromium, Firefox and WebKit on a GitHub runner \(small fixture genome\); the quantification in headless Chromium so far/);
  assert.deepEqual(errors, []);
});

test('catalog genome: download, transcripts and index reach "index ready" with the CLI counts', async (t) => {
  const { page, errors, context } = await open(t);
  await chooseYeast(page);
  await waitPhase(page, 'ready');
  assert.equal(await page.textContent('#ready-label'), 'Saccharomyces cerevisiae S288C');
  assert.equal(await page.textContent('#ready-accession'), YEAST);
  assert.equal(await page.textContent('#ready-genes'), String(EXPECTED.genes));
  assert.equal(await page.textContent('#ready-transcripts'), String(EXPECTED.transcripts));
  const s = await appState(page);
  assert.equal(s.indexReady, true);
  assert.equal(s.nTranscripts, EXPECTED.transcripts);
  assert.equal(s.nGenes, EXPECTED.genes);
  assert.equal(s.genome.accession, YEAST);
  // Step 2 (samples) opens; step 3 (results) is not built and stays locked.
  assert.equal(await page.getAttribute('#step-genome', 'data-state'), 'done');
  assert.equal(await page.getAttribute('#step-samples', 'data-state'), 'active');
  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'locked');
  // The index argv is the one of the Task 3 measurements (qrnamap-core toolcall, no -t).
  assert.deepEqual(await page.evaluate(() => window.qrnamapApp.state.indexArgs), ['index', '-i', 'transcripts.idx', 'transcripts.fa']);
  // The fixture sizes differ from the catalog sizes: a soft note in the state, not an error.
  assert.equal((await page.evaluate(() => window.qrnamapApp.state.downloadNotes)).length, 2);
  assert.deepEqual(errors, []);
  // The last genome is remembered (persist.js) and selected again after a reload; it is not downloaded again by itself.
  await page.waitForSelector('#step-genome[data-saved="catalog"]', { timeout: 10_000 });
  const again = await open(t, { context });
  assert.equal(await again.page.inputValue('#genome-menu'), YEAST);
  assert.match(await again.page.textContent('#genome-note'), /Last genome: Saccharomyces cerevisiae S288C \(GCF_000146045\.2\)/);
  assert.equal(await phase(again.page), 'choose');
  assert.deepEqual(again.ncbiRequests, []);
});

test('upload FASTA (gzip) and GFF3 (plain) reach "index ready"', async (t) => {
  const { page, errors } = await open(t);
  await upload(page, join(FIX, 'mito.fna.gz'), join(FIX, 'mito.gff3'));
  await waitPhase(page, 'ready');
  assert.equal(await page.textContent('#ready-transcripts'), String(EXPECTED.transcripts));
  assert.equal(await page.textContent('#ready-genes'), String(EXPECTED.genes));
  assert.equal(await page.textContent('#ready-accession'), 'uploaded files');
  assert.deepEqual(errors, []);
});

test('FASTA names that do not match the GFF3 give the build_transcriptome message', async (t) => {
  const { page } = await open(t);
  const gff = fixture('mito.gff3').toString().replaceAll('NC_001224.1\t', 'chrM\t');
  await upload(page, join(FIX, 'mito.fna'), fileArg('renamed.gff3', Buffer.from(gff)));
  const text = await errorText(page);
  assert.match(text, /The transcripts could not be built: gene '.+' is on contig 'chrM', which is not in the genome FASTA\. Check that the FASTA sequence names match the first column of the GFF3\./);
});

test('a GFF3 with no usable transcripts gives a clear message', async (t) => {
  const { page } = await open(t);
  const gff = fixture('mito.gff3').toString().split('\n').filter((l) => l.startsWith('#') || l.split('\t')[2] === 'region').join('\n') + '\n';
  await upload(page, join(FIX, 'mito.fna'), fileArg('regions.gff3', Buffer.from(gff)));
  assert.match(await errorText(page), /the GFF3 gives no usable transcripts \(no gene has an exon, CDS or UTR interval\)/);
});

test('a truncated gzip download gives the gzip message', async (t) => {
  const gz = fixture('mito.fna.gz');
  const { page } = await open(t, {
    ncbi: (route, url) => route.fulfill({ status: 200, headers: { 'access-control-allow-origin': '*' },
      body: url === FNA_URL ? gz.subarray(0, Math.floor(gz.length / 2)) : fixture('mito.gff3.gz') }),
  });
  await chooseYeast(page);
  assert.match(await errorText(page), /The genome FASTA \(GCF_000146045\.2_R64_genomic\.fna\.gz\) is not a complete gzip file \(unexpected EOF\)\. The file is damaged or cut short\. Download it again\./);
});

test('an empty file gives a clear message', async (t) => {
  const { page } = await open(t);
  await upload(page, fileArg('empty.fna', Buffer.alloc(0)), join(FIX, 'mito.gff3'));
  assert.match(await errorText(page), /The genome FASTA \(empty\.fna\) is empty \(0 bytes\)\./);
});

test('a catalog genome above the size limit is refused before download and names the command-line route', async (t) => {
  const { page, ncbiRequests } = await open(t, { query: '?genome_limit_bases=1000000' });
  await chooseYeast(page);
  const text = await errorText(page);
  assert.match(text, /has 12,071,326 bases, more than the page limit of 1,000,000 bases/);
  assert.match(text, /The largest genome tested in a browser is Fulvia fulva Race5_Kim \(GCF_020509005\.1, 67,169,167 bases/);
  assert.match(text, /Use the command-line tool instead .*qrnamap fetch --source ncbi --id GCF_000146045\.2 --out genome\//);
  assert.deepEqual(ncbiRequests, []);
});

test('an uploaded genome above the size limit is refused and names the command-line route', async (t) => {
  const { page } = await open(t, { query: '?genome_limit_bases=50000' });
  await upload(page, join(FIX, 'mito.fna.gz'), join(FIX, 'mito.gff3'));
  const text = await errorText(page);
  assert.match(text, /has more than 50,000 bases/);
  assert.match(text, /Use the command-line tool instead .*qrnamap transcripts --reference <your FASTA> --gff <your GFF3>/);
});

test('cancel during the download returns to "choose genome", and a new run works', async (t) => {
  let held = null;
  const { page, errors } = await open(t, {
    ncbi: (route, url) => {
      if (url === FNA_URL && !held) { held = route; return undefined; } // never answered
      return route.fulfill({ status: 200, headers: { 'access-control-allow-origin': '*' }, body: fixture(url === FNA_URL ? 'mito.fna.gz' : 'mito.gff3.gz') });
    },
  });
  await chooseYeast(page);
  await waitPhase(page, 'download', 10_000);
  await page.click('#genome-cancel');
  assert.equal(await phase(page), 'choose');
  assert.equal(await page.isHidden('#genome-busy'), true);
  assert.equal(await page.isHidden('#genome-error'), true);
  assert.match(await page.textContent('#genome-note'), /Cancelled/);
  assert.equal((await appState(page)).indexReady, false);
  await held?.fulfill({ status: 200, body: fixture('mito.fna.gz') }).catch(() => {}); // a late answer must change nothing
  await page.waitForTimeout(500);
  assert.equal(await phase(page), 'choose');
  await chooseYeast(page);
  await waitPhase(page, 'ready');
  assert.equal(await page.textContent('#ready-transcripts'), String(EXPECTED.transcripts));
  assert.deepEqual(errors, []);
});

/** A test-only wrapper of window.Worker: it records every construction and every terminate() call in window.__workerLog. */
function recordWorkers() {
  window.__workerLog = [];
  window.__replaced = [];
  window.addEventListener('qrnamap:worker-replaced', (e) => window.__replaced.push(e.detail.reason));
  const Native = window.Worker;
  let next = 0;
  window.Worker = class extends Native {
    constructor(url, options) {
      super(url, options);
      this.__id = ++next;
      this.__url = String(url);
      window.__workerLog.push({ id: this.__id, event: 'construct', url: this.__url });
    }
    terminate() {
      window.__workerLog.push({ id: this.__id, event: 'terminate', url: this.__url });
      return super.terminate();
    }
  };
}
const pageWorkerLog = (page) => page.evaluate(() => window.__workerLog.filter((e) => e.url.includes('/analyze/worker.js')));

test('cancel during the index terminates the worker, returns to "choose genome", and a new run works', async (t) => {
  const context = await browser.newContext();
  await context.addInitScript(recordWorkers);
  const { page, errors } = await open(t, { context });
  hooks.kallisto = 'hold'; // the kallisto module does not arrive until it is released
  await chooseYeast(page);
  await waitPhase(page, 'index', 30_000);
  const logBefore = await pageWorkerLog(page);
  assert.deepEqual(logBefore.map((e) => e.event), ['construct'], JSON.stringify(logBefore));
  const oldId = logBefore[0].id;
  // Diagnostic only (not asserted): the Playwright close event of the old page worker. CI run 45042394266: it did not come in Firefox within 5 s.
  const pwWorkers = page.workers().filter((w) => w.url().includes('/analyze/worker.js'));
  const closed = pwWorkers.length === 1
    ? Promise.race([new Promise((r) => pwWorkers[0].once('close', () => r('yes'))), new Promise((r) => setTimeout(() => r('no'), 5000))])
    : Promise.resolve(`not checked (${pwWorkers.length} page workers listed)`);

  await page.click('#genome-cancel');
  assert.equal(await phase(page), 'choose');
  assert.equal((await appState(page)).indexReady, false);
  // The contract: the app called terminate() once on the old page worker and made exactly one new page worker.
  const logAfter = await pageWorkerLog(page);
  const terminated = logAfter.filter((e) => e.event === 'terminate');
  const made = logAfter.filter((e) => e.event === 'construct');
  assert.deepEqual(terminated.map((e) => e.id), [oldId], JSON.stringify(logAfter));
  assert.equal(made.length, 2, JSON.stringify(logAfter));
  const newId = made[1].id;
  assert.notEqual(newId, oldId);
  assert.deepEqual(await page.evaluate(() => window.__replaced), ['cancel']);

  t.diagnostic(`${ENGINE}: Playwright close event of the old page worker within 5 s: ${await closed}`);
  t.diagnostic(`${ENGINE}: held kallisto.mjs request closed by the browser before release: ${hooks.heldAborted ? 'yes' : 'no'}`);
  hooks.kallisto = 'serve';
  hooks.release?.(); // the old run gets its kallisto now; it must not come back to life
  hooks.release = null;
  await page.waitForTimeout(1500);
  assert.equal(await phase(page), 'choose');
  assert.equal((await appState(page)).indexReady, false);
  t.diagnostic(`${ENGINE}: page workers listed by Playwright 1.5 s after cancel: ${page.workers().filter((w) => w.url().includes('/analyze/worker.js')).length}`);

  await chooseYeast(page);
  await waitPhase(page, 'ready');
  assert.equal(await page.textContent('#ready-transcripts'), String(EXPECTED.transcripts));
  // The new run used the new worker: no other page worker was made or terminated.
  assert.deepEqual((await pageWorkerLog(page)).map((e) => `${e.event}:${e.id}`), [`construct:${oldId}`, `terminate:${oldId}`, `construct:${newId}`]);
  assert.deepEqual(errors, []);
});

test('offline, 404 and 5xx downloads give a message, not a blank page', async (t) => {
  for (const [how, pattern] of [
    ['offline', /the genome FASTA: the download from ftp\.ncbi\.nlm\.nih\.gov failed \(.+\)\. Check the network connection and try again\./],
    ['404', /the genome FASTA: the server ftp\.ncbi\.nlm\.nih\.gov answered HTTP 404 \(file not found\)\./],
    ['503', /the genome FASTA: the server ftp\.ncbi\.nlm\.nih\.gov answered HTTP 503 \(a server error; try again later\)\./],
  ]) {
    const { page } = await open(t, {
      ncbi: (route) => (how === 'offline' ? route.abort('internetdisconnected')
        : route.fulfill({ status: Number(how), body: 'x', headers: { 'access-control-allow-origin': '*' } })),
    });
    await chooseYeast(page);
    assert.match(await errorText(page), pattern, how);
  }
});

test('a kallisto out-of-memory failure gives a message that names the command-line route', async (t) => {
  const { page, errors } = await open(t);
  hooks.kallisto = 'oom';
  await chooseYeast(page);
  const text = await errorText(page);
  assert.match(text, /kallisto index stopped \(RangeError: WebAssembly\.Memory\(\): could not allocate memory\)\./);
  assert.match(text, /Use the command-line tool instead .*qrnamap fetch --source ncbi --id GCF_000146045\.2/);
  hooks.kallisto = 'serve';
  await chooseYeast(page); // the page worker still works
  await waitPhase(page, 'ready');
  assert.deepEqual(errors, []);
});

test('a catalog download that is not gzip (an HTML error page with HTTP 200) gives a clear message', async (t) => {
  const { page } = await open(t, {
    ncbi: (route) => route.fulfill({ status: 200, body: '<!DOCTYPE html><html><body>Service unavailable</body></html>', headers: { 'access-control-allow-origin': '*', 'content-type': 'text/html' } }),
  });
  await chooseYeast(page);
  assert.match(await errorText(page), /the genome FASTA: the file from ftp\.ncbi\.nlm\.nih\.gov is not a gzip file from NCBI \(\d+ bytes; it starts with "<!DOCTYPE html>/);
});

test('qrnamap:worker-replaced fires on cancel and on a genome change, and nextId() never repeats an id of the genome step', async (t) => {
  const context = await browser.newContext();
  // Record the id of every message the page sends to a worker.
  await context.addInitScript(() => {
    window.sentIds = [];
    const post = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function (message, ...rest) { if (message && typeof message.id === 'number') window.sentIds.push(message.id); return post.call(this, message, ...rest); };
    window.replaced = [];
    window.addEventListener('qrnamap:worker-replaced', (e) => window.replaced.push(e.detail.reason));
  });
  const { page, errors } = await open(t, { context });
  const before = await page.evaluate(() => [window.qrnamapApp.nextId(), window.qrnamapApp.nextId()]);
  hooks.kallisto = 'hold';
  await chooseYeast(page);
  await waitPhase(page, 'index', 30_000);
  await page.click('#genome-cancel');
  assert.deepEqual(await page.evaluate(() => window.replaced), ['cancel']);
  hooks.kallisto = 'serve';
  hooks.release?.();
  hooks.release = null;
  await chooseYeast(page);
  await waitPhase(page, 'ready');
  await page.click('#genome-change');
  assert.deepEqual(await page.evaluate(() => window.replaced), ['cancel', 'change']);
  assert.equal((await appState(page)).indexReady, false);
  const after = await page.evaluate(() => window.qrnamapApp.nextId());
  const genomeIds = await page.evaluate(() => window.sentIds);
  assert.equal(genomeIds.length, 2);
  const all = [...before, ...genomeIds, after];
  assert.equal(new Set(all).size, all.length, `ids: ${all.join(',')}`);
  assert.deepEqual(errors, []);
});

test('a damaged stored last genome is ignored', async (t) => {
  for (const value of [{ source: 'catalog', accession: 42, label: { x: 1 } }, { source: 'upload', fasta: null, gff: ['a'] }, 'not an object']) {
    const context = await browser.newContext();
    const first = await open(t, { context });
    await first.page.evaluate((v) => new Promise((done, fail) => {
      const request = indexedDB.open('qrnamap', 1);
      request.onupgradeneeded = () => request.result.createObjectStore('kv');
      request.onerror = () => fail(request.error);
      request.onsuccess = () => {
        const tx = request.result.transaction('kv', 'readwrite');
        tx.objectStore('kv').put(v, 'last-genome');
        tx.oncomplete = () => { request.result.close(); done(); };
        tx.onerror = () => fail(tx.error);
      };
    }), value);
    const again = await open(t, { context });
    assert.equal(await again.page.textContent('#genome-note'), '', JSON.stringify(value));
    assert.equal(await again.page.inputValue('#genome-menu'), '');
    assert.deepEqual(again.errors, []);
  }
});

test('a genome_limit_bases value below 1 is ignored', async (t) => {
  const { page } = await open(t, { query: '?genome_limit_bases=0.9' });
  await chooseYeast(page);
  await waitPhase(page, 'ready');
});

test('real NCBI download of S. cerevisiae S288C (REAL_NCBI=1 only)', { skip: process.env.REAL_NCBI !== '1' }, async (t) => {
  const context = await browser.newContext();
  const page = await context.newPage();
  t.after(() => context.close());
  await page.goto(`${base}/analyze/`);
  await page.waitForSelector('#step-genome[data-ready="menu"]', { timeout: 30_000 });
  const started = Date.now();
  await chooseYeast(page);
  await page.waitForSelector('#step-genome[data-phase="ready"], #genome-error:not([hidden])', { timeout: 600_000 });
  const seconds = (Date.now() - started) / 1000;
  const s = await appState(page);
  console.log(`real NCBI S288C in ${ENGINE}: ${seconds.toFixed(1)} s; phase ${s.phase}; genes ${s.nGenes}; transcripts ${s.nTranscripts}; ` +
    `error: ${await page.isHidden('#genome-error') ? 'none' : await page.textContent('#genome-error')}; index time: ${await page.textContent('#ready-time')}`);
  assert.equal(s.phase, 'ready');
  assert.equal(s.nTranscripts, 6478); // docs/browser/page-v1-measurements.md, "Genomes"
  assert.equal(s.nGenes, 6477);
});
