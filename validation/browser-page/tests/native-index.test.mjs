// Browser tests of the desktop app's native index (genome step, site/analyze/genome.js prepareNativeIndex).
//   DIST=<absolute path of dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/native-index.test.mjs
// The Tauri bridge is a stand-in (window.__TAURI_INTERNALS__) that records each call and answers as a native kallisto that
// can, or cannot, read the page's kallisto.wasm index. kallisto index files are not portable between builds of 0.52.0
// (docs/browser/page-v1-measurements.md), so when the native kallisto cannot read it the page must hand over the transcript
// FASTA for a native build. The Rust side (check, build, cache) is tested in apps/desktop (cargo test).
// NCBI and S3 are never contacted: page.route answers them with the yeast mitochondrial fixtures.
import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import http from 'node:http';
import { gzipSync } from 'node:zlib';
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
const READY_MS = 120_000;
const fixture = (name) => readFileSync(join(FIX, name));

const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
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
  assert.ok(existsSync(join(root, 'analyze/index.html')), `no dist-site at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;
  browser = await { chromium, firefox, webkit }[ENGINE].launch();
});
after(async () => {
  await browser?.close();
  server.close();
});

/** The stand-in Tauri bridge. `mode.pageIndexLoads`: the native kallisto reads the page's index; `mode.cached`: a native index is kept. */
function fakeTauri(mode) {
  const calls = [];
  window.__nativeCalls = calls;
  window.__TAURI_INTERNALS__ = {
    invoke: async (cmd, args, options) => {
      const raw = args instanceof Uint8Array;
      calls.push({ cmd, headers: options?.headers || {}, first: raw ? String.fromCharCode(args[0]) : null, bytes: raw ? args.length : 0, args: raw ? null : args });
      switch (cmd) {
        case 'check_native_capabilities':
          return { available: true, index_compatible: true, kallisto_path: '/usr/bin/kallisto', kallisto_version: 'kallisto, version 0.52.0', num_cpus: 4, temp_dir: '/tmp', warning: null };
        case 'write_native_index': return { index_path: '/tmp/qrnamap_desktop/transcripts.idx', index_size: args.length };
        case 'check_native_index':
          return mode.pageIndexLoads ? { ok: true, detail: '' } : { ok: false, detail: '/usr/bin/kallisto cannot read this index (stopped by a signal: std::bad_array_new_length)' };
        case 'cached_native_index': return mode.cached ? { index_path: `/tmp/qrnamap_desktop/native/${args.key}.kallisto-0.52.0.idx`, index_size: 9 } : null;
        case 'build_native_index': return { index_path: `/tmp/qrnamap_desktop/native/${options.headers['x-qrnamap-key']}.kallisto-0.52.0.idx`, index_size: 9 };
        default: return null;
      }
    },
  };
}

/** Pre-built index files for the yeast accession (the index bytes are not read: the native side is the stand-in). */
const PREBUILT = {
  'idx.gz': gzipSync(Buffer.from('not a real index')),
  'tx2gene.tsv.gz': gzipSync(Buffer.from('rna-Q0045\tgene-Q0045\nrna-Q0050\tgene-Q0050\n')),
  'gene_info.tsv.gz': gzipSync(Buffer.from('gene_id\tname\tdescription\ngene-Q0045\tCOX1\tcytochrome c oxidase subunit 1\n')),
  'genes.tsv.gz': gzipSync(Buffer.from('gene_id\tname\ngene-Q0045\tCOX1\ngene-Q0050\tAI1\n')),
};

async function open(t, mode) {
  const context = await browser.newContext();
  t.after(() => context.close());
  await context.addInitScript(fakeTauri, mode);
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  const ncbi = [];
  const s3 = [];
  const cors = { 'access-control-allow-origin': '*' };
  await page.route('https://ftp.ncbi.nlm.nih.gov/**', (route) => {
    const url = route.request().url();
    ncbi.push(url);
    if (url === `${NCBI_DIR}GCF_000146045.2_R64_genomic.fna.gz`) return route.fulfill({ status: 200, body: fixture('mito.fna.gz'), headers: cors });
    if (url === `${NCBI_DIR}GCF_000146045.2_R64_genomic.gff.gz`) return route.fulfill({ status: 200, body: fixture('mito.gff3.gz'), headers: cors });
    return route.fulfill({ status: 404, body: 'not found', headers: cors });
  });
  await page.route('https://s3-west.nrp-nautilus.io/**', (route) => {
    const url = route.request().url();
    s3.push(url);
    const suffix = Object.keys(PREBUILT).find((s) => url.endsWith(`/${YEAST}.${s}`));
    return suffix ? route.fulfill({ status: 200, body: PREBUILT[suffix], headers: cors }) : route.fulfill({ status: 404, body: 'not found', headers: cors });
  });
  await page.goto(`${base}/analyze/`);
  await page.waitForSelector('#step-genome[data-ready="menu"]', { timeout: 30_000 });
  return { page, errors, ncbi, s3 };
}

async function upload(page) {
  await page.check('input[name="genome-source"][value="upload"]');
  await page.setInputFiles('#upload-fasta', join(FIX, 'mito.fna.gz'));
  await page.setInputFiles('#upload-gff', join(FIX, 'mito.gff3'));
  await page.click('#upload-load');
  await page.waitForSelector('#step-genome[data-phase="ready"]', { timeout: READY_MS });
}

async function chooseYeast(page) {
  await page.fill('#genome-filter', YEAST);
  await page.selectOption('#genome-menu', YEAST);
  await page.click('#genome-load');
  await page.waitForSelector('#step-genome[data-phase="ready"]', { timeout: READY_MS });
}

const calls = (page) => page.evaluate(() => window.__nativeCalls.filter((c) => c.cmd !== 'check_native_capabilities'));
const nativeState = (page) => page.evaluate(() => {
  const { nativeIndexPath, nativeIndexNote, nativeIndexError } = window.qrnamapApp.state;
  return { nativeIndexPath, nativeIndexNote, nativeIndexError };
});

test('a native kallisto that reads the page index uses it: nothing is built', async (t) => {
  const { page, errors } = await open(t, { pageIndexLoads: true });
  await upload(page);
  assert.deepEqual((await calls(page)).map((c) => c.cmd), ['write_native_index', 'check_native_index']);
  assert.deepEqual(await nativeState(page), { nativeIndexPath: '/tmp/qrnamap_desktop/transcripts.idx', nativeIndexNote: '', nativeIndexError: '' });
  assert.deepEqual(errors, []);
});

test('a native kallisto that cannot read the page index gets a native index built from the kept transcripts', async (t) => {
  const { page, errors, ncbi } = await open(t, { pageIndexLoads: false });
  await upload(page);
  const c = await calls(page);
  assert.deepEqual(c.map((x) => x.cmd), ['write_native_index', 'check_native_index', 'build_native_index']);
  const build = c.at(-1);
  assert.equal(build.first, '>', 'the body is the transcript FASTA');
  assert.ok(build.bytes > 1000);
  assert.equal(build.headers['x-qrnamap-key'], 'upload');
  const s = await nativeState(page);
  assert.equal(s.nativeIndexPath, '/tmp/qrnamap_desktop/native/upload.kallisto-0.52.0.idx');
  assert.match(await page.textContent('#genome-note'), /built its own index/);
  assert.deepEqual(ncbi, []);
  assert.deepEqual(errors, []);
});

test('a pre-built genome: the transcripts for the native index are built from the NCBI files, and kept under the accession', async (t) => {
  const { page, errors, ncbi, s3 } = await open(t, { pageIndexLoads: false });
  await chooseYeast(page);
  assert.match(await page.textContent('#ready-time'), /pre-built/);
  assert.ok(s3.some((u) => u.endsWith(`${YEAST}.idx.gz`)));
  const c = await calls(page);
  assert.deepEqual(c.map((x) => x.cmd), ['write_native_index', 'check_native_index', 'cached_native_index', 'build_native_index']);
  assert.equal(c[0].headers['x-qrnamap-gzip'], '1', 'the pre-built index goes to disk as its gzip');
  assert.equal(c[2].args.key, YEAST);
  assert.equal(c[3].first, '>');
  assert.equal(c[3].headers['x-qrnamap-key'], YEAST);
  assert.equal(ncbi.length, 2, 'the genome FASTA and the GFF3, once each');
  assert.equal((await nativeState(page)).nativeIndexPath, `/tmp/qrnamap_desktop/native/${YEAST}.kallisto-0.52.0.idx`);
  assert.deepEqual(errors, []);
});

test('a native index kept from an earlier run is used again: nothing downloaded or built', async (t) => {
  const { page, errors, ncbi } = await open(t, { pageIndexLoads: false, cached: true });
  await chooseYeast(page);
  assert.deepEqual((await calls(page)).map((c) => c.cmd), ['write_native_index', 'check_native_index', 'cached_native_index']);
  assert.deepEqual(ncbi, []);
  const s = await nativeState(page);
  assert.equal(s.nativeIndexPath, `/tmp/qrnamap_desktop/native/${YEAST}.kallisto-0.52.0.idx`);
  assert.match(s.nativeIndexNote, /built for this genome before/);
  assert.deepEqual(errors, []);
});
