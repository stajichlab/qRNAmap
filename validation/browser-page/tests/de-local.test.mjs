// DESeq2 on this computer from the page's results: (1) the page's downloads (counts, samples.tsv, de-command.txt) run
// through the real `qrnamap de`; (2) the desktop app's "Run DESeq2 (R) on this computer" button (stand-in Tauri bridge).
//   DIST=<dist-site> ENGINE=chromium node --test validation/browser-page/tests/de-local.test.mjs
// QRNAMAP_BIN: the qrnamap binary (default target/release/qrnamap; test 1 is skipped without one).
// QRNAMAP_RSCRIPT: an Rscript with DESeq2 (for example .pixi/envs/de/bin/Rscript). With it, test 1 runs DESeq2 and checks the
// results; without it, test 1 checks that `qrnamap de` accepts every input and stops only because R is missing.
import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import http from 'node:http';
import { spawnSync } from 'node:child_process';
import { readFileSync, writeFileSync, statSync, existsSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve, extname, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const here = dirname(fileURLToPath(import.meta.url));
const FIX = resolve(here, '../fixtures');
const ROOT = resolve(here, '../../..');
const root = resolve(process.env.DIST || join(ROOT, 'dist-site'));
const ENGINE = process.env.ENGINE || 'chromium';
const WAIT_MS = 120_000;
const QRNAMAP = process.env.QRNAMAP_BIN || join(ROOT, 'target/release/qrnamap');
const RSCRIPT = process.env.QRNAMAP_RSCRIPT || '';

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
  assert.ok(existsSync(join(root, 'analyze/app.js')), `no dist-site at ${root}: run scripts/build-site.sh`);
  await new Promise((r) => server.listen(0, '127.0.0.1', r));
  base = `http://127.0.0.1:${server.address().port}`;
  browser = await { chromium, firefox, webkit }[ENGINE].launch();
});
after(async () => {
  await browser?.close();
  server.close();
});

/** The stand-in desktop bridge: `mode.r` false = no R with DESeq2; `mode.csv` = what run_deseq2_local returns. */
function fakeTauri(mode) {
  const calls = [];
  window.__deCalls = calls;
  window.__TAURI_INTERNALS__ = {
    invoke: async (cmd, args) => {
      calls.push({ cmd, args });
      switch (cmd) {
        case 'check_native_capabilities': return { available: false, index_compatible: false, warning: 'test' };
        case 'check_deseq2':
          return mode.r ? { available: true, rscript: '/opt/R/bin/Rscript', version: '1.50.2', detail: '' }
            : { available: false, rscript: null, version: null, detail: 'no Rscript was found' };
        case 'run_deseq2_local': return { results_csv: mode.csv, deseq2_version: '1.50.2', rscript: '/opt/R/bin/Rscript', out_dir: '/tmp/qrnamap_desktop/deseq2-1' };
        default: return null;
      }
    },
  };
}

/** Load the fixture genome and count four samples, two per condition (DESeq2 needs replicates). */
async function runFourSamples(t, mode) {
  const context = await browser.newContext({ acceptDownloads: true });
  t.after(() => context.close());
  if (mode) await context.addInitScript(fakeTauri, mode);
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  await page.goto(`${base}/analyze/?prebuilt=0`);
  await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 30_000 });
  await page.check('input[name="genome-source"][value="upload"]');
  await page.setInputFiles('#upload-fasta', join(FIX, 'mito.fna.gz'));
  await page.setInputFiles('#upload-gff', join(FIX, 'mito.gff3'));
  await page.click('#upload-load');
  await page.waitForSelector('#step-genome[data-phase="ready"]', { timeout: WAIT_MS });
  const file = (name, from) => ({ name, mimeType: 'application/octet-stream', buffer: readFileSync(join(FIX, from)) });
  const rows = [
    { name: 'C1', condition: 'ctl', m1: file('c1_R1.fastq.gz', 's1_R1.fastq.gz'), m2: file('c1_R2.fastq.gz', 's1_R2.fastq.gz') },
    { name: 'C2', condition: 'ctl', m1: file('c2_R1.fastq', 's2_R1.fastq'), m2: file('c2_R2.fastq', 's2_R2.fastq') },
    { name: 'T1', condition: 'trt', m1: file('t1.fastq.gz', 's3.fastq.gz') },
    { name: 'T2', condition: 'trt', m1: file('t2_R1.fastq', 's2_R2.fastq'), m2: file('t2_R2.fastq', 's2_R1.fastq') },
  ];
  for (const [i, r] of rows.entries()) {
    while (await page.locator('#sample-table tbody tr').count() <= i) await page.click('#add-sample');
    const tr = page.locator('#sample-table tbody tr').nth(i);
    await tr.locator('.sample-name').fill(r.name);
    await tr.locator('.cov').first().fill(r.condition);
    await tr.locator('td[data-mate="1"] input').setInputFiles(r.m1);
    if (r.m2) await tr.locator('td[data-mate="2"] input').setInputFiles(r.m2);
  }
  await page.fill('#frag-mean', '200');
  await page.fill('#frag-sd', '25');
  await page.click('#run-samples');
  await page.waitForFunction(() => window.qrnamapApp?.state?.samplesPhase === 'done', null, { timeout: WAIT_MS });
  await page.waitForSelector('#volcano-container svg', { timeout: 30_000 });
  return { page, errors };
}

/** Split a shell command made by the page's de_command (words, or single-quoted words). */
function shellWords(cmd) {
  const words = [];
  const re = /'((?:[^']|'\\'')*)'|(\S+)/g;
  let m;
  while ((m = re.exec(cmd))) words.push(m[1] !== undefined ? m[1].replace(/'\\''/g, "'") : m[2]);
  return words;
}

test('the page downloads run through `qrnamap de` (DESeq2 with QRNAMAP_RSCRIPT, else every input accepted)', { skip: !existsSync(QRNAMAP) && `no qrnamap binary at ${QRNAMAP}` }, async (t) => {
  const { page, errors } = await runFourSamples(t);
  const { counts, samples, command } = await page.evaluate(() => ({
    counts: window.qrnamapResults.getCountsTsv(),
    samples: window.qrnamapResults.getSamplesTsv(),
    command: window.qrnamapResults.getDeCommand().trim(),
  }));
  const words = shellWords(command);
  assert.deepEqual(words.slice(0, 2), ['qrnamap', 'de']);
  const flag = (f) => words[words.indexOf(f) + 1];
  assert.equal(flag('--design'), '~ condition');
  assert.equal(flag('--contrast'), 'condition,trt,ctl');

  const dir = mkdtempSync(join(tmpdir(), 'qrnamap-de-handoff-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  writeFileSync(join(dir, flag('--counts')), counts);
  writeFileSync(join(dir, flag('--samples')), samples);
  // TMPDIR: qrnamap de needs a temporary folder ($SCRATCH or $TMPDIR), and CI runners set neither.
  const env = { ...process.env, TMPDIR: process.env.TMPDIR || dir, QRNAMAP_RSCRIPT: RSCRIPT || join(dir, 'no-such-dir', 'Rscript') };
  const run = spawnSync(QRNAMAP, words.slice(1), { cwd: dir, env, encoding: 'utf8' });

  if (RSCRIPT) {
    assert.equal(run.status, 0, `qrnamap de failed:\n${run.stderr}`);
    assert.match(run.stdout, /^DESeq2 \d+\.\d+\.\d+/);
    const results = readFileSync(join(dir, 'de.results.csv'), 'utf8').trim().split('\n');
    const nGenes = counts.trim().split('\n').length - 1;
    assert.equal(results.length - 1, nGenes, 'one DESeq2 row per gene of the page matrix');
    assert.match(results[0], /"log2FoldChange".*"padj"/);
    for (const f of ['de.normalized.csv', 'de.vst.csv', 'de.versions.txt']) assert.ok(existsSync(join(dir, f)), f);
  } else {
    // The sheet, the design, the contrast and the count matrix passed every check: only starting R failed.
    assert.notEqual(run.status, 0);
    assert.match(run.stderr, /cannot start 'Rscript'/, run.stderr);
  }
  assert.deepEqual(errors, []);
});

test('desktop app: "Run DESeq2 (R)" shows DESeq2 results next to the quick screen, with NA kept', async (t) => {
  const csv = '"gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"\n' +
    '"gene-COX1",500,3.2,0.4,8,1e-15,1e-13\n"gene-COX2",40,-2.5,0.5,-5,5e-7,2e-5\n"gene-ATP6",3,0.2,1.2,0.16,0.87,NA\n';
  const { page, errors } = await runFourSamples(t, { r: true, csv });
  await page.waitForSelector('#de-r-box:not([hidden])');
  await page.waitForFunction(() => /R with DESeq2 1\.50\.2 found/.test(document.getElementById('de-r-status').textContent));
  assert.match(await page.textContent('#de-comparison-label'), /Method: quick screen/);
  assert.equal(await page.isHidden('#de-method-choice'), true);

  await page.click('#de-run-r');
  await page.waitForFunction(() => /genes at padj < 0\.05/.test(document.getElementById('de-r-status').textContent));
  const call = await page.evaluate(() => window.__deCalls.find((c) => c.cmd === 'run_deseq2_local').args.req);
  assert.equal(call.design, '~ condition');
  assert.equal(call.contrast, 'condition,trt,ctl');
  assert.match(call.counts_tsv, /^gene\tC1\tC2\tT1\tT2\n/);
  assert.match(call.samples_tsv, /^sample\tfastq1\tfastq2\tcondition\n/);
  assert.match(await page.textContent('#de-r-status'), /2 genes at padj < 0\.05 of 3\. All outputs .* \/tmp\/qrnamap_desktop\/deseq2-1/);
  assert.match(await page.textContent('#de-comparison-label'), /Method: DESeq2 1\.50\.2 \(R on this computer\)/);
  assert.match(await page.textContent('#volcano-summary'), /1 upregulated, 1 downregulated, 1 unchanged out of 3 genes/);
  const tsv = await page.evaluate(() => window.qrnamapResults.getDeTsv());
  assert.match(tsv, /gene-ATP6\t[^\n]*\tNA\tno\n/, 'an NA padj is written as NA and is not significant');

  // Back to the quick screen, and the label says so.
  await page.check('input[name="de-method"][value="screen"]');
  assert.match(await page.textContent('#de-comparison-label'), /Method: quick screen/);
  assert.doesNotMatch(await page.textContent('#volcano-summary'), /out of 3 genes/);
  assert.deepEqual(errors, []);
});

test('desktop app without R and DESeq2: the button is off and the message says how to get them', async (t) => {
  const { page, errors } = await runFourSamples(t, { r: false, csv: '' });
  await page.waitForFunction(() => /was not found/.test(document.getElementById('de-r-status').textContent));
  assert.equal(await page.isDisabled('#de-run-r'), true);
  assert.match(await page.textContent('#de-r-status'), /no Rscript was found.*QRNAMAP_RSCRIPT/);
  assert.deepEqual(errors, []);
});

test('in a browser (not the desktop app) there is no DESeq2 button', async (t) => {
  const { page } = await runFourSamples(t);
  assert.equal(await page.isHidden('#de-r-box'), true);
});
