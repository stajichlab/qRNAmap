// Browser tests of step 2 (samples and quantification) of /analyze/.
//   DIST=<absolute path of dist-site> ENGINE=chromium|firefox|webkit node --test validation/browser-page/tests/samples.test.mjs
// A local static server serves dist-site. Each test uploads the fixture genome (the yeast mitochondrial genome, ../fixtures)
// in step 1 and waits for "index ready". The reads are the simulated fixtures of ../fixtures/make_reads.py:
// s1 (300 pairs, gzip), s2 (300 pairs, plain), s3 (300 single-end reads, gzip).
// Expected gene counts and mapping numbers: ../fixtures/expected_quant.tsv and expected_runinfo.tsv, made with
// `qrnamap quant` (kallisto 0.52.0) on the same reads and the same transcripts. The page must give the same gene counts.
// Test switches of the server for kallisto.mjs (after the index is built): `damage` serves a module that damages the
// abundance.tsv of the next n runs (the Firefox fault of the measurements); `hold` holds the next request (cancel test).
import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import http from 'node:http';
import { readFileSync, statSync, existsSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { gzipSync } from 'node:zlib';
import { join, resolve, extname, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const here = dirname(fileURLToPath(import.meta.url));
const FIX = resolve(here, '../fixtures');
const root = resolve(process.env.DIST || resolve(here, '../../../dist-site'));
const ENGINE = process.env.ENGINE || 'chromium';
const WAIT_MS = Number(process.env.WAIT_MS || 120_000);

const fixture = (name) => readFileSync(join(FIX, name));
const fileArg = (name, buffer) => ({ name, mimeType: 'application/octet-stream', buffer });

/** expected[case] = { counts: {gene: n}, pairs, n_processed, n_pseudoaligned } from the CLI. */
const expected = {};
for (const line of readFileSync(join(FIX, 'expected_runinfo.tsv'), 'utf8').trim().split('\n').slice(1)) {
  const [c, pairs, processed, pseudo] = line.split('\t');
  expected[c] = { counts: {}, pairs: Number(pairs), n_processed: Number(processed), n_pseudoaligned: Number(pseudo) };
}
for (const line of readFileSync(join(FIX, 'expected_quant.tsv'), 'utf8').trim().split('\n').slice(1)) {
  const [c, gene, n] = line.split('\t');
  expected[c].counts[gene] = Number(n);
}

// ------------------------------------------------------------------ static server with the kallisto.mjs switches
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm', '.json': 'application/json' };
const hooks = { damage: 0, damageWhat: 'abundance', hold: false, holdSkip: 0, release: null, kallistoRequests: 0 };
// A kallisto module that runs the real one, then writes the largest est_counts token of abundance.tsv three times
// (as Firefox wrote `9` as `999` in the measurements), or (WHAT = runinfo) cuts run_info.json in half.
const damaging = (what) => `const WHAT = '${what}';
import factory from '/vendor/kallisto/kallisto.mjs?real=1';
export default async function (options) {
  const m = await factory(options);
  const callMain = m.callMain;
  m.callMain = (args) => {
    try { return callMain(args); } finally {
      const p = '/work/kq/abundance.tsv';
      const info = '/work/kq/run_info.json';
      if (WHAT === 'runinfo' && m.FS.analyzePath(info).exists) {
        const t = m.FS.readFile(info, { encoding: 'utf8' });
        m.FS.writeFile(info, t.slice(0, Math.floor(t.length / 2)));
      } else if (m.FS.analyzePath(p).exists) {
        const lines = m.FS.readFile(p, { encoding: 'utf8' }).split('\\n');
        let best = 1;
        for (let i = 1; i < lines.length; i++) if (lines[i] && Number(lines[i].split('\\t')[3]) > Number(lines[best].split('\\t')[3])) best = i;
        const f = lines[best].split('\\t');
        f[3] = f[3].repeat(3);
        lines[best] = f.join('\\t');
        m.FS.writeFile(p, lines.join('\\n'));
      }
    }
  };
  return m;
}
`;
const server = http.createServer(async (req, res) => {
  const url = new URL(req.url, 'http://x');
  const path = decodeURIComponent(url.pathname);
  const headers = (type) => ({ 'content-type': type, 'cache-control': 'no-store' });
  if (path === '/vendor/kallisto/kallisto.mjs' && !url.searchParams.has('real')) {
    hooks.kallistoRequests++;
    if (hooks.damage > 0) {
      hooks.damage--;
      res.writeHead(200, headers('text/javascript')).end(damaging(hooks.damageWhat));
      return;
    }
    if (hooks.hold && hooks.holdSkip > 0) {
      hooks.holdSkip--;
    } else if (hooks.hold) {
      hooks.hold = false;
      await new Promise((r) => { hooks.release = r; });
    }
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
  assert.ok(existsSync(join(root, 'analyze/samples.js')), `no dist-site with the samples step at ${root}: run scripts/build-site.sh`);
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
/** A new page with the fixture genome index ready (upload path) and step 2 open. */
async function open(t, query = '') {
  const context = await browser.newContext();
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  t.after(async () => {
    hooks.damage = 0;
    hooks.damageWhat = 'abundance';
    hooks.hold = false;
    hooks.holdSkip = 0;
    hooks.release?.();
    hooks.release = null;
    await context.close();
  });
  await page.goto(analyzeUrl(query));
  await page.waitForSelector('#step-samples[data-ready="table"]', { timeout: 30_000 });
  assert.equal(await page.getAttribute('#step-samples', 'data-state'), 'locked');
  await page.check('input[name="genome-source"][value="upload"]');
  await page.setInputFiles('#upload-fasta', join(FIX, 'mito.fna.gz'));
  await page.setInputFiles('#upload-gff', join(FIX, 'mito.gff3'));
  await page.click('#upload-load');
  await page.waitForSelector('#step-genome[data-phase="ready"]', { timeout: WAIT_MS });
  assert.equal(await page.getAttribute('#step-samples', 'data-state'), 'active');
  await page.evaluate(() => {
    window.__samplesDone = [];
    window.__replaced = [];
    document.addEventListener('qrnamap:samples-done', (e) => window.__samplesDone.push(e.detail.samples.map((s) => s.status)));
    window.addEventListener('qrnamap:worker-replaced', (e) => window.__replaced.push(e.detail.reason));
  });
  return { page, errors };
}

const row = (page, i) => page.locator('#sample-table tbody tr').nth(i);

/** Fill table row i (adds rows as needed). mate1/mate2: a fixture name or a {name, buffer} object. */
async function fillRow(page, i, { name, condition, mate1, mate2 }) {
  while (await page.locator('#sample-table tbody tr').count() <= i) await page.click('#add-sample');
  const r = row(page, i);
  await r.locator('.sample-name').fill(name);
  if (condition !== undefined) await r.locator('.cov').first().fill(condition);
  const arg = (m) => (typeof m === 'string' ? resolve(FIX, m) : m); // an absolute path stays as it is
  if (mate1) await r.locator('td[data-mate="1"] input').setInputFiles(arg(mate1));
  if (mate2) await row(page, i).locator('td[data-mate="2"] input').setInputFiles(arg(mate2));
}

/** Press "Count reads" and wait for the end of the run. Returns the samples of the page state. */
async function run(page) {
  await page.evaluate(() => { window.qrnamapApp.state.samplesPhase = 'test-wait'; }); // so a phase left by an earlier run does not end the wait
  await page.click('#run-samples');
  await page.waitForFunction(() => ['done', 'cancelled', 'lost'].includes(window.qrnamapApp.state.samplesPhase), null, { timeout: WAIT_MS });
  return page.evaluate(() => window.qrnamapApp.state.samples);
}

const noErrors = (page) => page.$eval('#samples-errors', (e) => e.hidden);
const tableSnapshot = (page) => page.$eval('#sample-table', (t) => [...t.querySelectorAll('input[type="text"]')].map((i) => i.value).join('|') + `#${t.querySelectorAll('tbody tr').length}`);

function assertCounts(sample, c) {
  assert.equal(sample.status, 'done', sample.error);
  assert.deepEqual(sample.geneCounts, expected[c].counts, `gene counts of ${sample.name} against the CLI (${c})`);
  assert.equal(sample.nPairs, expected[c].pairs);
  assert.equal(sample.processed, expected[c].n_processed);
  assert.equal(sample.pseudoaligned, expected[c].n_pseudoaligned);
  assert.equal(sample.fraction, expected[c].n_pseudoaligned / expected[c].pairs);
}

// ------------------------------------------------------------------ tests
test('manual entry: three samples (gzip, plain, single-end) give the gene counts of the CLI', async (t) => {
  const { page, errors } = await open(t);
  assert.equal(await page.$('#preview-label'), null, 'no preview label before results exist');
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'S2', condition: 'trt', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  await fillRow(page, 2, { name: 'S3', condition: 'trt', mate1: 's3.fastq.gz' });
  // A single-end sample needs the fragment length: first without it.
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /Single-end samples \(no mate 2 file\) need the fragment length mean and SD/);
  await page.fill('#frag-mean', '200');
  await page.fill('#frag-sd', '20');
  const samples = await run(page);
  assert.equal(await noErrors(page), true);
  assertCounts(samples[0], 's1');
  assertCounts(samples[1], 's2');
  assertCounts(samples[2], 's3');
  assert.deepEqual(samples.map((s) => s.covariates), [{ condition: 'ctl' }, { condition: 'trt' }, { condition: 'trt' }]);
  assert.deepEqual(samples.map((s) => s.allPairs), [true, true, true]);
  // kallisto argv from kallisto_args: threads 1; single-end with -l and -s.
  assert.deepEqual(samples[0].args, ['quant', '-i', 'transcripts.idx', '-o', 'kq', '--plaintext', '-t', '1', 'r1.fastq', 'r2.fastq']);
  assert.deepEqual(samples[2].args, ['quant', '-i', 'transcripts.idx', '-o', 'kq', '--plaintext', '-t', '1', '--single', '-l', '200', '-s', '20', 'r1.fastq']);
  // Mapping numbers on the row, "all pairs (300)" because the files have fewer than N pairs.
  assert.match(await row(page, 0).locator('.status').textContent(), new RegExp(`^300 pairs processed, ${expected.s1.n_pseudoaligned} pseudoaligned \\(${(100 * expected.s1.n_pseudoaligned / 300).toFixed(1)}%\\)\\. All pairs \\(300\\)\\. Reads of 50 bases: N = 200,000\\.$`));
  assert.match(await row(page, 2).locator('.status').textContent(), /^300 reads processed, .* All pairs \(300\)\. Reads of 50 bases: N = 200,000\.$/);
  assert.match(await page.textContent('#preview-label'), /preview of the first 200,000 pairs of each sample/);
  assert.match(await page.textContent('#construction-notice'), /Genome step and quantification tested in headless Chromium, Firefox and WebKit on a CI runner \(small fixture genome and reads\); memory on full-size reads is not measured yet; results view is not built yet\./);
  assert.deepEqual(samples.map((s) => s.warning), ['', '', '']);
  assert.deepEqual(await page.evaluate(() => window.qrnamapApp.state.preview), { n: 200000, label: 'preview of the first 200,000 pairs' });
  assert.deepEqual(await page.evaluate(() => window.__samplesDone), [['done', 'done', 'done']]);
  assert.equal(await page.getAttribute('#step-samples', 'data-state'), 'done');
  assert.equal(await page.getAttribute('#step-results', 'data-state'), 'active');
  assert.equal(await page.textContent('#samples-summary'), '3 of 3 samples counted.');
  assert.deepEqual(errors, []);
});

test('one sample, cut at N = 100 pairs: exactly 100 pairs are read, counts equal the CLI on the first 100 pairs', async (t) => {
  const { page, errors } = await open(t);
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await page.fill('#preview-pairs', '100');
  const [s] = await run(page);
  assertCounts(s, 's1_first100');
  assert.equal(s.allPairs, false);
  assert.match(await row(page, 0).locator('.status').textContent(), /100 pairs processed, .* Preview of the first 100 pairs\. Reads of 50 bases: N = 100\.$/);
  assert.match(await page.textContent('#preview-label'), /preview of the first 100 pairs/);
  // N cannot be raised above the hard limit.
  await page.fill('#preview-pairs', '200001');
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /Read pairs per sample must be a whole number from 1 to 200,000/);
  assert.deepEqual(errors, []);
});

const SHEET = '# a comment line\r\n\r\nsample\tfastq1\tfastq2\tcondition\tdose:numeric\r\n' +
  'S1\t/data/run 1/s1_R1.fastq.gz\t/data/run 1/s1_R2.fastq.gz\tctl\t0\r\n' +
  'S2\treads/s2_R1.fastq\treads/s2_R2.fastq\ttrt\t1.5\r\n';

test('sheet upload: CRLF, comments, blank lines, folders and spaces; files matched by name; the table is filled; counts equal the CLI', async (t) => {
  const { page, errors } = await open(t);
  await page.setInputFiles('#sheet-fastqs', ['s2_R2.fastq', 's1_R1.fastq.gz', 's1_R2.fastq.gz', 's2_R1.fastq'].map((n) => join(FIX, n)));
  await page.setInputFiles('#sheet-file', fileArg('samples.tsv', Buffer.from(SHEET)));
  await page.waitForSelector('#samples-panel[data-sheet="samples.tsv"]');
  assert.equal(await page.$eval('#sheet-error', (e) => e.hidden), true);
  assert.deepEqual(await page.$$eval('#sample-table .col-name', (es) => es.map((e) => e.value)), ['condition', 'dose:numeric']);
  assert.deepEqual(await page.$$eval('#sample-table .sample-name', (es) => es.map((e) => e.value)), ['S1', 'S2']);
  assert.deepEqual(await page.$$eval('#sample-table .cov', (es) => es.map((e) => e.value)), ['ctl', '0', 'trt', '1.5']);
  assert.deepEqual(await page.$$eval('#sample-table tbody tr', (es) => es.map((e) => e.dataset.missing)), ['false', 'false']);
  assert.match(await row(page, 0).locator('td[data-mate="1"] .mate-name').textContent(), /^s1_R1\.fastq\.gz \(matched to the sheet\)$/);
  const samples = await run(page);
  assertCounts(samples[0], 's1');
  assertCounts(samples[1], 's2');
  assert.deepEqual(samples.map((s) => s.covariates), [{ condition: 'ctl', dose: '0' }, { condition: 'trt', dose: '1.5' }]);
  assert.deepEqual(errors, []);
});

test('fill from FASTQ files: mates paired by name, one row per sample; skipped files and lanes reported; counts equal the CLI', async (t) => {
  const { page, errors } = await open(t);
  await fillRow(page, 0, { name: 'kept', condition: 'ctl' }); // a row with something in it stays
  await page.setInputFiles('#fill-fastqs', [
    ...['s2_R2.fastq', 's1_R1.fastq.gz', 's3.fastq.gz', 's1_R2.fastq.gz', 's2_R1.fastq'].map((n) => fileArg(n, fixture(n))), // paths and buffers cannot mix
    fileArg('notes.txt', Buffer.from('x')),
    fileArg('L_L001_R1.fq.gz', Buffer.from('x')), fileArg('L_L002_R1.fq.gz', Buffer.from('x')),
    fileArg('L_L001_R2.fq.gz', Buffer.from('x')), fileArg('L_L002_R2.fq.gz', Buffer.from('x')),
  ]);
  assert.deepEqual(await page.$$eval('#sample-table .sample-name', (es) => es.map((e) => e.value)), ['kept', 's1', 's2', 's3']);
  assert.match(await page.textContent('#fill-status'), /^3 sample\(s\) added \(2 paired-end, 1 single-end\)/);
  const notes = await page.textContent('#sheet-error');
  assert.match(notes, /notes\.txt was left out: not a FASTQ file name/);
  assert.match(notes, /sample "L" was not added: it has 2 lanes per mate/);
  assert.match(await row(page, 1).locator('td[data-mate="1"] .mate-name').textContent(), /^s1_R1\.fastq\.gz$/);
  assert.match(await row(page, 1).locator('td[data-mate="2"] .mate-name').textContent(), /^s1_R2\.fastq\.gz$/);
  assert.match(await row(page, 3).locator('td[data-mate="2"] .mate-name').textContent(), /none \(single-end\)/);

  // A second fill does not add a sample twice.
  await page.setInputFiles('#fill-fastqs', ['s1_R1.fastq.gz', 's1_R2.fastq.gz'].map((n) => join(FIX, n)));
  assert.match(await page.textContent('#sheet-error'), /sample "s1" was not added: the table already has a row with that name/);
  assert.equal(await page.locator('#sample-table tbody tr').count(), 4);

  await page.click('#sample-table tbody tr:first-child .remove-row');
  for (const [i, c] of [[0, 'ctl'], [1, 'trt'], [2, 'trt']]) await row(page, i).locator('.cov').first().fill(c);
  await page.fill('#frag-mean', '200');
  await page.fill('#frag-sd', '20');
  const samples = await run(page);
  assertCounts(samples[0], 's1');
  assertCounts(samples[1], 's2');
  assertCounts(samples[2], 's3');
  assert.deepEqual(samples.map((s) => s.name), ['s1', 's2', 's3']);
  assert.deepEqual(errors, []);
});

test('fill from FASTQ files: _1/_2 mates can be turned off', async (t) => {
  const { page } = await open(t);
  const names = ['SRR9_1.fastq.gz', 'SRR9_2.fastq.gz'];
  await page.setInputFiles('#fill-fastqs', names.map((n) => fileArg(n, fixture('s1_R1.fastq.gz'))));
  assert.deepEqual(await page.$$eval('#sample-table .sample-name', (es) => es.map((e) => e.value)), ['SRR9']);
  await page.click('#sample-table tbody tr:first-child .remove-row');
  await page.uncheck('#fill-numeric-mates');
  await page.setInputFiles('#fill-fastqs', names.map((n) => fileArg(n, fixture('s1_R1.fastq.gz'))));
  assert.deepEqual(await page.$$eval('#sample-table .sample-name', (es) => es.map((e) => e.value)), ['SRR9_1', 'SRR9_2']);
});

test('sheet upload: a row with no matching file is flagged and gets a picker; picking the file clears the flag', async (t) => {
  const { page } = await open(t);
  await page.setInputFiles('#sheet-file', fileArg('samples.tsv', Buffer.from(SHEET)));
  await page.setInputFiles('#sheet-fastqs', ['s1_R1.fastq.gz', 's1_R2.fastq.gz'].map((n) => join(FIX, n)));
  await page.waitForFunction(() => document.querySelectorAll('#sample-table tbody tr[data-missing="true"]').length === 1);
  assert.deepEqual(await page.$$eval('#sample-table tbody tr', (es) => es.map((e) => e.dataset.missing)), ['false', 'true']);
  assert.match(await row(page, 1).locator('td[data-mate="1"] .mate-name').textContent(), /missing: the sheet names "s2_R1\.fastq"\. Pick that file\./);
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /Row 2 \(S2\): the sheet names "s2_R1\.fastq" for mate 1, but no picked file has that name\. Pick the file\./);
  assert.equal(await row(page, 1).getAttribute('data-invalid'), 'true');
  await row(page, 1).locator('td[data-mate="1"] input').setInputFiles(join(FIX, 's2_R1.fastq'));
  await row(page, 1).locator('td[data-mate="2"] input').setInputFiles(join(FIX, 's2_R2.fastq'));
  assert.equal(await row(page, 1).getAttribute('data-missing'), 'false');
});

test('two picked files with the same name, and one file name used by two rows, give clear errors', async (t) => {
  const { page } = await open(t);
  await page.setInputFiles('#sheet-fastqs', [fileArg('s1_R1.fastq.gz', fixture('s1_R1.fastq.gz')), fileArg('s1_R1.fastq.gz', fixture('s2_R1.fastq'))]);
  assert.match(await page.textContent('#sheet-error'), /two picked files have the same name "s1_R1\.fastq\.gz"\. The page matches files by name only: pick only one of them\./);
  const before = await tableSnapshot(page);
  const twice = 'sample\tfastq1\tfastq2\nA\tx/r1.fq\tx/r2.fq\nB\ty/r1.fq\ty/r3.fq\n';
  await page.setInputFiles('#sheet-file', fileArg('twice.tsv', Buffer.from(twice)));
  await page.waitForSelector('#sheet-error:not([hidden])');
  assert.match(await page.textContent('#sheet-error'), /the file name "r1\.fq" is used twice \(sample "A" fastq1 and sample "B" fastq1\)/);
  assert.equal(await tableSnapshot(page), before, 'a refused sheet changes nothing');
});

test('a bad sheet shows the validate_sheet message and changes nothing', async (t) => {
  const { page } = await open(t);
  await fillRow(page, 0, { name: 'Kept', condition: 'x' });
  const before = await tableSnapshot(page);
  const cases = [
    ['sample\tfastq1\nA\ta.fq\nA\tb.fq\n', /duplicate sample name "A"/],
    ['sample\tfastq2\nA\ta.fq\n', /required column "fastq1" is missing/],
    ['name\tfastq1\nA\ta.fq\n', /required column "sample" is missing/],
    ['sample\tfastq1\tfastq2\nA\ta1.fq,a2.fq\tb1.fq\n', /sample "A": fastq1 and fastq2 list different numbers of files/],
    ['sample\tfastq1\nbad name!\ta.fq\n', /name or value "bad name!" contains forbidden characters or has a bad form/],
    ['sample\tfastq1\nindex\ta.fq\n', /sample 'index': the name is reserved: qrnamap writes 'index' in the output directory\. Rename the sample in the sheet/],
    ['sample\tfastq1\nmapping_summary.tsv.zst.partial\ta.fq\n', /the name is reserved/],
  ];
  for (const [i, [text, message]] of cases.entries()) {
    await page.setInputFiles('#sheet-file', fileArg(`bad${i}.tsv`, Buffer.from(text)));
    await page.waitForFunction((n) => document.getElementById('sheet-error').textContent.includes(`bad${n}.tsv was not loaded`), i);
    assert.match(await page.textContent('#sheet-error'), message);
    assert.equal(await tableSnapshot(page), before, `case ${i} changed the table`);
  }
});

test('several lanes in one cell: "one file per mate in this version"', async (t) => {
  const { page } = await open(t);
  const before = await tableSnapshot(page);
  await page.setInputFiles('#sheet-file', fileArg('lanes.tsv', Buffer.from('sample\tfastq1\tfastq2\nA\tL1_1.fq.gz,L2_1.fq.gz\tL1_2.fq.gz,L2_2.fq.gz\n')));
  await page.waitForSelector('#sheet-error:not([hidden])');
  assert.match(await page.textContent('#sheet-error'), /sample "A": the fastq1 cell lists 2 files \(lanes\)\. This version of the page takes one file per mate\./);
  assert.equal(await tableSnapshot(page), before);
});

test('manual entry: duplicate, unsafe and reserved names are refused next to the rows', async (t) => {
  const { page } = await open(t);
  await fillRow(page, 0, { name: 'A', condition: 'x', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'A', condition: 'x', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /duplicate sample name "A"/);
  assert.deepEqual(await page.$$eval('#sample-table tbody tr', (es) => es.map((e) => e.dataset.invalid)), ['true', 'true']);
  await row(page, 1).locator('.sample-name').fill('slurm');
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /sample 'slurm': the name is reserved/);
  assert.match(await row(page, 1).locator('.status').textContent(), /the name is reserved/);
  await row(page, 1).locator('.sample-name').fill('2bad');
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /name or value "2bad" contains forbidden characters/);
  // The same file picked in two rows: one name used twice.
  await row(page, 1).locator('.sample-name').fill('B');
  await row(page, 1).locator('td[data-mate="1"] input').setInputFiles(join(FIX, 's1_R1.fastq.gz'));
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /the file name "s1_R1\.fastq\.gz" is used twice/);
  assert.equal(await page.evaluate(() => window.qrnamapApp.state.samplesPhase), 'edit', 'no run started');
});

test('mates with different record counts, a broken FASTQ record and an empty file give clear errors; the other samples still run', async (t) => {
  const { page, errors } = await open(t);
  const short = fixture('s2_R2.fastq').toString().split('\n').slice(0, 250 * 4).join('\n') + '\n';
  const broken = fixture('s2_R1.fastq').toString().split('\n');
  broken[4 * 9 + 2] = 'not a plus line';
  await fillRow(page, 0, { name: 'Short', condition: 'x', mate1: 's2_R1.fastq', mate2: fileArg('short_R2.fastq', Buffer.from(short)) });
  await fillRow(page, 1, { name: 'Broken', condition: 'x', mate1: fileArg('broken_R1.fastq', Buffer.from(broken.join('\n'))), mate2: 's2_R2.fastq' });
  await fillRow(page, 2, { name: 'Empty', condition: 'x', mate1: fileArg('empty_R1.fastq', Buffer.alloc(0)), mate2: fileArg('empty_R2.fastq', Buffer.from('@a\nA\n+\nI\n')) });
  await fillRow(page, 3, { name: 'Good', condition: 'x', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  const samples = await run(page);
  assert.deepEqual(samples.map((s) => s.status), ['error', 'error', 'error', 'done']);
  assert.match(samples[0].error, /The mates of sample Short have different record counts: mate 1 \(s2_R1\.fastq\) has 300 records, mate 2 \(short_R2\.fastq\) has 250 records\./);
  assert.match(samples[1].error, /Mate 1 of sample Broken \(broken_R1\.fastq\) is not a FASTQ file: record 10 \(line 39\) has no "\+" line\./);
  assert.match(samples[2].error, /Mate 1 of sample Empty \(empty_R1\.fastq\) is empty \(0 bytes\)\./);
  assertCounts(samples[3], 's1');
  assert.match(await row(page, 0).locator('.status').textContent(), /^Error: The mates of sample Short/);
  assert.deepEqual(errors, []);
});

test('a damaged abundance.tsv (the Firefox fault) is found: one re-run gives the right counts; damaged twice gives an error', async (t) => {
  const { page, errors } = await open(t);
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  hooks.damage = 1;
  let [s] = await run(page);
  assertCounts(s, 's1');
  assert.equal(s.reruns, 1);
  assert.match(await row(page, 0).locator('.status').textContent(), /kallisto ran twice: the first output failed the check\./);
  hooks.damage = 2;
  [s] = await run(page);
  assert.equal(s.status, 'error');
  assert.equal(s.kind, 'output');
  assert.match(s.error, /kallisto wrote a damaged output twice for sample S1 \(first run: the est_counts add up to .*, but run_info\.json says \d+ reads pseudoaligned; second run: /);
  assert.equal(hooks.damage, 0);
  // A run_info.json that cannot be read also starts the second run.
  hooks.damage = 1;
  hooks.damageWhat = 'runinfo';
  [s] = await run(page);
  assertCounts(s, 's1');
  assert.equal(s.reruns, 1);
  assert.equal(hooks.damage, 0);
  assert.deepEqual(errors, []);
});

test('cancel during quant stops the sample; the index stays, and a new run needs no new index', async (t) => {
  const { page, errors } = await open(t);
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'S2', condition: 'trt', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  hooks.hold = true;
  await page.click('#run-samples');
  await page.waitForFunction(() => /Running kallisto quant/.test(document.getElementById('samples-status').textContent), null, { timeout: WAIT_MS });
  assert.equal(await page.$eval('#samples-edit', (f) => f.disabled), true, 'the table is locked during a run');
  await page.click('#samples-cancel');
  await page.waitForFunction(() => window.qrnamapApp.state.samplesPhase === 'cancelled', null, { timeout: WAIT_MS });
  hooks.release?.();
  hooks.release = null;
  let samples = await page.evaluate(() => window.qrnamapApp.state.samples);
  assert.deepEqual(samples.map((s) => s.status), ['cancelled', 'not run']);
  assert.equal(await page.evaluate(() => window.qrnamapApp.state.indexReady), true);
  assert.deepEqual(await page.evaluate(() => window.__replaced), [], 'the page worker was not replaced');
  assert.deepEqual(await page.evaluate(() => window.__samplesDone), []);
  assert.equal(await page.$('#preview-label'), null);
  assert.equal(await page.$eval('#samples-edit', (f) => f.disabled), false);
  assert.match(await page.textContent('#samples-summary'), /^Cancelled\./);
  samples = await run(page);
  assertCounts(samples[0], 's1');
  assertCounts(samples[1], 's2');
  assert.deepEqual(errors, []);
});

test('choosing another genome clears the counts and locks the samples step', async (t) => {
  const { page } = await open(t);
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await run(page);
  assert.ok(await page.$('#preview-label'));
  await page.click('#genome-change');
  await page.waitForFunction(() => window.qrnamapApp.state.samples.length === 0);
  assert.equal(await page.$('#preview-label'), null);
  assert.equal(await page.getAttribute('#step-samples', 'data-state'), 'locked');
  assert.match(await page.textContent('#samples-errors'), /the counts were cleared/);
  assert.equal(await row(page, 0).locator('.sample-name').inputValue(), 'S1', 'the table itself is kept');
});

test('covariate columns can be added, renamed and removed; an empty column is left out of the sheet', async (t) => {
  const { page } = await open(t);
  await page.click('#add-column');
  await page.locator('#sample-table .col-name').nth(1).fill('batch');
  await fillRow(page, 0, { name: 'S1', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await row(page, 0).locator('.cov').nth(1).fill('b1');
  await page.click('#add-column'); // a third column, left empty: not in the sheet
  let [s] = await run(page);
  assert.deepEqual(s.covariates, { batch: 'b1' }, 'the empty condition column and the empty third column are left out');
  await page.locator('#sample-table .remove-col').first().click();
  assert.deepEqual(await page.$$eval('#sample-table .col-name', (es) => es.map((e) => e.value)), ['batch', '']);
  await page.locator('#sample-table .col-name').nth(1).fill('if');
  await row(page, 0).locator('.cov').nth(1).fill('v');
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /column "if" is a reserved R name/);
  await page.click('.remove-row');
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /Add at least one sample\./);
});

test('mates whose first read names differ give a warning, not a refusal; a file name with a comma is refused', async (t) => {
  const { page, errors } = await open(t);
  const renamed = fixture('s2_R2.fastq').toString().replace(/^@s2_/gm, '@other_');
  await fillRow(page, 0, { name: 'S2', condition: 'x', mate1: 's2_R1.fastq', mate2: fileArg('other_R2.fastq', Buffer.from(renamed)) });
  await fillRow(page, 1, { name: 'C', condition: 'x', mate1: fileArg('a,b_R1.fastq', fixture('s2_R1.fastq')), mate2: 's2_R2.fastq' });
  await page.click('#run-samples');
  assert.match(await page.textContent('#samples-errors'), /Row 2 \(C\): the mate 1 file name "a,b_R1\.fastq" has a comma/);
  assert.equal(await row(page, 1).getAttribute('data-invalid'), 'true');
  await page.locator('#sample-table .remove-row').nth(1).click();
  const [s] = await run(page);
  assertCounts(s, 's2');
  assert.equal(s.warning, 'the first read names of the two mates differ ("s2_000001" and "other_000001"). Check that the two files belong together.');
  assert.match(await row(page, 0).locator('.status').textContent(), /Warning: the first read names of the two mates differ/);
  assert.deepEqual(errors, []);
});

test('cancel at sample 2 keeps the counts of sample 1', async (t) => {
  const { page, errors } = await open(t);
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'S2', condition: 'trt', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  hooks.hold = true;
  hooks.holdSkip = 1; // the kallisto run of sample 1 goes through; the one of sample 2 is held
  await page.click('#run-samples');
  await page.waitForFunction(() => window.qrnamapApp.state.samples[0]?.status === 'done' && /Running kallisto quant/.test(document.getElementById('samples-status').textContent), null, { timeout: WAIT_MS });
  await page.click('#samples-cancel');
  await page.waitForFunction(() => window.qrnamapApp.state.samplesPhase === 'cancelled', null, { timeout: WAIT_MS });
  const samples = await page.evaluate(() => window.qrnamapApp.state.samples);
  assert.deepEqual(samples.map((s) => s.status), ['done', 'cancelled']);
  assertCounts(samples[0], 's1');
  assert.ok(await page.$('#preview-label'), 'the counts of sample 1 are results: the preview label is shown');
  assert.deepEqual(await page.evaluate(() => window.__samplesDone), []);
  assert.deepEqual(errors, []);
});

/**
 * Wait for `fn` in the page. On a timeout, throw an error that names the wait and holds the page state: samplesPhase, the
 * samples-status text, the sample statuses and errors, the kallisto.mjs requests of the test server, page errors, and the history
 * of the samples-status text with times (recorded by a MutationObserver, see watchStatus).
 */
async function waitOrDump(page, fn, what, errors, arg) {
  try {
    await page.waitForFunction(fn, arg, { timeout: WAIT_MS, polling: 10 });
  } catch (error) {
    const state = await page.evaluate(() => ({
      samplesPhase: window.qrnamapApp.state.samplesPhase,
      status: document.getElementById('samples-status').textContent,
      samples: window.qrnamapApp.state.samples.map((s) => ({ status: s.status, error: s.error, nPairs: s.nPairs, previewN: s.previewN })),
      history: window.__statusHistory || [],
    })).catch((e) => ({ evaluateFailed: e.message }));
    throw new Error(`wait "${what}" failed (${error.message.split('\n')[0]}). kallisto.mjs requests: ${hooks.kallistoRequests}. ` +
      `page errors: ${JSON.stringify(errors)}. state: ${JSON.stringify(state)}`);
  }
}

/** Record every change of the samples-status text with the ms since the call (window.__statusHistory). */
const watchStatus = (page) => page.evaluate(() => {
  const t0 = performance.now();
  const box = document.getElementById('samples-status');
  window.__statusHistory = [];
  new MutationObserver(() => window.__statusHistory.push(`${Math.round(performance.now() - t0)} ms: ${box.textContent}`))
    .observe(box, { childList: true, characterData: true, subtree: true });
  document.getElementById('samples-cancel').addEventListener('click', () => window.__statusHistory.push(`${Math.round(performance.now() - t0)} ms: CANCEL CLICKED`));
});

test('cancel while the reads are read (before kallisto runs)', async (t) => {
  const { page, errors } = await open(t);
  // 200,000 records of 150 bases (gzip level 1, in a temporary file): reading them takes long enough to cancel.
  const dir = mkdtempSync(join(tmpdir(), 'qrnamap-prep-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  const [big1, big2] = writeBig(dir, 200000);
  await fillRow(page, 0, { name: 'Big', condition: 'x', mate1: big1, mate2: big2 });
  const before = hooks.kallistoRequests;
  await watchStatus(page);
  await page.click('#run-samples');
  // Wait until the bulk read has started (a record count is shown), not only the read of the first record.
  await waitOrDump(page, () => /Reading mate \d \(.*\): [\d,]+ records/.test(document.getElementById('samples-status').textContent), 'status shows a record count', errors);
  await page.click('#samples-cancel');
  await waitOrDump(page, () => window.qrnamapApp.state.samplesPhase === 'cancelled', 'samplesPhase is cancelled', errors);
  const history = await page.evaluate(() => window.__statusHistory);
  t.diagnostic(`status history: ${history.join(' | ')}`);
  assert.equal(hooks.kallistoRequests, before, `kallisto did not start: the cancel came during the read. History: ${history.join(' | ')}`);
  assert.deepEqual(await page.evaluate(() => window.qrnamapApp.state.samples.map((s) => s.status)), ['cancelled']);
  assert.equal(await page.evaluate(() => window.qrnamapApp.state.indexReady), true);
  assert.deepEqual(await page.evaluate(() => window.__replaced), []);
  assert.deepEqual(errors, []);
});

/**
 * Two gzip files (mate 1 and 2) of `n` pairs of 150-base reads, in `dir`: windows of the fixture genome (mito.fna) at a fixed stride;
 * mate 2 is the reverse complement of the window 200 bases further on. Many pairs fall in transcripts, so some pseudoalign.
 */
function writeBig(dir, n) {
  const genome = fixture('mito.fna').toString().split('\n').filter((l) => !l.startsWith('>')).join('').toUpperCase();
  const comp = { A: 'T', C: 'G', G: 'C', T: 'A' };
  const rc = (x) => [...x].reverse().map((c) => comp[c] || 'N').join('');
  const qual = 'I'.repeat(150);
  const m1 = [];
  const m2 = [];
  for (let i = 0; i < n; i++) {
    const p = (i * 997) % (genome.length - 400);
    m1.push(`@big_${i}/1\n${genome.slice(p, p + 150)}\n+\n${qual}\n`);
    m2.push(`@big_${i}/2\n${rc(genome.slice(p + 200, p + 350))}\n+\n${qual}\n`);
  }
  const files = [join(dir, 'big_R1.fastq.gz'), join(dir, 'big_R2.fastq.gz')];
  writeFileSync(files[0], gzipSync(m1.join(''), { level: 1 }));
  writeFileSync(files[1], gzipSync(m2.join(''), { level: 1 }));
  return files;
}

test('the read length sets N: 150-base reads use 75,000 pairs, 50-base reads 200,000; the label shows both', async (t) => {
  const { page, errors } = await open(t);
  const dir = mkdtempSync(join(tmpdir(), 'qrnamap-long-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  const [big1, big2] = writeBig(dir, 75010);
  await fillRow(page, 0, { name: 'S1', condition: 'x', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'Long', condition: 'y', mate1: big1, mate2: big2 });
  const samples = await run(page);
  assertCounts(samples[0], 's1');
  assert.equal(samples[0].previewN, 200000);
  assert.equal(samples[0].readLength, 50);
  assert.equal(samples[1].status, 'done', samples[1].error);
  assert.deepEqual([samples[1].readLength, samples[1].previewN, samples[1].nPairs, samples[1].processed, samples[1].allPairs], [150, 75000, 75000, 75000, false]);
  assert.match(await row(page, 1).locator('.status').textContent(), /^75,000 pairs processed, .* Preview of the first 75,000 pairs\. Reads of 150 bases: N = 75,000\.$/);
  assert.deepEqual(await page.evaluate(() => window.qrnamapApp.state.preview),
    { n: null, label: 'preview of the first 75,000 to 200,000 pairs (the N of each sample is in its row)' });
  assert.match(await page.textContent('#preview-label'), /preview of the first 75,000 to 200,000 pairs/);
  assert.deepEqual(errors, []);
});

test('measurement option ?preview_pairs replaces the read-length limit (never above 200,000)', async (t) => {
  const { page, errors } = await open(t, '?preview_pairs=100005');
  const dir = mkdtempSync(join(tmpdir(), 'qrnamap-over-'));
  t.after(() => rmSync(dir, { recursive: true, force: true }));
  const [big1, big2] = writeBig(dir, 100010);
  await fillRow(page, 0, { name: 'Long', condition: 'y', mate1: big1, mate2: big2 });
  const [s] = await run(page);
  assert.equal(s.status, 'done', s.error);
  assert.deepEqual([s.previewN, s.nPairs], [100005, 100005]);
  assert.deepEqual(await page.evaluate(() => window.qrnamapApp.state.preview), { n: 100005, label: 'preview of the first 100,005 pairs' });
  assert.deepEqual(errors, []);
});

test('a sample where no read pseudoaligns gives a clear message', async (t) => {
  const { page, errors } = await open(t);
  const reads = (m) => Array.from({ length: 300 }, (_, i) => `@z_${i}/${m}\n${'ACGT'.repeat(25)}\n+\n${'I'.repeat(100)}\n`).join('');
  await fillRow(page, 0, { name: 'Zero', condition: 'x', mate1: fileArg('z_R1.fastq', Buffer.from(reads(1))), mate2: fileArg('z_R2.fastq', Buffer.from(reads(2))) });
  const [s] = await run(page);
  assert.equal(s.status, 'error');
  assert.equal(s.kind, 'input');
  assert.match(s.error, /No read pseudoaligned to the transcripts of this genome \(kallisto: "zero reads pseudoaligned"\)/);
  assert.deepEqual(errors, []);
});

// ------------------------------------------------------------------ all reads, streamed from the File objects
test('all reads (streamed): three samples (gzip, plain, single-end) give the gene counts of the CLI', async (t) => {
  const { page, errors } = await open(t);
  assert.equal(await page.isVisible('#preview-pairs-p'), true);
  await page.check('#depth-all');
  assert.equal(await page.isVisible('#preview-pairs-p'), false, 'the N field is hidden for all reads');
  await fillRow(page, 0, { name: 'S1', condition: 'ctl', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await fillRow(page, 1, { name: 'S2', condition: 'trt', mate1: 's2_R1.fastq', mate2: 's2_R2.fastq' });
  await fillRow(page, 2, { name: 'S3', condition: 'trt', mate1: 's3.fastq.gz' });
  await page.fill('#frag-mean', '200');
  await page.fill('#frag-sd', '20');
  const samples = await run(page);
  assert.equal(await noErrors(page), true);
  assertCounts(samples[0], 's1');
  assertCounts(samples[1], 's2');
  assertCounts(samples[2], 's3');
  assert.deepEqual(samples.map((s) => [s.allReads, s.allPairs, s.previewN]), [[true, true, 0], [true, true, 0], [true, true, 0]]);
  // kallisto reads the files itself (gzip included): the names keep the .gz, and no file was copied into the page.
  assert.deepEqual(samples[0].args, ['quant', '-i', 'transcripts.idx', '-o', 'kq', '--plaintext', '-t', '1', 'r1.fastq.gz', 'r2.fastq.gz']);
  assert.deepEqual(samples[1].args, ['quant', '-i', 'transcripts.idx', '-o', 'kq', '--plaintext', '-t', '1', 'r1.fastq', 'r2.fastq']);
  assert.match(await row(page, 0).locator('.status').textContent(), /^300 pairs processed, .* All pairs \(300\), streamed from the files\.$/);
  assert.deepEqual(await page.evaluate(() => window.qrnamapApp.state.preview), { n: null, full: true, label: 'all reads' });
  assert.match(await page.textContent('#preview-label'), /all reads of each sample/);
  // A full-depth count table is not named a preview
  const [dl] = await Promise.all([page.waitForEvent('download'), page.click('#dl-counts')]);
  assert.equal(dl.suggestedFilename(), 'counts.full.tsv');
  assert.deepEqual(errors, []);
});

test('all reads (streamed): mates with different record counts and a cut gzip give clear errors', async (t) => {
  const { page, errors } = await open(t);
  await page.check('#depth-all');
  const lines = fixture('s2_R2.fastq').toString().split('\n');
  const short = { name: 'short_R2.fastq', mimeType: 'text/plain', buffer: Buffer.from(lines.slice(0, 200 * 4).join('\n') + '\n') };
  const gz = fixture('s1_R2.fastq.gz');
  const cut = { name: 'cut_R2.fastq.gz', mimeType: 'application/gzip', buffer: gz.subarray(0, gz.length - 20) };
  await fillRow(page, 0, { name: 'A', condition: 'x', mate1: 's2_R1.fastq', mate2: short });
  const b1 = { name: 'b_R1.fastq.gz', mimeType: 'application/gzip', buffer: fixture('s1_R1.fastq.gz') };
  await fillRow(page, 1, { name: 'B', condition: 'x', mate1: b1, mate2: cut });
  await fillRow(page, 2, { name: 'C', condition: 'y', mate1: 's1_R1.fastq.gz', mate2: 's1_R2.fastq.gz' });
  await watchStatus(page);
  await page.click('#run-samples');
  assert.equal(await noErrors(page), true, await page.textContent('#samples-errors'));
  await waitOrDump(page, () => ['done', 'cancelled', 'lost'].includes(window.qrnamapApp.state.samplesPhase), 'all-reads error run', errors);
  const samples = await page.evaluate(() => window.qrnamapApp.state.samples);
  assert.equal(samples[0].status, 'error');
  assert.match(samples[0].error, /different record counts: mate 1 \(s2_R1\.fastq\) has 300, mate 2 \(short_R2\.fastq\) has 200/);
  assert.equal(samples[1].status, 'error');
  assert.match(samples[1].error, /cut_R2\.fastq\.gz.* (not a complete gzip file|not a valid gzip file)/);
  assertCounts(samples[2], 's1');
  assert.deepEqual(errors, []);
});
