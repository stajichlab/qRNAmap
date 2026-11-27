// Unit tests (Node, no browser) of the page modules that need no browser: the abundance.tsv check, the sheet helpers and
// the FASTQ cut. node --test validation/browser-page/tests/units.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { checkAbundance } from '../../../site/analyze/abundance-check.js';
import { baseName, pageRules, filesByName, sheetText } from '../../../site/analyze/sheet.js';
import { firstRecords, checkFastq, firstNameStem, firstReadLength, streamFastqRecords } from '../../../site/analyze/fastq.js';
import { StreamRingBuffer } from '../../../site/analyze/stream-ring.js';
import { geneExternalLinks, genePopupHtml } from '../../../site/analyze/gene-links.js';

const GOOD = 'target_id\tlength\teff_length\test_counts\ttpm\n' +
  'tx1\t1491\t1327.5\t9\t31.1145\n' +
  'tx2\t800\t636.2\t40\t204.998\n' +
  'tx3\t75\t26.5\t0.5\t10.2\n' +
  'tx4\t300\t140.1\t0\t0\n';
const OK = { nTranscripts: 4, nPseudoaligned: 50 }; // 9 + 40 + 0.5 = 49.5: within the rounding

test('checkAbundance: a good file passes', () => {
  assert.deepEqual(checkAbundance(GOOD, OK), []);
});

test('checkAbundance: `9` written as `999` (the measured Firefox fault) fails the sum check', () => {
  const p = checkAbundance(GOOD.replace('\t9\t31.1145', '\t999\t31.1145'), OK);
  assert.deepEqual(p, ['the est_counts add up to 1039.500, but run_info.json says 50 reads pseudoaligned']);
});

test('checkAbundance: a repeated decimal token fails the number check', () => {
  const p = checkAbundance(GOOD.replace('\t204.998\n', '\t204.998204.998204.998\n'), OK);
  assert.equal(p.length, 1);
  assert.match(p[0], /^line 3 is not a name and 4 numbers/);
});

test('checkAbundance: a missing line, a bad header and a wrong sum are named', () => {
  const lines = GOOD.split('\n');
  assert.deepEqual(checkAbundance([lines[0], lines[1], lines[2], lines[4], ''].join('\n'), OK), ['it has 3 transcript lines, the index has 4 transcripts']);
  assert.deepEqual(checkAbundance(GOOD.replace('est_counts', 'counts'), OK), ['the first line is not the kallisto header']);
  assert.deepEqual(checkAbundance(GOOD, { nTranscripts: 4, nPseudoaligned: 60 }), ['the est_counts add up to 49.500, but run_info.json says 60 reads pseudoaligned']);
  assert.deepEqual(checkAbundance(GOOD, { nTranscripts: 4, nPseudoaligned: undefined }), ['run_info.json has no n_pseudoaligned number']);
});

test('checkAbundance: the tolerance follows the printed digits; nan is allowed only for tpm', () => {
  // 6 significant digits: 123457 can be 123456.5 to 123457.5, so 0.5 of rounding per value, plus 1.
  const big = 'target_id\tlength\teff_length\test_counts\ttpm\na\t1\t1\t123457\t1\nb\t1\t1\t123457\t1\n';
  assert.deepEqual(checkAbundance(big, { nTranscripts: 2, nPseudoaligned: 246916 }), []);
  assert.equal(checkAbundance(big, { nTranscripts: 2, nPseudoaligned: 246911 }).length, 1);
  const zero = 'target_id\tlength\teff_length\test_counts\ttpm\na\t1\t1\t0\t-nan\n';
  assert.deepEqual(checkAbundance(zero, { nTranscripts: 1, nPseudoaligned: 0 }), []);
  assert.equal(checkAbundance(zero.replace('\t0\t-nan', '\tnan\t0'), { nTranscripts: 1, nPseudoaligned: 0 }).length, 1);
});

test('sheet helpers: base names, lanes, a name used twice, duplicate picked files, the sheet text', () => {
  assert.equal(baseName('/data/run 1/a_R1.fastq.gz'), 'a_R1.fastq.gz');
  assert.equal(baseName('C:\\reads\\a.fq'), 'a.fq');
  assert.equal(baseName('a.fq'), 'a.fq');
  assert.deepEqual(pageRules([{ sample: 'A', fastq1: ['x/a.fq'], fastq2: ['b.fq'] }, { sample: 'B', fastq1: ['c.fq'], fastq2: [] }]), []);
  const lanes = pageRules([{ sample: 'A', fastq1: ['a1.fq', 'a2.fq'], fastq2: ['b1.fq', 'b2.fq'] }]);
  assert.equal(lanes.length, 2);
  assert.match(lanes[0], /sample "A": the fastq1 cell lists 2 files \(lanes\)/);
  assert.deepEqual(pageRules([{ sample: 'A', fastq1: ['x/a.fq'], fastq2: ['y/a.fq'] }]),
    ['the file name "a.fq" is used twice (sample "A" fastq1 and sample "A" fastq2). The page matches files by name only, so each file name can be used once.']);
  const { byName, errors } = filesByName([{ name: 'a.fq' }, { name: 'b.fq' }, { name: 'a.fq' }]);
  assert.deepEqual([...byName.keys()], ['a.fq', 'b.fq']);
  assert.equal(errors.length, 1);
  assert.equal(sheetText(['condition', 'empty', 'dose:numeric'], [
    { name: 'S1', values: ['ctl', '', '1'], fastq1: 's1_1.fq', fastq2: 's1_2.fq' },
    { name: 'S2', values: [' trt ', '', '2'], fastq1: 's2.fq', fastq2: '' },
  ]), 'sample\tfastq1\tfastq2\tcondition\tdose:numeric\nS1\ts1_1.fq\ts1_2.fq\tctl\t1\nS2\ts2.fq\t\ttrt\t2\n');
});

const fastq = (n, from = 0) => Array.from({ length: n }, (_, i) => `@r${from + i}\nACGT\n+\nIIII\n`).join('');

test('firstRecords: cuts plain and gzip files at N records, keeps all of a short file, and names a bad end', async () => {
  const text = fastq(10);
  let r = await firstRecords(new Blob([text]), 4, 'p');
  assert.deepEqual([r.records, r.wholeFile, new TextDecoder().decode(r.bytes)], [4, false, fastq(4)]);
  // Two gzip members (as bgzip writes), cut inside the second member.
  const two = new Blob([gzipSync(fastq(3)), gzipSync(fastq(7, 3))]);
  r = await firstRecords(two, 6, 'g');
  assert.deepEqual([r.records, r.wholeFile, new TextDecoder().decode(r.bytes)], [6, false, fastq(6)]);
  r = await firstRecords(two, 100, 'g');
  assert.deepEqual([r.records, r.wholeFile], [10, true]);
  // No newline at the end of the file: the last line still counts.
  r = await firstRecords(new Blob([text.slice(0, -1)]), 100, 'p');
  assert.deepEqual([r.records, r.wholeFile, r.bytes.length], [10, true, text.length]);
  await assert.rejects(firstRecords(new Blob([fastq(2) + '@x\nAC\n']), 100, 'Mate 1'), /Mate 1 ends inside a record: it has 10 lines, which is not a multiple of 4/);
  const gz = gzipSync(text);
  await assert.rejects(firstRecords(new Blob([gz.subarray(0, gz.length - 5)]), 100, 'Mate 2'), /Mate 2 is not a complete gzip file/);
  await assert.rejects(firstRecords(new Blob([]), 1, 'Mate 1'), /Mate 1 is empty \(0 bytes\)/);
  // An aborted read stops.
  const controller = new AbortController();
  controller.abort();
  await assert.rejects(firstRecords(new Blob([text]), 4, 'p', { signal: controller.signal }), { name: 'AbortError' });
});

test('checkFastq: names the first bad record', () => {
  const enc = (s) => new TextEncoder().encode(s);
  checkFastq(enc(fastq(3)), 'm');
  checkFastq(enc(fastq(2).replaceAll('\n', '\r\n')), 'm');
  assert.throws(() => checkFastq(enc(fastq(1) + '>r\nACGT\n+\nIIII\n'), 'm'), /m is not a FASTQ file: record 2 \(line 5\) does not start with "@"/);
  assert.throws(() => checkFastq(enc(fastq(1) + '@r\nACGT\n-\nIIII\n'), 'm'), /record 2 \(line 7\) has no "\+" line/);
  assert.throws(() => checkFastq(enc('@r\nACGT\n+\nIII\n'), 'm'), /record 1 \(line 4\) has 4 bases but 3 quality letters/);
});

test('firstRecords: a blank line gives its own message (not the line-count message)', async () => {
  await assert.rejects(firstRecords(new Blob([fastq(2) + '\n']), 100, 'Mate 1'), /Mate 1 has a blank line \(line 9\)\. A FASTQ file has no blank lines/);
  await assert.rejects(firstRecords(new Blob(['\n' + fastq(2)]), 100, 'Mate 1'), /Mate 1 has a blank line \(line 1\)/);
  await assert.rejects(firstRecords(new Blob([gzipSync(fastq(2).replaceAll('\n', '\r\n') + '\r\n')]), 100, 'Mate 2'), /Mate 2 has a blank line \(line 9\)/);
  // A blank line after the cut is not read, so it does not matter.
  const r = await firstRecords(new Blob([fastq(2) + '\n']), 2, 'Mate 1');
  assert.equal(r.records, 2);
});

test('firstNameStem: first word, without @ and without /1 or /2', () => {
  const enc = (s) => new TextEncoder().encode(s);
  assert.equal(firstNameStem(enc('@r1/1\nA\n+\nI\n')), 'r1');
  assert.equal(firstNameStem(enc('@SRR1.1 1 length=101\nA\n+\nI\n')), 'SRR1.1');
  assert.equal(firstNameStem(enc('@x/2\r\nA\r\n')), 'x');
});

test('previewLimit: 200,000 pairs up to 101 bases, 75,000 above; never above the hard limit', async () => {
  const { previewLimit, PREVIEW_TABLE, HARD_LIMIT } = await import('../../../site/analyze/preview-limit.js');
  assert.equal(previewLimit(50), 200000);
  assert.equal(previewLimit(101), 200000);
  assert.equal(previewLimit(102), 75000);
  assert.equal(previewLimit(150), 75000);
  assert.equal(previewLimit(250), 75000);
  assert.throws(() => previewLimit(0), RangeError);
  assert.throws(() => previewLimit(undefined), RangeError);
  assert.ok(PREVIEW_TABLE.every((r) => r.pairs <= HARD_LIMIT));
  assert.ok(PREVIEW_TABLE.every((r) => !/PROVISIONAL/.test(r.measured)));
});

test('firstReadLength: the sequence length of the first record (CRLF ok)', () => {
  const enc = (s) => new TextEncoder().encode(s);
  assert.equal(firstReadLength(enc('@r\nACGTA\n+\nIIIII\n@s\nAC\n+\nII\n')), 5);
  assert.equal(firstReadLength(enc('@r\r\nACG\r\n+\r\nIII\r\n')), 3);
  assert.equal(firstReadLength(enc('@r\n')), 0);
});

test('geneExternalLinks: SGD first, then NCBI Gene (RefSeq yeast Dbxref)', () => {
  assert.deepEqual(geneExternalLinks('Q0130', { geneId: '854584', dbxref: ['GeneID:854584', 'SGD:S000007274'] }), [
    { name: 'SGD', url: 'https://www.yeastgenome.org/locus/S000007274' },
    { name: 'NCBI Gene', url: 'https://www.ncbi.nlm.nih.gov/gene/854584' },
  ]);
});

test('geneExternalLinks: FungiDb xref before NCBI Gene (Af293-style RefSeq Dbxref)', () => {
  assert.deepEqual(geneExternalLinks('AFUA_1G00100', { geneId: '3507995', dbxref: ['FungiDB:AFUA_1G00100', 'GeneID:3507995'] }), [
    { name: 'FungiDB', url: 'https://fungidb.org/gene/AFUA_1G00100' },
    { name: 'NCBI Gene', url: 'https://www.ncbi.nlm.nih.gov/gene/3507995' },
  ]);
});

test('geneExternalLinks: an AFU locus tag links FungiDB without a Dbxref (A1163 GenBank)', () => {
  assert.deepEqual(geneExternalLinks('AFUB_102960', { geneId: '', dbxref: [] }), [
    { name: 'FungiDB', url: 'https://fungidb.org/gene/AFUB_102960' },
  ]);
  // the xref and the label give the same URL once, not twice
  assert.deepEqual(geneExternalLinks('AFUA_1G00100', { dbxref: ['FungiDB:AFUA_1G00100'] }), [
    { name: 'FungiDB', url: 'https://fungidb.org/gene/AFUA_1G00100' },
  ]);
});

test('geneExternalLinks: SGD beats FungiDB when both are present', () => {
  assert.equal(geneExternalLinks('Q0130', { dbxref: ['FungiDB:S000007274', 'SGD:S000007274'] })[0].name, 'SGD');
});

test('geneExternalLinks: no dbxref but a numeric gene_id still links NCBI; nothing else links nothing', () => {
  assert.deepEqual(geneExternalLinks('YAL068C', { geneId: '854598', dbxref: [] }), [
    { name: 'NCBI Gene', url: 'https://www.ncbi.nlm.nih.gov/gene/854598' },
  ]);
  assert.deepEqual(geneExternalLinks('YAL068C', { geneId: '', dbxref: [] }), []);
  assert.deepEqual(geneExternalLinks('tP(UGG)Q', { geneId: '', dbxref: ['FungiDB:', 'SGD:'] }), []);
  assert.deepEqual(geneExternalLinks('anything', undefined), []);
});

test('StreamRingBuffer and streamFastqRecords: streams gzip and plain FASTQ records through ring buffer', async () => {
  const { Worker } = await import('node:worker_threads');
  const text = fastq(20);
  const ring = new StreamRingBuffer(null, 16384); // 16 KB ring buffer
  const gzBlob = new Blob([gzipSync(text)]);

  const worker = new Worker(`
    import { parentPort, workerData } from 'node:worker_threads';
    import { StreamRingBuffer } from './site/analyze/stream-ring.js';
    const ring = new StreamRingBuffer(workerData.sab);
    const buf = new Uint8Array(1024);
    const chunks = [];
    while (true) {
      const n = ring.read(buf, 0, buf.length);
      if (n === 0) break;
      chunks.push(buf.slice(0, n));
    }
    const total = chunks.reduce((a, b) => a + b.length, 0);
    const joined = new Uint8Array(total);
    let off = 0;
    for (const c of chunks) { joined.set(c, off); off += c.length; }
    parentPort.postMessage({ total, text: new TextDecoder().decode(joined) });
  `, { eval: true, type: 'module', workerData: { sab: ring.sab } });

  const [streamRes, workerRes] = await Promise.all([
    streamFastqRecords(gzBlob, 15, 'test-mate', ring),
    new Promise(res => worker.on('message', res))
  ]);

  assert.equal(streamRes.records, 15);
  assert.equal(streamRes.wholeFile, false);
  assert.equal(workerRes.text, fastq(15));
});


test('mountBlobFile: a File-backed file in the real kallisto.mjs FS has the blob size and its bytes, read through a small window', async () => {
  const { mountBlobFile } = await import('../../../site/vendor/kallisto/runtime/execute.mjs');
  const factory = (await import('../../../site/vendor/kallisto/kallisto.mjs')).default;
  const data = new Uint8Array(100_003).map((_, i) => (i * 7 + 3) & 255);
  const blob = new Blob([data]);
  const reads = [];
  // Node has no FileReaderSync: the test reader takes the bytes from `data` and records each window read.
  const readSync = (b, start, end) => { assert.equal(b, blob); reads.push([start, end]); return data.slice(start, end); };
  const m = await factory({ noInitialRun: true, print: () => {}, printErr: () => {} });
  m.FS.mkdirTree('/work');
  mountBlobFile(m.FS, '/work/r1.fastq.gz', blob, { readSync, window: 4096 });
  assert.equal(m.FS.stat('/work/r1.fastq.gz').size, data.length);
  assert.ok(m.FS.isFile(m.FS.stat('/work/r1.fastq.gz').mode));
  assert.deepEqual(m.FS.readFile('/work/r1.fastq.gz'), data);
  assert.ok(reads.every(([s, e]) => e - s <= Math.max(4096, data.length)), 'each read asks for a window, not the file');
  // Small sequential reads and a seek back are served from the window; a read past the end gives 0 bytes.
  const s = m.FS.open('/work/r1.fastq.gz', 'r');
  const buf = new Uint8Array(10);
  reads.length = 0;
  assert.equal(m.FS.read(s, buf, 0, 10), 10);
  assert.equal(m.FS.read(s, buf, 0, 10), 10);
  assert.deepEqual(buf, data.slice(10, 20));
  assert.ok(reads.length <= 1, 'two small reads, at most one window read');
  m.FS.llseek(s, 0, 2);
  assert.equal(m.FS.read(s, buf, 0, 10), 0);
  m.FS.close(s);
  assert.throws(() => m.FS.writeFile('/work/r1.fastq.gz', 'x'));
});

// In-page differential expression (site/analyze/de-stats.js).
import { computeDifferentialExpression, tTwoSided, bhAdjust } from '../../../site/analyze/de-stats.js';

/** A matrix of nGenes genes, 3 control + 3 treated samples, with gene 0..nUp-1 up fourfold in treated. */
function deMatrix({ nGenes = 400, nUp = 20, depth = 1, seed = 7 } = {}) {
  let s = seed;
  const rnd = () => ((s = (s * 1103515245 + 12345) % 2147483648) / 2147483648);
  const names = ['c1', 'c2', 'c3', 't1', 't2', 't3'];
  const lines = ['gene\t' + names.join('\t')];
  for (let g = 0; g < nGenes; g++) {
    const mu = depth * (50 + 2000 * rnd());
    const vals = names.map((n) => {
      const fold = n.startsWith('t') && g < nUp ? 4 : 1;
      return Math.round(mu * fold * Math.exp(0.15 * (rnd() + rnd() + rnd() - 1.5))); // ~15% biological noise
    });
    lines.push(`g${g}\t${vals.join('\t')}`);
  }
  const samples = names.map((name) => ({ name, covariates: { condition: name[0] === 't' ? 'trt' : 'ctl' } }));
  return { tsv: lines.join('\n'), samples, contrast: { factor: 'condition', level: 'trt', reference: 'ctl' } };
}

test('tTwoSided and bhAdjust match reference values', () => {
  assert.ok(Math.abs(tTwoSided(2.776445, 4) - 0.05) < 1e-4); // qt(0.975, 4)
  assert.ok(Math.abs(tTwoSided(1.959964, 1e6) - 0.05) < 1e-4);
  assert.equal(tTwoSided(0, 5), 1);
  assert.deepEqual(bhAdjust([0.01, 0.04, 0.03, 0.5]).map((x) => +x.toFixed(4)), [0.04, 0.0533, 0.0533, 0.5]); // p.adjust(..., "BH")
});

test('computeDifferentialExpression: a fourfold change is found at any depth (the SE is on the log2 scale)', () => {
  for (const depth of [1, 100, 10000]) { // 10000x is a deep run: the old raw-count SE found nothing here
    const { tsv, samples, contrast } = deMatrix({ depth });
    const res = computeDifferentialExpression(tsv, samples, contrast);
    const up = res.slice(0, 20);
    assert.ok(up.every((r) => r.padj < 0.05 && Math.abs(r.log2FoldChange - 2) < 0.4), `depth ${depth}: all 20 up genes significant`);
    const falsePos = res.slice(20).filter((r) => r.padj < 0.05).length;
    assert.ok(falsePos <= 5, `depth ${depth}: ${falsePos} false positives of 380`);
  }
});

test('computeDifferentialExpression: no replicates means no test (p = 1); missing groups give no rows', () => {
  const tsv = 'gene\ta\tb\n' + Array.from({ length: 20 }, (_, i) => `g${i}\t${10 + i}\t${40 + i}`).join('\n');
  const samples = [{ name: 'a', covariates: { condition: 'ctl' } }, { name: 'b', covariates: { condition: 'trt' } }];
  const res = computeDifferentialExpression(tsv, samples, { factor: 'condition', level: 'trt', reference: 'ctl' });
  assert.equal(res.length, 20);
  assert.ok(res.every((r) => r.pvalue === 1 && r.padj === 1));
  assert.deepEqual(computeDifferentialExpression(tsv, samples, { factor: 'condition', level: 'x', reference: 'ctl' }), []);
});

test('genePopupHtml: function text, IDs and the volcano link-outs, all escaped', () => {
  const html = genePopupHtml('AFUA_2G00200', {
    symbol: 'catA<x>', geneId: '3510001', description: 'catalase A & "B"', dbxref: ['FungiDB:AFUA_2G00200'], contig: 'NC_1', strand: '+',
  });
  assert.match(html, /<h4>AFUA_2G00200 <span class="small">\(catA&lt;x&gt;\)<\/span><\/h4>/);
  assert.match(html, /<p>catalase A &amp; &quot;B&quot;<\/p>/);
  assert.match(html, /<dt>NCBI GeneID<\/dt><dd>3510001<\/dd>/);
  assert.match(html, /<dt>Location<\/dt><dd>NC_1 \(\+\)<\/dd>/);
  const hrefs = [...html.matchAll(/href="([^"]+)" target="_blank" rel="noopener noreferrer"/g)].map((m) => m[1]);
  assert.deepEqual(hrefs, geneExternalLinks('AFUA_2G00200', { geneId: '3510001', dbxref: ['FungiDB:AFUA_2G00200'] }).map((l) => l.url));
  const bare = genePopupHtml('<g>', undefined);
  assert.match(bare, /<h4>&lt;g&gt;<\/h4>/);
  assert.match(bare, /No function text/);
  assert.match(bare, /No database link/);
});
