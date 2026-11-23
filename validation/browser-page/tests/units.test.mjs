// Unit tests (Node, no browser) of the page modules that need no browser: the abundance.tsv check, the sheet helpers and
// the FASTQ cut. node --test validation/browser-page/tests/units.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync } from 'node:zlib';
import { checkAbundance } from '../../../site/analyze/abundance-check.js';
import { baseName, pageRules, filesByName, sheetText } from '../../../site/analyze/sheet.js';
import { firstRecords, checkFastq, firstNameStem, firstReadLength } from '../../../site/analyze/fastq.js';

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

test('previewLimit: 200,000 pairs up to 101 bases, 100,000 above; never above the hard limit', async () => {
  const { previewLimit, PREVIEW_TABLE, HARD_LIMIT } = await import('../../../site/analyze/preview-limit.js');
  assert.equal(previewLimit(50), 200000);
  assert.equal(previewLimit(101), 200000);
  assert.equal(previewLimit(102), 100000);
  assert.equal(previewLimit(150), 100000);
  assert.equal(previewLimit(250), 100000);
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
