// The page worker (a module worker). After a genome is ready it holds the kallisto index (with the genome label and
// accession). It never sends the index to the page. It also reads the preview reads of one sample at a time and runs
// kallisto quant on them (quant.js sends the requests). The message protocol is in site/analyze/PROTOCOL.md.
// kallisto runs through the vendored ToolRunner (vendor/kallisto/runtime/index.mjs), which starts one more worker per
// command and terminates it at the end, so a crash of kallisto does not end this worker.
import init, { kallisto_args } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { prepareTranscripts, StepError } from './prep.js?v=__BUILD__';
import { VerifyingGunzip } from '../gz.js?v=__BUILD__';
import { ToolRunner } from '../vendor/kallisto/runtime/index.mjs?v=__BUILD__';
import { firstRecords, checkFastq, firstNameStem, firstReadLength, streamFastqRecords, countRecords, ReadError } from './fastq.js?v=__BUILD__';
import { StreamRingBuffer } from './stream-ring.js?v=__BUILD__';
import { previewLimit, HARD_LIMIT } from './preview-limit.js?v=__BUILD__';
import { checkAbundance } from './abundance-check.js?v=__BUILD__';

const ready = init({ module_or_path: new URL('../pkg-web/qrnamap_wasm_bg.wasm?v=__BUILD__', import.meta.url) });
// Stamped like the page files, so a cached old kallisto cannot meet a new page. The runtime worker imports execute.mjs without a stamp
// (a vendored file, not changed).
const KALLISTO = new URL('../vendor/kallisto/kallisto.mjs?v=__BUILD__', import.meta.url);
const RUNTIME_WORKER = new URL('../vendor/kallisto/runtime/worker.mjs?v=__BUILD__', import.meta.url);
/** The genome of the page, or null. */
let genome = null;

const post = (message) => self.postMessage(message);

function memoryLike(text) {
  return /out of memory|memory access out of bounds|could not allocate|cannot enlarge memory|allocation failed|Array buffer allocation|RangeError|OOM/i.test(text);
}

/** The kallisto index argument list, from qrnamap-core (the same list as the Task 3 measurements). */
const indexArgs = () => JSON.parse(kallisto_args(JSON.stringify({ kind: 'index', fasta: 'transcripts.fa', index: 'transcripts.idx' })));

/**
 * Read the genome files and build the transcripts in a short-lived worker and end it, so the memory of the genome bytes is
 * released before the kallisto index is built (default; `?prep_worker=0` builds them in this worker, for measurements). In
 * WebKit this lowered the peak of H99 from 1,586,569 to 1,491,448 kB (docs/browser/page-v1-measurements.md).
 */
function prepareInWorker(message, progress) {
  return new Promise((resolve, reject) => {
    const prep = new Worker(new URL('./prep-worker.js?v=__BUILD__', import.meta.url), { type: 'module' });
    const done = (fn, value) => { prep.terminate(); fn(value); };
    prep.onmessage = ({ data }) => {
      if (data.type === 'progress') progress(data.phase, data.text, data.done, data.total);
      else if (data.type === 'prepared') { data.built.transcripts = new Uint8Array(data.built.transcripts.buffer || data.built.transcripts); done(resolve, data.built); }
      else done(reject, new StepError(data.kind, data.message));
    };
    prep.onerror = (event) => done(reject, new StepError('internal', `The prep worker failed: ${event.message || 'unknown error'}`));
    const { fasta, gff, limitBases, limitText } = message;
    message.fasta = null;
    message.gff = null;
    prep.postMessage({ fasta, gff, limitBases, limitText });
  });
}

async function buildGenome(message) {
  const { id, label, accession } = message;
  genome = null; // a new genome replaces the old one at once, so no step can use a stale index
  const progress = (phase, text, done = 0, total = 0) => post({ type: 'progress', id, phase, text, done, total });
  const built = message.prepWorker ? await prepareInWorker(message, progress) : await prepareTranscripts(message, progress);

  progress('index', `Building the kallisto index of ${built.nTranscripts.toLocaleString('en-US')} transcripts. This can take one to two minutes for a large genome.`);
  const args = indexArgs();
  const runner = new ToolRunner({ moduleUrl: KALLISTO, workerUrl: RUNTIME_WORKER, waitForExit: false });
  let run;
  const started = performance.now();
  // The desktop app keeps the transcripts: its native kallisto may not read this kallisto.wasm index and then builds its own.
  const kept = message.keepTranscripts ? built.transcripts : null;
  try {
    run = await runner.run(args, { files: [{ path: 'transcripts.fa', data: built.transcripts }], outputs: ['transcripts.idx'] });
  } catch (error) {
    const text = `${error?.name || 'Error'}: ${error?.message || error}`;
    throw new StepError(memoryLike(text) ? 'memory' : 'kallisto', `kallisto index stopped (${text}).`);
  } finally {
    built.transcripts = null; // not kept: Task 6 needs only the index (PROTOCOL.md)
  }
  const stderr = new TextDecoder().decode(run.stderr);
  if (run.exitCode !== 0 || run.missingOutputs.length) {
    const tail = stderr.trim().split('\n').slice(-3).join(' ');
    throw new StepError(memoryLike(stderr) ? 'memory' : 'kallisto', `kallisto index failed (exit ${run.exitCode}${run.missingOutputs.length ? ', no index file' : ''}). ${tail}`);
  }
  const index = run.files.find((f) => f.path === '/work/transcripts.idx').data;
  // The worker keeps only the index (and, in the desktop app, the transcripts). tx2gene and the gene table go to the page, which keeps them.
  genome = { label, accession, index, transcripts: kept };
  post({ type: 'genome-ready', id, label, accession, nTranscripts: built.nTranscripts, nGenes: built.nGenes, bases: built.bases,
    indexBytes: index.length, indexSeconds: (performance.now() - started) / 1000, indexArgs: args, tx2gene: built.tx2gene, genesTsv: built.genesTsv, geneInfoTsv: built.geneInfoTsv });
}

function decompressGzip(bytes) {
  return new Promise((resolve, reject) => {
    const chunks = [];
    try {
      const gunzip = new VerifyingGunzip((chunk) => chunks.push(chunk));
      gunzip.push(new Uint8Array(bytes), true);
      const total = chunks.reduce((acc, c) => acc + c.length, 0);
      const out = new Uint8Array(total);
      let off = 0;
      for (const c of chunks) {
        out.set(c, off);
        off += c.length;
      }
      resolve(out);
    } catch (err) {
      reject(err);
    }
  });
}

/** Load a pre-computed kallisto index and annotation tables directly into memory. */
async function loadPrebuiltGenome(message) {
  const { id, label, accession, idxGz, tx2geneGz, geneInfoGz, genesGz, bases, nTranscripts, nGenes } = message;
  genome = null;
  const started = performance.now();

  const [index, tx2geneBytes, geneInfoBytes, genesBytes] = await Promise.all([
    decompressGzip(idxGz),
    decompressGzip(tx2geneGz),
    decompressGzip(geneInfoGz),
    decompressGzip(genesGz),
  ]);

  const dec = new TextDecoder();
  const tx2gene = dec.decode(tx2geneBytes);
  const geneInfoTsv = dec.decode(geneInfoBytes);
  const genesTsv = dec.decode(genesBytes);

  const calcTranscripts = nTranscripts || (tx2gene.trim().split('\n').length);
  const calcGenes = nGenes || (genesTsv.trim().split('\n').length - 1);

  genome = { label, accession, index };
  const args = indexArgs();
  post({
    type: 'genome-ready',
    id,
    label,
    accession,
    nTranscripts: calcTranscripts,
    nGenes: calcGenes,
    bases: bases || 0,
    indexBytes: index.length,
    indexSeconds: (performance.now() - started) / 1000,
    indexArgs: args,
    tx2gene,
    genesTsv,
    geneInfoTsv,
  });
}

// ------------------------------------------------------------------ quantification of one sample

const fmt = (n) => Number(n).toLocaleString('en-US');

/**
 * One kallisto quant run. Gives { abundance (text), info (run_info.json numbers, or null when the file cannot be read),
 * infoError }. `files` is emptied as soon as ToolRunner has the run, so only ToolRunner refers to the reads. `signal` stops it.
 */
async function kallistoQuant(args, files, signal, devices = []) {
  const runner = new ToolRunner({ moduleUrl: KALLISTO, workerUrl: RUNTIME_WORKER, waitForExit: false });
  let run;
  try {
    const running = runner.run(args, { files: [...files], devices: [...devices], outputs: ['kq/abundance.tsv', 'kq/run_info.json'], signal });
    files.length = 0; // ToolRunner has posted (copied) the files; it keeps its own reference until the run returns
    run = await running;
  } catch (error) {
    if (error?.name === 'AbortError') throw error;
    const text = `${error?.name || 'Error'}: ${error?.message || error}`;
    throw new StepError(memoryLike(text) ? 'memory' : 'kallisto', `kallisto quant stopped (${text}).`);
  }
  const stderr = new TextDecoder().decode(run.stderr);
  if (run.exitCode !== 0 && /zero reads pseudoaligned/.test(stderr)) {
    throw new StepError('input', 'No read pseudoaligned to the transcripts of this genome (kallisto: "zero reads pseudoaligned"). Check that the reads come from this organism and that they are RNA-seq reads.');
  }
  if (run.exitCode !== 0 || run.missingOutputs.length) {
    const tail = stderr.trim().split('\n').slice(-3).join(' ');
    throw new StepError(memoryLike(stderr) ? 'memory' : 'kallisto', `kallisto quant failed (exit ${run.exitCode}${run.missingOutputs.length ? ', an output file is missing' : ''}). ${tail}`);
  }
  const text = (path) => new TextDecoder().decode(run.files.find((f) => f.path === `/work/${path}`).data);
  let info = null;
  let infoError = '';
  try {
    const raw = JSON.parse(text('kq/run_info.json'));
    info = { n_processed: raw.n_processed, n_pseudoaligned: raw.n_pseudoaligned, p_pseudoaligned: raw.p_pseudoaligned };
  } catch (error) {
    infoError = `run_info.json cannot be read (${error.message})`;
  }
  return { abundance: text('kq/abundance.tsv'), info, infoError };
}

/** Problems of one kallisto output: run_info.json, n_processed, and the abundance.tsv check. Empty when the output can be used. */
function outputProblems(result, nPairs, nTranscripts, unit) {
  if (!result.info) return [result.infoError];
  const problems = [];
  if (result.info.n_processed !== nPairs) problems.push(`kallisto processed ${fmt(result.info.n_processed)} ${unit}, but the page gave it ${fmt(nPairs)}`);
  return problems.concat(checkAbundance(result.abundance, { nTranscripts, nPseudoaligned: result.info.n_pseudoaligned }));
}

/**
 * Read the first `pairs` records of each mate, check them, run kallisto quant with the index of the genome, check the
 * output (outputProblems; one more run when it fails), and answer `quant-done`. For the second run the files are read again
 * instead of keeping the arrays. Reasoning, not a measurement: the vendored ToolRunner keeps its input until the run returns,
 * so the peak during a run is the same either way; re-reading only means that no reads are held between the two runs.
 * Nothing of the sample is kept after the answer.
 */
async function quantSample(message, signal) {
  const { id, name, fragment, nTranscripts, limitOverride } = message;
  const { mate1, mate2 } = message;
  message.mate1 = null; // the request must not keep the files
  message.mate2 = null;
  if (!genome) throw new StepError('no-index', 'The kallisto index is not in the page worker any more. Build the genome index again (step 1).');
  const progress = (phase, text) => post({ type: 'progress', id, phase, text });
  // N for this sample: the read length (the longest first record of the two mates) gives the limit (preview-limit.js), or the
  // measurement option `limitOverride` replaces it (never above HARD_LIMIT); `message.pairs` (the page field) can lower it.
  progress('prep', `Reading the first record of each mate to find the read length`);
  let readLength = 0;
  let first1 = null;
  let first2 = null;
  for (const [file, mate] of [[mate1, 1], [mate2, 2]]) {
    if (!file) continue;
    const label = `Mate ${mate} of sample ${name} (${file.name})`;
    const first = await firstRecords(file, 1, label, { signal });
    checkFastq(first.bytes, label);
    readLength = Math.max(readLength, firstReadLength(first.bytes));
    if (mate === 1) first1 = first;
    else first2 = first;
  }
  if (readLength < 1) throw new StepError('input', `The first record of sample ${name} has no bases.`);
  let warning = '';
  if (first2) {
    const a = firstNameStem(first1.bytes);
    const b = firstNameStem(first2.bytes);
    if (a !== b) warning = `the first read names of the two mates differ ("${a}" and "${b}"). Check that the two files belong together.`;
  }
  const unit = mate2 ? 'pairs' : 'reads';
  if (message.allReads) return quantAllReads({ id, name, fragment, nTranscripts, mate1, mate2, readLength, warning, unit }, signal, progress);
  const limit = limitOverride ? Math.min(limitOverride, HARD_LIMIT) : previewLimit(readLength);
  const pairs = Math.min(message.pairs, limit);
  await ready;
  const args = JSON.parse(kallisto_args(JSON.stringify({
    kind: 'quant', sample: { name, fastq1: ['r1.fastq'], fastq2: mate2 ? ['r2.fastq'] : [] },
    index: 'transcripts.idx', out: 'kq', threads: 1, ...(mate2 ? {} : { fragment_len: fragment }),
  })));

  const canStream = typeof SharedArrayBuffer !== 'undefined';
  let result;
  let nPairs = 0;
  let wholeFile = false;
  let firstProblems = null;

  if (canStream) {
    for (let attempt = 1; attempt <= 2; attempt++) {
      if (attempt === 2) {
        progress('prep', `The kallisto output failed the check (${firstProblems.join('; ')}). Streaming the files again for one more run.`);
      }
      progress('quant', `Streaming kallisto quant on up to ${fmt(pairs)} ${unit}${attempt === 2 ? ' (second run)' : ''}`);
      const ring1 = new StreamRingBuffer();
      const ring2 = mate2 ? new StreamRingBuffer() : null;
      const devices = [{ path: 'r1.fastq', ring: ring1 }];
      if (ring2) devices.push({ path: 'r2.fastq', ring: ring2 });
      const files = [{ path: 'transcripts.idx', data: genome.index }];

      const quantPromise = kallistoQuant(args, files, signal, devices);
      const stream1 = streamFastqRecords(mate1, pairs, `Mate 1 of sample ${name} (${mate1.name})`, ring1, {
        signal,
        onRecords: (n) => progress('prep', `Streaming mate 1 (${mate1.name}): ${fmt(n)} ${unit}`)
      });
      const stream2 = ring2 ? streamFastqRecords(mate2, pairs, `Mate 2 of sample ${name} (${mate2.name})`, ring2, {
        signal,
        onRecords: (n) => progress('prep', `Streaming mate 2 (${mate2.name}): ${fmt(n)} ${unit}`)
      }) : Promise.resolve(null);

      let streamResults;
      try {
        [result, streamResults] = await Promise.all([
          quantPromise,
          Promise.all([stream1, stream2])
        ]);
      } catch (err) {
        ring1.abort();
        if (ring2) ring2.abort();
        throw err;
      }

      const [r1, r2] = streamResults;
      nPairs = r1.records;
      wholeFile = r1.wholeFile && (!r2 || r2.wholeFile);

      if (r2 && r1.records !== r2.records) {
        const count = (r) => `${fmt(r.records)} records${r.wholeFile ? '' : ` in its first ${fmt(pairs)}`}`;
        throw new StepError('input', `The mates of sample ${name} have different record counts: mate 1 (${mate1.name}) has ${count(r1)}, ` +
          `mate 2 (${mate2.name}) has ${count(r2)}. kallisto pairs the records by position and stops at the shorter file without an error, so the page stops here. ` +
          'Check that the two files belong together and are complete.');
      }

      const problems = outputProblems(result, nPairs, nTranscripts, unit);
      if (!problems.length) break;
      if (attempt === 2) {
        throw new StepError('output', `kallisto wrote a damaged output twice for sample ${name} (first run: ${firstProblems.join('; ')}; second run: ${problems.join('; ')}). ` +
          'The page does not use it. A damaged abundance.tsv was seen in Firefox in the measurements (docs/browser/page-v1-measurements.md). Try again, or use another browser or the command-line tool.');
      }
      firstProblems = problems;
    }
  } else {
    const read = async (file, mate) => {
      const label = `Mate ${mate} of sample ${name} (${file.name})`;
      progress('prep', `Reading the first ${fmt(pairs)} records of mate ${mate} (${file.name})`);
      const r = await firstRecords(file, pairs, label, { signal, onRecords: (n) => progress('prep', `Reading mate ${mate} (${file.name}): ${fmt(n)} records`) });
      checkFastq(r.bytes, label);
      return r;
    };
    const prepare = async () => {
      const r1 = await read(mate1, 1);
      const r2 = mate2 ? await read(mate2, 2) : null;
      if (r2 && r1.records !== r2.records) {
        const count = (r) => `${fmt(r.records)} records${r.wholeFile ? '' : ` in its first ${fmt(pairs)}`}`;
        throw new StepError('input', `The mates of sample ${name} have different record counts: mate 1 (${mate1.name}) has ${count(r1)}, ` +
          `mate 2 (${mate2.name}) has ${count(r2)}. kallisto pairs the records by position and stops at the shorter file without an error, so the page stops here. ` +
          'Check that the two files belong together and are complete.');
      }
      const files = [{ path: 'transcripts.idx', data: genome.index }, { path: 'r1.fastq', data: r1.bytes }];
      if (r2) files.push({ path: 'r2.fastq', data: r2.bytes });
      return { files, nPairs: r1.records, wholeFile: r1.wholeFile && (!r2 || r2.wholeFile) };
    };
    let prepared = await prepare();
    nPairs = prepared.nPairs;
    wholeFile = prepared.wholeFile;
    for (let attempt = 1; attempt <= 2; attempt++) {
      if (attempt === 2) {
        progress('prep', `The kallisto output failed the check (${firstProblems.join('; ')}). Reading the files again for one more run.`);
        prepared = await prepare();
      }
      progress('quant', `Running kallisto quant on ${fmt(nPairs)} ${unit}${attempt === 2 ? ' (second run)' : ''}`);
      const files = prepared.files;
      prepared = null;
      result = await kallistoQuant(args, files, signal);
      const problems = outputProblems(result, nPairs, nTranscripts, unit);
      if (!problems.length) break;
      if (attempt === 2) {
        throw new StepError('output', `kallisto wrote a damaged output twice for sample ${name} (first run: ${firstProblems.join('; ')}; second run: ${problems.join('; ')}). ` +
          'The page does not use it. A damaged abundance.tsv was seen in Firefox in the measurements (docs/browser/page-v1-measurements.md). Try again, or use another browser or the command-line tool.');
      }
      firstProblems = problems;
    }
  }
  post({ type: 'quant-done', id, name, nPairs, wholeFile, previewN: pairs, readLength, limitSource: limitOverride ? 'override' : 'read length', args, reruns: firstProblems ? 1 : 0, warning,
    abundance: result.abundance, runInfo: result.info });
}

/**
 * All reads of a sample. kallisto gets the picked File objects as files that it reads on demand (runtime execute.mjs,
 * mountBlobFile): it gunzips them itself, and only a read-ahead window of each file is in memory at a time. The reads are
 * never copied into the page, the worker message, or MEMFS, so memory does not grow with the read count (the measurements
 * found the kallisto heap flat in the read count; the growth was the copies). At the same time this worker reads each file
 * once with countRecords (nothing kept) to check the gzip and FASTQ format and to count the records: the mates must have the
 * same count, and kallisto must report that count as n_processed. A file that fails the check stops kallisto.
 */
async function quantAllReads({ id, name, fragment, nTranscripts, mate1, mate2, readLength, warning, unit }, signal, progress) {
  await ready;
  const ext = (f) => (/\.gz$/i.test(f.name) ? '.fastq.gz' : '.fastq');
  const r1 = `r1${ext(mate1)}`;
  const r2 = mate2 ? `r2${ext(mate2)}` : null;
  const args = JSON.parse(kallisto_args(JSON.stringify({
    kind: 'quant', sample: { name, fastq1: [r1], fastq2: mate2 ? [r2] : [] },
    index: 'transcripts.idx', out: 'kq', threads: 1, ...(mate2 ? {} : { fragment_len: fragment }),
  })));
  const local = new AbortController();
  const stop = () => local.abort(signal.reason);
  signal.addEventListener('abort', stop, { once: true });
  let result;
  let records = 0;
  let firstProblems = null;
  try {
    for (let attempt = 1; attempt <= 2; attempt++) {
      if (attempt === 2) progress('prep', `The kallisto output failed the check (${firstProblems.join('; ')}). Running kallisto again on the files.`);
      const started = performance.now();
      const counted = [0, 0];
      const report = () => {
        const s = Math.round((performance.now() - started) / 1000);
        progress('quant', `kallisto quant on all ${unit} of ${name}${attempt === 2 ? ' (second run)' : ''}: ` +
          `${fmt(counted[0])} ${unit} checked so far, ${Math.floor(s / 60)}m ${String(s % 60).padStart(2, '0')}s`);
      };
      report();
      const devices = [{ path: r1, blob: mate1 }];
      if (mate2) devices.push({ path: r2, blob: mate2 });
      const quant = kallistoQuant(args, [{ path: 'transcripts.idx', data: genome.index }], local.signal, devices);
      // The count is needed once: the files do not change between the two runs (a change is caught by n_processed).
      const counting = attempt === 1
        ? Promise.all([mate1, mate2].map((file, i) => file && countRecords(file, `Mate ${i + 1} of sample ${name} (${file.name})`, {
          signal: local.signal, onRecords: (n) => { counted[i] = n; report(); },
        })))
        : Promise.resolve(null);
      // Whichever fails first stops the other (a damaged file should not wait for kallisto, and the other way round).
      const [res, counts] = await Promise.all([
        quant.catch((e) => { local.abort(e); throw e; }),
        counting.catch((e) => { local.abort(e); throw e; }),
      ]);
      result = res;
      if (counts) {
        const [c1, c2] = counts;
        if (c2 && c1.records !== c2.records) {
          throw new StepError('input', `The mates of sample ${name} have different record counts: mate 1 (${mate1.name}) has ${fmt(c1.records)}, ` +
            `mate 2 (${mate2.name}) has ${fmt(c2.records)}. kallisto pairs the records by position and stops at the shorter file without an error, so the page stops here. ` +
            'Check that the two files belong together and are complete.');
        }
        records = c1.records;
      }
      const problems = outputProblems(result, records, nTranscripts, unit);
      if (!problems.length) break;
      if (attempt === 2) {
        throw new StepError('output', `kallisto wrote a damaged output twice for sample ${name} (first run: ${firstProblems.join('; ')}; second run: ${problems.join('; ')}). ` +
          'The page does not use it. Try again, or use another browser or the command-line tool.');
      }
      firstProblems = problems;
    }
  } catch (error) {
    if (signal.aborted) throw signal.reason ?? error;
    throw error;
  } finally {
    signal.removeEventListener('abort', stop);
  }
  post({ type: 'quant-done', id, name, nPairs: records, wholeFile: true, previewN: 0, readLength, limitSource: 'all reads', args,
    reruns: firstProblems ? 1 : 0, warning, abundance: result.abundance, runInfo: result.info });
}

/** Quant requests run one at a time, in the order they came. `quant-cancel` stops the one with that id. */
const quantRuns = new Map();
let quantQueue = Promise.resolve();

async function runQuantMessage(data, signal) {
  try {
    await quantSample(data, signal);
  } catch (error) {
    let kind = 'internal';
    if (error?.name === 'AbortError' || signal.aborted) kind = 'cancelled';
    else if (error instanceof StepError) kind = error.kind;
    else if (error instanceof ReadError) kind = 'input';
    else if (memoryLike(String(error?.message || error))) kind = 'memory';
    post({ type: 'error', id: data.id, kind, message: kind === 'cancelled' ? 'Cancelled.' : error?.message || String(error) });
  } finally {
    quantRuns.delete(data.id);
  }
}

self.onmessage = async ({ data }) => {
  if (data?.type === 'genome') {
    try {
      await buildGenome(data);
    } catch (error) {
      genome = null;
      const kind = error instanceof StepError ? error.kind : (memoryLike(String(error?.message || error)) ? 'memory' : 'internal');
      post({ type: 'error', id: data.id, kind, message: error?.message || String(error) });
    }
  } else if (data?.type === 'genome-prebuilt') {
    try {
      await loadPrebuiltGenome(data);
    } catch (error) {
      genome = null;
      const kind = error instanceof StepError ? error.kind : (memoryLike(String(error?.message || error)) ? 'memory' : 'internal');
      post({ type: 'error', id: data.id, kind, message: error?.message || String(error) });
    }
  } else if (data?.type === 'quant') {
    const controller = new AbortController();
    quantRuns.set(data.id, controller);
    quantQueue = quantQueue.then(() => runQuantMessage(data, controller.signal));
  } else if (data?.type === 'quant-cancel') {
    quantRuns.get(data.id)?.abort();
  } else if (data?.type === 'status') {
    post({ type: 'status', id: data.id, genome: genome && { label: genome.label, accession: genome.accession, indexBytes: genome.index.length } });
  } else if (data?.type === 'get-index') {
    post({ type: 'index-data', id: data.id, index: genome?.index || null });
  } else if (data?.type === 'transcripts') {
    // The transcript FASTA for the desktop app's native index: the kept one, or built from `fasta` and `gff` (a pre-built genome).
    try {
      if (genome?.transcripts) {
        post({ type: 'transcripts-ready', id: data.id, transcripts: genome.transcripts });
      } else if (data.fasta && data.gff) {
        const progress = (phase, text, done = 0, total = 0) => post({ type: 'progress', id: data.id, phase, text, done, total });
        const built = data.prepWorker === false ? await prepareTranscripts(data, progress) : await prepareInWorker(data, progress);
        self.postMessage({ type: 'transcripts-ready', id: data.id, transcripts: built.transcripts }, [built.transcripts.buffer]);
      } else {
        throw new StepError('no-index', 'The transcripts of this genome are not in the page worker. Choose the genome again.');
      }
    } catch (error) {
      post({ type: 'error', id: data.id, kind: error instanceof StepError ? error.kind : 'internal', message: error?.message || String(error) });
    }
  }
};
