// The page worker (a module worker). After a genome is ready it holds the kallisto index (with the genome label and
// accession). It never sends the index to the page. It also reads the preview reads of one sample at a time and runs
// kallisto quant on them (quant.js sends the requests). The message protocol is in site/analyze/PROTOCOL.md.
// kallisto runs through the vendored ToolRunner (vendor/kallisto/runtime/index.mjs), which starts one more worker per
// command and terminates it at the end, so a crash of kallisto does not end this worker.
import init, { kallisto_args } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { prepareTranscripts, StepError } from './prep.js?v=__BUILD__';
import { ToolRunner } from '../vendor/kallisto/runtime/index.mjs?v=__BUILD__';
import { firstRecords, checkFastq, firstNameStem, firstReadLength, ReadError } from './fastq.js?v=__BUILD__';
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
 * Page option `?prep_worker=1` (a measurement option, off by default): read the genome files and build the transcripts in a
 * short-lived worker and end it, so the memory of the genome bytes is released before the kallisto index is built.
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
  // The worker keeps only the index. tx2gene and the gene table go to the page, which keeps them.
  genome = { label, accession, index };
  post({ type: 'genome-ready', id, label, accession, nTranscripts: built.nTranscripts, nGenes: built.nGenes, bases: built.bases,
    indexBytes: index.length, indexSeconds: (performance.now() - started) / 1000, indexArgs: args, tx2gene: built.tx2gene, genesTsv: built.genesTsv });
}

// ------------------------------------------------------------------ quantification of one sample

const fmt = (n) => Number(n).toLocaleString('en-US');

/**
 * One kallisto quant run. Gives { abundance (text), info (run_info.json numbers, or null when the file cannot be read),
 * infoError }. `files` is emptied as soon as ToolRunner has the run, so only ToolRunner refers to the reads. `signal` stops it.
 */
async function kallistoQuant(args, files, signal) {
  const runner = new ToolRunner({ moduleUrl: KALLISTO, workerUrl: RUNTIME_WORKER, waitForExit: false });
  let run;
  try {
    const running = runner.run(args, { files: [...files], outputs: ['kq/abundance.tsv', 'kq/run_info.json'], signal });
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
  for (const [file, mate] of [[mate1, 1], [mate2, 2]]) {
    if (!file) continue;
    const label = `Mate ${mate} of sample ${name} (${file.name})`;
    const first = await firstRecords(file, 1, label, { signal });
    checkFastq(first.bytes, label);
    readLength = Math.max(readLength, firstReadLength(first.bytes));
  }
  if (readLength < 1) throw new StepError('input', `The first record of sample ${name} has no bases.`);
  const limit = limitOverride ? Math.min(limitOverride, HARD_LIMIT) : previewLimit(readLength);
  const pairs = Math.min(message.pairs, limit);
  const read = async (file, mate) => {
    const label = `Mate ${mate} of sample ${name} (${file.name})`;
    progress('prep', `Reading the first ${fmt(pairs)} records of mate ${mate} (${file.name})`);
    const r = await firstRecords(file, pairs, label, { signal, onRecords: (n) => progress('prep', `Reading mate ${mate} (${file.name}): ${fmt(n)} records`) });
    checkFastq(r.bytes, label);
    return r;
  };
  /** The kallisto input files of this sample, read from the File objects. */
  const prepare = async () => {
    const r1 = await read(mate1, 1);
    const r2 = mate2 ? await read(mate2, 2) : null;
    if (r2 && r1.records !== r2.records) {
      const count = (r) => `${fmt(r.records)} records${r.wholeFile ? '' : ` in its first ${fmt(pairs)}`}`;
      throw new StepError('input', `The mates of sample ${name} have different record counts: mate 1 (${mate1.name}) has ${count(r1)}, ` +
        `mate 2 (${mate2.name}) has ${count(r2)}. kallisto pairs the records by position and stops at the shorter file without an error, so the page stops here. ` +
        'Check that the two files belong together and are complete.');
    }
    let warning = '';
    if (r2) {
      const a = firstNameStem(r1.bytes);
      const b = firstNameStem(r2.bytes);
      if (a !== b) warning = `the first read names of the two mates differ ("${a}" and "${b}"). Check that the two files belong together.`;
    }
    const files = [{ path: 'transcripts.idx', data: genome.index }, { path: 'r1.fastq', data: r1.bytes }];
    if (r2) files.push({ path: 'r2.fastq', data: r2.bytes });
    return { files, nPairs: r1.records, wholeFile: r1.wholeFile && (!r2 || r2.wholeFile), warning };
  };
  let prepared = await prepare();
  const { nPairs, wholeFile, warning } = prepared;
  const unit = mate2 ? 'pairs' : 'reads';
  await ready;
  const args = JSON.parse(kallisto_args(JSON.stringify({
    kind: 'quant', sample: { name, fastq1: ['r1.fastq'], fastq2: mate2 ? ['r2.fastq'] : [] },
    index: 'transcripts.idx', out: 'kq', threads: 1, ...(mate2 ? {} : { fragment_len: fragment }),
  })));
  let result;
  let firstProblems = null;
  for (let attempt = 1; attempt <= 2; attempt++) {
    if (attempt === 2) {
      progress('prep', `The kallisto output failed the check (${firstProblems.join('; ')}). Reading the files again for one more run.`);
      prepared = await prepare();
    }
    progress('quant', `Running kallisto quant on ${fmt(nPairs)} ${unit}${attempt === 2 ? ' (second run)' : ''}`);
    const files = prepared.files;
    prepared = null; // the reads are referenced only by `files`, which kallistoQuant empties
    result = await kallistoQuant(args, files, signal);
    const problems = outputProblems(result, nPairs, nTranscripts, unit);
    if (!problems.length) break;
    if (attempt === 2) {
      throw new StepError('output', `kallisto wrote a damaged output twice for sample ${name} (first run: ${firstProblems.join('; ')}; second run: ${problems.join('; ')}). ` +
        'The page does not use it. A damaged abundance.tsv was seen in Firefox in the measurements (docs/browser/page-v1-measurements.md). Try again, or use another browser or the command-line tool.');
    }
    firstProblems = problems;
  }
  post({ type: 'quant-done', id, name, nPairs, wholeFile, previewN: pairs, readLength, limitSource: limitOverride ? 'override' : 'read length', args, reruns: firstProblems ? 1 : 0, warning,
    abundance: result.abundance, runInfo: result.info });
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
  } else if (data?.type === 'quant') {
    const controller = new AbortController();
    quantRuns.set(data.id, controller);
    quantQueue = quantQueue.then(() => runQuantMessage(data, controller.signal));
  } else if (data?.type === 'quant-cancel') {
    quantRuns.get(data.id)?.abort();
  } else if (data?.type === 'status') {
    post({ type: 'status', id: data.id, genome: genome && { label: genome.label, accession: genome.accession, indexBytes: genome.index.length } });
  }
};
