// The page worker (a module worker). It holds the genome state of the page: the transcript FASTA, tx2gene, the gene table
// and the kallisto index. It never sends the index to the page. The message protocol is in site/analyze/PROTOCOL.md.
// kallisto runs through the vendored ToolRunner (vendor/kallisto/runtime/index.mjs), which starts one more worker per
// command and terminates it at the end, so a crash of kallisto does not end this worker.
import init, { build_transcriptome, kallisto_args } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { VerifyingGunzip } from '../gz.js?v=__BUILD__';
import { ToolRunner } from '../vendor/kallisto/runtime/index.mjs?v=__BUILD__';

const ready = init({ module_or_path: new URL('../pkg-web/qrnamap_wasm_bg.wasm?v=__BUILD__', import.meta.url) });
// Stamped like the page files, so a cached old kallisto cannot meet a new page. The runtime worker imports execute.mjs without a stamp
// (a vendored file, not changed).
const KALLISTO = new URL('../vendor/kallisto/kallisto.mjs?v=__BUILD__', import.meta.url);
const RUNTIME_WORKER = new URL('../vendor/kallisto/runtime/worker.mjs?v=__BUILD__', import.meta.url);
// A guard, not a measurement: a GFF3 larger than this after gunzip is refused, so a bad file cannot fill the memory of the tab.
// The largest GFF3 in the Task 3 genome runs is 64,659,993 bytes after gunzip (Laccaria bicolor, GCF_000143565.1;
// validation/browser-page/measure/results/genomes_runs.tsv, column gff_bytes). The guard is about 4 times that.
const MAX_GFF_BYTES = 260e6;

/** The genome of the page, or null. */
let genome = null;

class StepError extends Error {
  constructor(kind, message) { super(message); this.kind = kind; }
}

const post = (message) => self.postMessage(message);

/**
 * Read a File or Blob whole. A gzip file (first bytes 1f 8b) is unpacked through gz.js (fflate plus the CRC-32 and length
 * check of every member). `onChunk(plain)` sees every plain chunk (it can throw to stop). `progress(read, total)` gets the
 * compressed bytes read so far.
 */
async function readWhole(file, label, { onChunk = () => {}, progress = () => {}, maxBytes = Infinity } = {}) {
  if (file.size === 0) throw new StepError('input', `${label} is empty (0 bytes).`);
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  const gz = head.length === 2 && head[0] === 0x1f && head[1] === 0x8b;
  const parts = [];
  let size = 0;
  const keep = (chunk) => {
    size += chunk.length;
    if (size > maxBytes) throw new StepError('limit', `${label} unpacks to more than ${Math.round(maxBytes / 1e6)} MB, which the page does not load.`);
    onChunk(chunk);
    parts.push(chunk);
  };
  const reader = file.stream().getReader();
  let read = 0;
  try {
    const gunzip = gz ? new VerifyingGunzip(keep) : null;
    for (;;) {
      const { done, value } = await reader.read();
      if (value) read += value.length;
      if (gunzip) {
        try {
          gunzip.push(value || new Uint8Array(0), done);
        } catch (error) {
          if (error instanceof StepError) throw error;
          throw new StepError('gzip', `${label} is not a complete gzip file (${error.message}). The file is damaged or cut short. Download it again.`);
        }
      } else if (value) {
        keep(value.slice());
      }
      progress(read, file.size);
      if (done) break;
    }
  } finally {
    reader.cancel().catch(() => {});
  }
  const out = new Uint8Array(size);
  let offset = 0;
  for (const part of parts) { out.set(part, offset); offset += part.length; }
  if (size === 0) throw new StepError('input', `${label} is empty after gunzip.`);
  return out;
}

/** A counter of sequence letters in FASTA chunks (header lines and white space are not counted). */
function baseCounter(limit, label, limitText) {
  let bases = 0;
  let header = false;
  let lineStart = true;
  return {
    get bases() { return bases; },
    feed(chunk) {
      for (let i = 0; i < chunk.length; i++) {
        const b = chunk[i];
        if (b === 10) { header = false; lineStart = true; continue; }
        if (lineStart && b === 62) header = true; // '>'
        lineStart = false;
        if (!header && b !== 13 && b !== 32 && b !== 9) bases++;
      }
      if (bases > limit) throw new StepError('limit', `${label} has more than ${limit.toLocaleString('en-US')} bases. ${limitText}`);
    },
  };
}

function memoryLike(text) {
  return /out of memory|memory access out of bounds|could not allocate|cannot enlarge memory|allocation failed|Array buffer allocation|RangeError|OOM/i.test(text);
}

/**
 * Read both files and build the transcripts. Everything this function reads (the File objects, the genome bytes, the GFF3
 * bytes and text) is unreachable when it returns: only the transcript FASTA and the small results leave it.
 */
async function prepareTranscripts(message, progress) {
  const { fasta, gff, limitBases, limitText } = message;
  message.fasta = null; // the request object must not keep the files while kallisto runs
  message.gff = null;
  const mbText = (n) => `${(n / 1e6).toFixed(1)} MB`;
  progress('read', `Reading the genome FASTA (${fasta.name || 'FASTA'})`);
  const counter = baseCounter(limitBases, `The genome FASTA (${fasta.name || 'FASTA'})`, limitText);
  const fastaBytes = await readWhole(fasta, `The genome FASTA (${fasta.name || 'FASTA'})`, {
    onChunk: (c) => counter.feed(c),
    progress: (r, t) => progress('read', `Reading the genome FASTA: ${mbText(r)} of ${mbText(t)}`, r, t),
  });
  progress('read', `Reading the GFF3 (${gff.name || 'GFF3'})`);
  const gffText = new TextDecoder().decode(await readWhole(gff, `The GFF3 (${gff.name || 'GFF3'})`, {
    maxBytes: MAX_GFF_BYTES,
    progress: (r, t) => progress('read', `Reading the GFF3: ${mbText(r)} of ${mbText(t)}`, r, t),
  }));

  progress('build', `Building the transcripts (${counter.bases.toLocaleString('en-US')} bases)`);
  await ready;
  let tr;
  try {
    tr = build_transcriptome(fastaBytes, gffText);
  } catch (error) {
    throw new StepError('input', `The transcripts could not be built: ${typeof error === 'string' ? error : error?.message || error}`);
  }
  try {
    return { transcripts: tr.fasta, tx2gene: tr.tx2gene, genesTsv: tr.genes_tsv, nTranscripts: tr.n_transcripts, nGenes: tr.n_genes, bases: counter.bases };
  } finally {
    tr.free();
  }
}

/** The kallisto index argument list, from qrnamap-core (the same list as the Task 3 measurements). */
const indexArgs = () => JSON.parse(kallisto_args(JSON.stringify({ kind: 'index', fasta: 'transcripts.fa', index: 'transcripts.idx' })));

async function buildGenome(message) {
  const { id, label, accession } = message;
  genome = null; // a new genome replaces the old one at once, so no step can use a stale index
  const progress = (phase, text, done = 0, total = 0) => post({ type: 'progress', id, phase, text, done, total });
  const built = await prepareTranscripts(message, progress);

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

self.onmessage = async ({ data }) => {
  if (data?.type === 'genome') {
    try {
      await buildGenome(data);
    } catch (error) {
      genome = null;
      const kind = error instanceof StepError ? error.kind : (memoryLike(String(error?.message || error)) ? 'memory' : 'internal');
      post({ type: 'error', id: data.id, kind, message: error?.message || String(error) });
    }
  } else if (data?.type === 'status') {
    post({ type: 'status', id: data.id, genome: genome && { label: genome.label, accession: genome.accession, indexBytes: genome.index.length } });
  }
};
