// Quantification of the samples, one at a time, in the page worker (PROTOCOL.md, messages `quant` and `quant-cancel`).
// For each sample: the worker reads the preview reads, runs kallisto quant and checks its output; the page then sums the
// transcript counts per gene with `gene_counts` of qrnamap-wasm (the same code as the command-line tool) and the tx2gene
// table of the genome step. The page keeps only the gene counts and the mapping numbers of each sample.
// Cancel stops the sample in the worker (the reading stops and ToolRunner terminates the kallisto worker). The page worker
// and its index stay, so a new run needs no new index.
import { gene_counts } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';

/** The largest preview size: N = hard limit = 200,000 read pairs per sample (docs/browser/page-v1-measurements.md). */
export const MAX_PAIRS = 200000;

const fmt = (n) => Number(n).toLocaleString('en-US');

/** The label of every result made from a preview of `n` pairs. */
export const previewLabel = (n) => `preview of the first ${fmt(n)} pairs`;

/**
 * `app`: { state, getWorker, nextId } (window.qrnamapApp). `onChange(sample)` is called when a sample changes.
 * Returns { run(jobs, pairs), cancel(), busy() }. A job: { name, covariates, mate1 (File), mate2 (File or null), fragment ([mean, sd] or null) }.
 */
export function createQuant(app, { onChange = () => {}, onProgress = () => {} } = {}) {
  let current = null; // { id, worker, resolve, reject, listener }
  let runNumber = 0;

  const stopListening = () => {
    if (current) current.worker.removeEventListener('message', current.listener);
  };

  function request(message) {
    return new Promise((resolve, reject) => {
      const worker = app.getWorker();
      const id = app.nextId();
      const listener = ({ data }) => {
        if (data.id !== id) return;
        if (data.type === 'progress') onProgress(data);
        else if (data.type === 'quant-done') { stopListening(); current = null; resolve(data); }
        else if (data.type === 'error') { stopListening(); current = null; reject(Object.assign(new Error(data.message), { kind: data.kind })); }
      };
      worker.addEventListener('message', listener);
      current = { id, worker, resolve, reject, listener };
      worker.postMessage({ ...message, type: 'quant', id });
    });
  }

  /** The page worker was replaced (cancel or change of the genome, or a crash): its index is gone, and so is any open request. */
  window.addEventListener('qrnamap:worker-replaced', ({ detail }) => {
    if (!current) return;
    const { reject } = current;
    stopListening();
    current = null;
    reject(Object.assign(new Error(`The page worker was replaced (${detail?.reason || 'unknown reason'}), so the kallisto index is gone. Build the genome index again (step 1), then count the reads again.`), { kind: 'lost' }));
  });

  async function run(jobs, pairs) {
    const mine = ++runNumber;
    const state = app.state;
    state.preview = { n: pairs, label: previewLabel(pairs) };
    state.samples = jobs.map((j) => ({ name: j.name, covariates: { ...j.covariates }, status: 'waiting', nPairs: 0, processed: 0, pseudoaligned: 0,
      fraction: 0, allPairs: false, geneCounts: null, error: '', warning: '', reruns: 0 }));
    state.samplesPhase = 'running';
    state.samples.forEach(onChange);
    for (const [i, job] of jobs.entries()) {
      const sample = state.samples[i];
      sample.status = 'running';
      onChange(sample);
      try {
        const done = await request({ name: job.name, mate1: job.mate1, mate2: job.mate2, pairs, fragment: job.fragment, nTranscripts: state.nTranscripts });
        const counts = JSON.parse(gene_counts(done.abundance, state.tx2gene));
        Object.assign(sample, {
          status: 'done', nPairs: done.nPairs, processed: done.runInfo.n_processed, pseudoaligned: done.runInfo.n_pseudoaligned,
          fraction: done.nPairs ? done.runInfo.n_pseudoaligned / done.nPairs : 0, allPairs: done.wholeFile, geneCounts: counts.counts,
          reruns: done.reruns, args: done.args, warning: done.warning || '',
        });
      } catch (error) {
        if (mine !== runNumber) return; // a newer run took over
        const message = typeof error === 'string' ? `The gene counts could not be made: ${error}` : error.message; // gene_counts throws a string
        Object.assign(sample, { status: error.kind === 'cancelled' ? 'cancelled' : 'error', error: message, kind: error.kind || 'internal' });
        if (error.kind === 'lost' || error.kind === 'cancelled') {
          state.samples.slice(i + 1).forEach((s) => { s.status = 'not run'; onChange(s); });
          state.samplesPhase = error.kind === 'lost' ? 'lost' : 'cancelled';
          onChange(sample);
          return;
        }
      } finally {
        job.mate1 = null; // the page keeps no reference to the files of a finished sample
        job.mate2 = null;
      }
      onChange(sample);
      if (mine !== runNumber) return;
    }
    state.samplesPhase = 'done';
    document.dispatchEvent(new CustomEvent('qrnamap:samples-done', { detail: { samples: state.samples, preview: state.preview } }));
  }

  /** Stop the running sample. The worker stops reading or terminates its kallisto worker and answers `error` kind `cancelled`. */
  function cancel() {
    if (!current) return;
    current.worker.postMessage({ type: 'quant-cancel', id: current.id });
  }

  return { run, cancel, busy: () => app.state.samplesPhase === 'running' };
}
