// Quantification of the samples, one at a time, in the page worker (PROTOCOL.md, messages `quant` and `quant-cancel`).
// For each sample: the worker reads the preview reads, runs kallisto quant and checks its output; the page then sums the
// transcript counts per gene with `gene_counts` of qrnamap-wasm (the same code as the command-line tool) and the tx2gene
// table of the genome step. The page keeps only the gene counts and the mapping numbers of each sample.
// Cancel stops the sample in the worker (the reading stops and ToolRunner terminates the kallisto worker). The page worker
// and its index stay, so a new run needs no new index.
import { gene_counts } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { HARD_LIMIT } from './preview-limit.js?v=__BUILD__';
import { isTauri, runNativeQuant, prepareNativeIndex } from './tauri-bridge.js?v=__BUILD__';

/** The largest value of the page field "Read pairs per sample". The N of each sample also depends on its read length (preview-limit.js). */
export const MAX_PAIRS = HARD_LIMIT;

const fmt = (n) => Number(n).toLocaleString('en-US');

/** The label of every result made from a preview of `n` pairs. */
export const previewLabel = (n) => `preview of the first ${fmt(n)} pairs`;

/**
 * The preview of the counted samples: `{ n, label }`. `n` is the N that every counted sample used, or null when the samples used
 * different N (then the label gives the range, and each sample has its own `previewN`). Empty label when no sample is counted.
 */
export function previewOf(samples) {
  const doneSamples = samples.filter((s) => s.status === 'done');
  if (!doneSamples.length) return { n: null, label: '' };
  const allFull = doneSamples.every((s) => s.allPairs || s.previewN === null || s.previewN === 0);
  if (allFull) return { n: null, label: 'full depth (all read pairs directly from disk)' };
  const used = [...new Set(doneSamples.map((s) => s.previewN).filter((n) => typeof n === 'number' && n > 0))].sort((a, b) => a - b);
  if (!used.length) return { n: null, label: 'full depth (all read pairs directly from disk)' };
  if (used.length === 1) return { n: used[0], label: previewLabel(used[0]) };
  return { n: null, label: `preview of the first ${fmt(used[0])} to ${fmt(used[used.length - 1])} pairs (the N of each sample is in its row)` };
}

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
        else if (data.type === 'quant-done' && current?.cancelRequested) {
          // The worker finished the sample before it read the cancel: the user asked to cancel, so the counts are not used.
          stopListening(); current = null; reject(Object.assign(new Error('Cancelled.'), { kind: 'cancelled' }));
        } else if (data.type === 'quant-done') { stopListening(); current = null; resolve(data); }
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

  /** `pairs`: the page field (an upper bound). `limitOverride`: the measurement option ?preview_pairs (replaces the read-length limit), or 0. */
  async function run(jobs, pairs, limitOverride = 0) {
    const mine = ++runNumber;
    const state = app.state;
    state.preview = { n: null, label: '' };
    state.samples = jobs.map((j) => ({
      name: j.name,
      covariates: { ...j.covariates },
      fastq1: j.mate1?.name || (j.paths1 ? j.paths1.map((p) => p.split(/[\\/]/).pop()).join(',') : '') || j.fastq1 || '',
      fastq2: j.mate2?.name || (j.paths2 ? j.paths2.map((p) => p.split(/[\\/]/).pop()).join(',') : '') || j.fastq2 || '',
      status: 'waiting', nPairs: 0, processed: 0, pseudoaligned: 0,
      fraction: 0, allPairs: false, previewN: 0, readLength: 0, geneCounts: null, error: '', warning: '', reruns: 0,
    }));
    state.samplesPhase = 'running';
    state.samples.forEach(onChange);

    const isDesktopNative = isTauri() &&
      (document.getElementById('mode-native')?.checked ?? true) &&
      jobs.some((j) => j.paths1 && j.paths1.length > 0);

    if (isDesktopNative) {
      if (!state.nativeIndexPath) {
        onProgress({ text: 'Preparing native kallisto index on disk...' });
        const idxBytes = await app.getIndexBytes?.();
        if (idxBytes) {
          const prep = await prepareNativeIndex({ index_bytes: Array.from(idxBytes) });
          state.nativeIndexPath = prep.index_path;
        } else {
          throw new Error('Native kallisto index is not available. Please build the genome index in step 1 first.');
        }
      }
    }

    current = { id: 0, cancelRequested: false };

    for (const [i, job] of jobs.entries()) {
      const sample = state.samples[i];
      sample.status = 'running';
      onChange(sample);
      if (current?.cancelRequested) {
        sample.status = 'cancelled';
        onChange(sample);
        state.samples.slice(i + 1).forEach((s) => { s.status = 'not run'; onChange(s); });
        state.samplesPhase = 'cancelled';
        return;
      }
      try {
        if (isDesktopNative && job.paths1 && job.paths1.length > 0) {
          onProgress({ text: `Running native kallisto on ${job.name} (full depth)...` });
          const resp = await runNativeQuant({
            name: job.name,
            fastq1: job.paths1,
            fastq2: job.paths2 || [],
            index_path: state.nativeIndexPath,
            tx2gene: state.tx2gene,
            strandedness: app.state.strandedness || 'unstranded',
            fragment_mean: job.fragment ? job.fragment[0] : null,
            fragment_sd: job.fragment ? job.fragment[1] : null,
            threads: null,
            out_dir: null,
          });
          Object.assign(sample, {
            status: 'done',
            nPairs: resp.n_processed,
            processed: resp.n_processed,
            pseudoaligned: resp.n_pseudoaligned,
            fraction: resp.p_pseudoaligned || (resp.n_processed ? resp.n_pseudoaligned / resp.n_processed : 0),
            allPairs: true,
            geneCounts: resp.gene_counts,
            reruns: 0,
            args: resp.args,
            warning: resp.warning || '',
            previewN: 0,
            readLength: 0,
          });
        } else {
          const done = await request({ name: job.name, mate1: job.mate1, mate2: job.mate2, pairs, limitOverride, fragment: job.fragment, nTranscripts: state.nTranscripts });
          const counts = JSON.parse(gene_counts(done.abundance, state.tx2gene));
          Object.assign(sample, {
            status: 'done', nPairs: done.nPairs, processed: done.runInfo.n_processed, pseudoaligned: done.runInfo.n_pseudoaligned,
            fraction: done.nPairs ? done.runInfo.n_pseudoaligned / done.nPairs : 0, allPairs: done.wholeFile, geneCounts: counts.counts,
            reruns: done.reruns, args: done.args, warning: done.warning || '', previewN: done.previewN, readLength: done.readLength,
          });
        }
      } catch (error) {
        if (mine !== runNumber) return; // a newer run took over
        const message = typeof error === 'string' ? error : error.message;
        Object.assign(sample, { status: error.kind === 'cancelled' || current?.cancelRequested ? 'cancelled' : 'error', error: message, kind: error.kind || 'internal' });
        if (error.kind === 'lost' || error.kind === 'cancelled' || current?.cancelRequested) {
          state.samples.slice(i + 1).forEach((s) => { s.status = 'not run'; onChange(s); });
          state.samplesPhase = error.kind === 'lost' ? 'lost' : 'cancelled';
          onChange(sample);
          return;
        }
      } finally {
        job.mate1 = null; // the page keeps no reference to the files of a finished sample
        job.mate2 = null;
      }
      state.preview = previewOf(state.samples);
      onChange(sample);
      if (mine !== runNumber) return;
    }
    current = null;
    state.samplesPhase = 'done';
    document.dispatchEvent(new CustomEvent('qrnamap:samples-done', { detail: { samples: state.samples, preview: state.preview } }));
  }

  /** Stop the running sample. The worker stops reading or terminates its kallisto worker and answers `error` kind `cancelled`. */
  function cancel() {
    if (!current) return;
    current.cancelRequested = true;
    if (current.worker && current.id) {
      current.worker.postMessage({ type: 'quant-cancel', id: current.id });
    }
  }

  return { run, cancel, busy: () => app.state.samplesPhase === 'running' };
}
