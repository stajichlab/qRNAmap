// Step 1 of /analyze/: choose a catalog genome (or upload a FASTA and a GFF3), download it, and give it to the page worker,
// which builds the transcripts and the kallisto index. The page keeps only the state below; the index stays in the worker.
// The worker protocol is in PROTOCOL.md. Cancel terminates the worker and starts a new one, so no work goes on after a cancel.
import init, { catalog } from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { checkUrl, mb } from '../reference-panel.js?v=__BUILD__';
import { readStore, writeStore } from '../persist.js?v=__BUILD__';
import { isTauri, prepareNativeIndex } from './tauri-bridge.js?v=__BUILD__';

/**
 * The largest genome that finished download, transcript building and kallisto index in Chromium, Firefox and WebKit in the
 * Task 3 measurements (docs/browser/page-v1-measurements.md, section "Genomes": Fulvia fulva Race5_Kim, 67,169,167 bases,
 * index peak up to 1,166.8 MiB in Firefox and 1,728.9 MiB in WebKit). No larger genome was tested, so the page refuses a larger one.
 * The count is the sequence letters of the FASTA (the catalog genome_size for a catalog genome). `?genome_limit_bases=N` can
 * lower the limit (the tests use it); it cannot raise it.
 */
export const LARGEST_TESTED = { label: 'Fulvia fulva Race5_Kim', accession: 'GCF_020509005.1', bases: 67169167 };
const askedLimit = Math.floor(Number(new URLSearchParams(location.search).get('genome_limit_bases')));
export const LIMIT_BASES = Number.isFinite(askedLimit) && askedLimit >= 1 ? Math.min(askedLimit, LARGEST_TESTED.bases) : LARGEST_TESTED.bases;

const fmt = (n) => Number(n).toLocaleString('en-US');
const STORE_KEY = 'last-genome';

const RELEASE_INDEX_URL = 'https://github.com/stajichlab/qRNAmap/releases/download/indexes-v1';

let prebuiltAccessions = new Set();
try {
  const prebuiltRes = await fetch('./prebuilt-indexes.json?v=__BUILD__');
  if (prebuiltRes.ok) {
    const list = await prebuiltRes.json();
    prebuiltAccessions = new Set(list);
  }
} catch (_) {}

/** What the command-line tool runs instead of this step. */
export function cliRoute(entry) {
  const fetchCmd = entry ? `qrnamap fetch --source ncbi --id ${entry.accession} --out genome/` : null;
  const ref = entry ? `genome/${entry.fasta.name}` : '<your FASTA>';
  const gff = entry ? `genome/${entry.gff.name}` : '<your GFF3>';
  return `Use the command-line tool instead (README, "Install and build"): ${fetchCmd ? `${fetchCmd}, then ` : ''}` +
    `qrnamap transcripts --reference ${ref} --gff ${gff}, then qrnamap quant.`;
}

function limitText(entry) {
  return `The largest genome tested in a browser is ${LARGEST_TESTED.label} (${LARGEST_TESTED.accession}, ${fmt(LARGEST_TESTED.bases)} bases; ` +
    `docs/browser/page-v1-measurements.md, section "Genomes")${LIMIT_BASES < LARGEST_TESTED.bases ? `, and this page now uses a lower limit of ${fmt(LIMIT_BASES)} bases` : ''}. ` +
    cliRoute(entry);
}

/** The state that the later steps read. `indexReady` is true only while the worker holds the index of `genome`. */
export const state = { phase: 'choose', genome: null, nTranscripts: 0, nGenes: 0, indexReady: false, tx2gene: '', genesTsv: '', geneInfoTsv: '', indexArgs: null, downloadNotes: [], nativeIndexPath: null };

let worker = null;
let pending = null; // { id, resolve, reject, onProgress }
let messageId = 0;

/** A message id for the page worker. Every step takes its ids from here, so no two messages share an id. */
export const nextId = () => ++messageId;

export function getIndexBytes() {
  return new Promise((resolve) => {
    const id = nextId();
    const w = getWorker();
    if (!w) return resolve(null);
    const handler = ({ data }) => {
      if (data.id === id && data.type === 'index-data') {
        w.removeEventListener('message', handler);
        resolve(data.index);
      }
    };
    w.addEventListener('message', handler);
    w.postMessage({ type: 'get-index', id });
  });
}

function startWorker() {
  worker = new Worker(new URL(`worker.js?v=__BUILD__`, import.meta.url), { type: 'module' });
  worker.onmessage = ({ data }) => {
    if (!pending || data.id !== pending.id) return; // an answer to a run that was cancelled
    if (data.type === 'progress') pending.onProgress(data);
    else if (data.type === 'genome-ready') { const p = pending; pending = null; p.resolve(data); }
    else if (data.type === 'error') { const p = pending; pending = null; p.reject(Object.assign(new Error(data.message), { kind: data.kind })); }
  };
  const crashed = (text) => {
    const message = `The page worker stopped (${text || 'no message from the browser'}). This often means the browser ran out of memory.`;
    if (pending) { const p = pending; pending = null; p.reject(Object.assign(new Error(message), { kind: 'memory' })); }
    else if (state.indexReady) hooks.lost?.(message);
    restartWorker('crash');
  };
  worker.onerror = (event) => { event.preventDefault?.(); crashed(event.message); };
  worker.onmessageerror = () => crashed('a message could not be read');
}

/**
 * Terminate the page worker and start a new one. Its genome state (the index) is lost. `reason` is 'cancel', 'crash' or
 * 'change'. The event `qrnamap:worker-replaced` on window tells the other steps (detail: { reason }).
 */
function restartWorker(reason) {
  state.nativeIndexPath = null;
  worker?.terminate();
  startWorker();
  window.dispatchEvent(new CustomEvent('qrnamap:worker-replaced', { detail: { reason } }));
}

/** The page worker. Task 6 sends its own message types to the same worker (PROTOCOL.md). */
export const getWorker = () => worker;

function runInWorker(message, onProgress) {
  return new Promise((resolve, reject) => {
    const id = nextId();
    pending = { id, resolve, reject, onProgress };
    worker.postMessage({ ...message, id });
  });
}

const hooks = {};

/**
 * Set up step 1. `$` finds an element by id. `shell` has markDone(name) and resetFrom(name) (app.js).
 * Returns { state, getWorker, cancel }.
 */
export async function createGenomeStep({ $, shell }) {
  startWorker();
  const section = $('step-genome');
  let entries = [];
  let run = 0; // the number of the current run; a run that sees another number stops
  let controller = null;

  const setPhase = (phase) => { state.phase = phase; section.dataset.phase = phase; };
  const busy = (flag) => {
    $('genome-choose').hidden = flag;
    $('genome-busy').hidden = !flag;
  };
  const status = (text, done = 0, total = 0) => {
    $('genome-status').textContent = text;
    const bar = $('genome-progress');
    if (total > 0) { bar.max = total; bar.value = Math.min(done, total); } else { bar.removeAttribute('value'); }
  };
  const showError = (text) => { const e = $('genome-error'); e.textContent = text; e.hidden = !text; };
  const note = (text) => { $('genome-note').textContent = text; };
  const changed = () => document.dispatchEvent(new CustomEvent('qrnamap:genome', { detail: { ...state } }));

  function clearGenome() {
    Object.assign(state, { genome: null, nTranscripts: 0, nGenes: 0, indexReady: false, tx2gene: '', genesTsv: '', geneInfoTsv: '', indexArgs: null, downloadNotes: [] });
    $('genome-ready').hidden = true;
    shell.resetFrom('genome');
  }

  function backToChoose() {
    setPhase('choose');
    busy(false);
    $('genome-ready').hidden = true;
    changed();
  }

  hooks.lost = (message) => {
    clearGenome();
    backToChoose();
    showError(`${message} The index is lost. Choose the genome again. ${cliRoute(null)}`);
  };

  // ------------------------------------------------------------------ the menu
  try {
    await init({ module_or_path: new URL('../pkg-web/qrnamap_wasm_bg.wasm?v=__BUILD__', import.meta.url) });
    entries = JSON.parse(catalog());
  } catch (error) {
    showError(`The genome catalog could not be loaded (${error?.message || error}). Upload a FASTA and a GFF3 instead.`);
  }
  // NCBI RefSeq first, then GenBank; by label in each group.
  entries.sort((a, b) => (a.source === 'RefSeq' ? 0 : 1) - (b.source === 'RefSeq' ? 0 : 1) || a.label.localeCompare(b.label));
  $('genome-count').textContent = `${entries.length} genomes with gene models (${entries.filter((e) => e.source === 'RefSeq').length} NCBI RefSeq, ` +
    `${entries.filter((e) => e.source !== 'RefSeq').length} GenBank). The page downloads the files from NCBI.`;

  function fillMenu() {
    const filter = $('genome-filter').value.trim().toLowerCase();
    const menu = $('genome-menu');
    const selected = menu.value;
    menu.replaceChildren();
    for (const source of ['RefSeq', 'GenBank']) {
      const group = document.createElement('optgroup');
      group.label = `NCBI ${source}`;
      for (const e of entries) {
        if ((e.source === 'RefSeq') !== (source === 'RefSeq')) continue;
        if (filter && !`${e.label} ${e.organism} ${e.accession}`.toLowerCase().includes(filter)) continue;
        group.appendChild(new Option(`${e.label} | ${e.accession} | ${(e.genome_size / 1e6).toFixed(1)} Mb | ${fmt(e.genes)} genes`, e.accession));
      }
      if (group.children.length) menu.appendChild(group);
    }
    if (selected && [...menu.options].some((o) => o.value === selected)) menu.value = selected;
    describe();
  }

  function describe() {
    const e = entries.find((x) => x.accession === $('genome-menu').value);
    if (!e) { $('genome-info').textContent = 'Choose a genome in the list.'; return; }
    const download = e.fasta.bytes + (e.gff?.bytes || 0);
    let text = `${e.label} (${e.accession}, ${e.source}, ${e.level}): ${fmt(e.genome_size)} bases, ${fmt(e.genes)} genes. ` +
      `Download about ${mb(download)} from NCBI.`;
    if (e.genome_size > LIMIT_BASES) text += ` This genome is larger than the page limit (${fmt(LIMIT_BASES)} bases).`;
    else if (e.genome_size > 50e6) text += ' In the measurements, the index of the three largest catalog genomes (63 to 67 Mb) took 48 to 120 s.';
    $('genome-info').textContent = text;
  }
  $('genome-filter').addEventListener('input', fillMenu);
  $('genome-menu').addEventListener('change', describe);
  fillMenu();

  for (const radio of document.querySelectorAll('input[name="genome-source"]')) {
    radio.addEventListener('change', () => {
      const upload = document.querySelector('input[name="genome-source"]:checked').value === 'upload';
      $('genome-catalog').hidden = upload;
      $('genome-upload').hidden = !upload;
    });
  }

  // The last genome (persist.js): select it in the menu again. The page does not download it again by itself.
  try {
    const raw = (await readStore([STORE_KEY]))[STORE_KEY];
    const text = (v) => typeof v === 'string' && v.length > 0 && v.length < 1000;
    // A stored value can be damaged or come from another version of the page: use it only when every field is a string.
    const saved = raw && typeof raw === 'object' && ((raw.source === 'catalog' && text(raw.accession) && text(raw.label))
      || (raw.source === 'upload' && text(raw.fasta) && text(raw.gff))) ? raw : null;
    if (saved?.source === 'catalog' && entries.some((e) => e.accession === saved.accession)) {
      $('genome-menu').value = saved.accession;
      describe();
      note(`Last genome: ${saved.label} (${saved.accession}). Press the button to build its index again.`);
    } else if (saved?.source === 'upload') {
      note(`Last genome: uploaded files ${saved.fasta} and ${saved.gff}. Choose the files again to use them.`);
    }
  } catch (error) { /* no browser store: the page works without it */ }

  // ------------------------------------------------------------------ download
  async function fetchFile(file, label, signal, onProgress) {
    const url = checkUrl(file.url, label);
    let response;
    try {
      response = await fetch(url, { signal });
    } catch (error) {
      if (signal.aborted) throw error;
      throw new Error(`${label}: the download from ${url.host} failed (${error.message}). Check the network connection and try again.`);
    }
    if (!response.ok) {
      const why = response.status === 404 ? ' (file not found)' : response.status >= 500 ? ' (a server error; try again later)' : '';
      throw new Error(`${label}: the server ${url.host} answered HTTP ${response.status}${why}.`);
    }
    const total = Number(response.headers.get('Content-Length')) || file.bytes || 0;
    const reader = response.body.getReader();
    const parts = [];
    let got = 0;
    try {
      for (;;) {
        const { done, value } = await reader.read();
        if (done) break;
        parts.push(value);
        got += value.length;
        onProgress(got, total);
      }
    } catch (error) {
      if (signal.aborted) throw error;
      throw new Error(`${label}: the download from ${url.host} stopped after ${mb(got)} (${error.message}). Check the network connection and try again.`);
    }
    const head = parts.length ? parts[0] : new Uint8Array(0);
    if (!(head.length >= 2 && head[0] === 0x1f && head[1] === 0x8b)) {
      const start = new TextDecoder().decode(head.subarray(0, 40)).replace(/[^\x20-\x7e]/g, '?');
      throw new Error(`${label}: the file from ${url.host} is not a gzip file from NCBI (${got} bytes; it starts with "${start}"). The server may have sent an error page. Try again later.`);
    }
    // A soft check: the catalog size was read on the catalog date. A different size is noted in the state, not refused.
    if (file.bytes && got !== file.bytes) state.downloadNotes.push(`${file.name}: ${got} bytes, the catalog says ${file.bytes}`);
    return new File(parts, file.name);
  }

  // ------------------------------------------------------------------ one run
  async function build({ entry, fasta, gff, label, accession }) {
    const mine = ++run;
    controller = new AbortController();
    const { signal } = controller;
    clearGenome();
    showError('');
    note('');
    busy(true);
    changed();
    try {
      let isPrebuilt = entry && accession && prebuiltAccessions.has(accession);
      let prebuiltPayload = null;

      if (isPrebuilt) {
        setPhase('download');
        try {
          status(`Downloading pre-built index for ${label}...`);
          const fetchAsset = async (suffix) => {
            const url = `${RELEASE_INDEX_URL}/${accession}.${suffix}`;
            const res = await fetch(url, { signal });
            if (!res.ok) throw new Error(`HTTP ${res.status}`);
            return new Uint8Array(await res.arrayBuffer());
          };

          const [idxGz, tx2geneGz, geneInfoGz, genesGz] = await Promise.all([
            fetchAsset('idx.gz'),
            fetchAsset('tx2gene.tsv.gz'),
            fetchAsset('gene_info.tsv.gz'),
            fetchAsset('genes.tsv.gz'),
          ]);

          prebuiltPayload = {
            type: 'genome-prebuilt',
            label,
            accession,
            bases: entry.genome_size,
            nTranscripts: entry.protein_coding_genes || 0,
            nGenes: entry.genes || 0,
            idxGz,
            tx2geneGz,
            geneInfoGz,
            genesGz,
          };
        } catch (e) {
          // If GitHub release asset download fails or is rate-limited, fall back transparently to in-browser build
          console.warn('Pre-built index download unavailable, falling back to in-browser indexing:', e);
          isPrebuilt = false;
        }
      }

      if (entry && !isPrebuilt) {
        if (entry.genome_size > LIMIT_BASES) {
          throw Object.assign(new Error(`${entry.label} (${entry.accession}) has ${fmt(entry.genome_size)} bases, more than the page limit of ${fmt(LIMIT_BASES)} bases. ${limitText(entry)}`), { kind: 'limit' });
        }
        setPhase('download');
        const files = [];
        for (const [slot, name] of [[entry.fasta, 'the genome FASTA'], [entry.gff, 'the GFF3']]) {
          files.push(await fetchFile(slot, name, signal, (got, total) => {
            if (mine === run) status(`Downloading ${name} from NCBI: ${mb(got)}${total ? ` of ${mb(total)}` : ''}`, got, total);
          }));
          if (mine !== run) return;
        }
        [fasta, gff] = files;
        files.length = 0;
      }
      if (mine !== run) return;
      setPhase('read');
      const saved = entry ? { source: 'catalog', accession, label } : { source: 'upload', fasta: fasta.name, gff: gff.name };
      const request = isPrebuilt
        ? prebuiltPayload
        : { type: 'genome', fasta, gff, label, accession, limitBases: LIMIT_BASES, limitText: limitText(entry), prepWorker: new URLSearchParams(location.search).get('prep_worker') !== '0' };
      fasta = null; // the page keeps no reference to the genome bytes while the worker runs
      gff = null;
      const answer = runInWorker(request, (p) => {
        if (mine !== run) return;
        setPhase(p.phase);
        status(p.text, p.done, p.total);
      });
      request.fasta = null;
      request.gff = null;
      const ready = await answer;
      if (mine !== run) return;
      Object.assign(state, {
        genome: { label, accession, source: entry ? 'catalog' : 'upload', bases: ready.bases },
        nTranscripts: ready.nTranscripts, nGenes: ready.nGenes, indexReady: true, tx2gene: ready.tx2gene, genesTsv: ready.genesTsv, geneInfoTsv: ready.geneInfoTsv, indexArgs: ready.indexArgs,
      });

      if (isTauri()) {
        try {
          if (prebuiltPayload?.idxGz) {
            const prep = await prepareNativeIndex({ index_gz_bytes: Array.from(prebuiltPayload.idxGz) });
            state.nativeIndexPath = prep.index_path;
          } else {
            const idxBytes = await getIndexBytes();
            if (idxBytes) {
              const prep = await prepareNativeIndex({ index_bytes: Array.from(idxBytes) });
              state.nativeIndexPath = prep.index_path;
            }
          }
        } catch (e) {
          console.warn('Native index preparation note:', e);
        }
      }

      setPhase('ready');
      busy(false);
      $('ready-label').textContent = label;
      $('ready-accession').textContent = accession || 'uploaded files';
      $('ready-genes').textContent = fmt(ready.nGenes);
      $('ready-transcripts').textContent = fmt(ready.nTranscripts);
      $('ready-time').textContent = `${ready.indexSeconds.toFixed(1)} s${isPrebuilt ? ' (pre-built)' : ''}`;
      $('genome-ready').hidden = false;
      shell.markDone('genome');
      changed();
      writeStore({ [STORE_KEY]: { ...saved, saved: Date.now() } }).then(() => { if (mine === run) section.dataset.saved = saved.source; }).catch(() => {});
    } catch (error) {
      if (mine !== run) return; // cancelled: cancel() has reset the step
      clearGenome();
      backToChoose();
      const route = ['memory', 'kallisto'].includes(error.kind) ? ` ${cliRoute(entry)}` : '';
      showError(`Error: ${error.message}${route}`);
    } finally {
      if (mine === run) controller = null;
    }
  }

  /** Stop the download or the worker at once and return to "choose genome", with nothing kept. */
  function cancel() {
    run++;
    controller?.abort();
    controller = null;
    pending = null;
    restartWorker('cancel');
    clearGenome();
    backToChoose();
    showError('');
    note('Cancelled. Choose a genome.');
  }

  $('genome-cancel').addEventListener('click', cancel);
  $('genome-load').addEventListener('click', () => {
    const entry = entries.find((e) => e.accession === $('genome-menu').value);
    if (!entry) { showError('Choose a genome in the list first.'); return; }
    build({ entry, label: entry.label, accession: entry.accession });
  });
  $('upload-load').addEventListener('click', () => {
    const fasta = $('upload-fasta').files[0];
    const gff = $('upload-gff').files[0];
    if (!fasta || !gff) { showError('Choose a genome FASTA file and a GFF3 file first.'); return; }
    build({ fasta, gff, label: `${fasta.name} + ${gff.name}`, accession: '' });
  });
  $('genome-change').addEventListener('click', () => {
    run++;
    restartWorker('change'); // drops the index of the old genome
    clearGenome();
    backToChoose();
  });

  async function loadCatalogAccession(accession) {
    const entry = entries.find((x) => x.accession === accession);
    if (!entry) throw new Error(`Genome accession ${accession} not in catalog`);
    const radio = document.querySelector('input[name="genome-source"][value="catalog"]');
    if (radio) radio.checked = true;
    $('genome-catalog').hidden = false;
    $('genome-upload').hidden = true;
    $('genome-filter').value = accession;
    fillMenu();
    $('genome-menu').value = accession;
    return build({ entry, label: entry.label, accession });
  }

  setPhase('choose');
  section.dataset.ready = 'menu'; // the tests wait for this
  return { state, getWorker, cancel, nextId, loadCatalogAccession, getIndexBytes };
}
