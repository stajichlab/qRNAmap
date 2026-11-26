// Tool shell: the three steps and the order in which they unlock. Step 1 is genome.js, step 2 is samples.js (with quant.js), step 3 is results.js.
import { createGenomeStep } from './genome.js?v=__BUILD__';
import { createSamplesStep } from './samples.js?v=__BUILD__';
import { createResultsStep } from './results.js?v=__BUILD__';
import { createTestsets } from './testsets.js?v=__BUILD__';
import { clearStore } from '../persist.js?v=__BUILD__';
import * as wasm from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';
import { isTauri } from './tauri-bridge.js?v=__BUILD__';

// The desktop window is wide and has no phone layout: style.css widens the page under html.desktop.
if (isTauri()) document.documentElement.classList.add('desktop');

const ORDER = ['genome', 'samples', 'results'];
const REASONS = {
  genome: 'Choose the genome first. The tool builds the transcript index from it.',
  samples: 'Locked. Finish step 1 (the genome and its index) to add samples.',
  results: 'Locked. Finish step 2 (at least one sample counted) to see results.',
};

function setState(name, state) {
  const section = document.getElementById(`step-${name}`);
  section.dataset.state = state;
  document.getElementById(`fs-${name}`).disabled = state === 'locked';
  document.getElementById(`reason-${name}`).hidden = state !== 'locked' && name !== 'genome';
}

/** Mark `name` done and unlock the next step. */
export function completeStep(name) {
  const i = ORDER.indexOf(name);
  setState(name, 'done');
  if (i + 1 < ORDER.length) setState(ORDER[i + 1], 'active');
}

/** Mark `name` done and keep the next step locked. */
export function markDone(name) {
  completeStep(name);
}

/** Make `name` the active step again and lock every step after it (for example when the genome changes). */
export function resetFrom(name) {
  const i = ORDER.indexOf(name);
  ORDER.forEach((n, k) => {
    if (k < i) return;
    setState(n, k === i ? 'active' : 'locked');
    document.getElementById(`reason-${n}`).textContent = REASONS[n];
  });
}

window.qrnamapShell = { completeStep, markDone, resetFrom };
window.qrnamapWasm = wasm;

const $ = (id) => document.getElementById(id);
const genomeStep = await createGenomeStep({ $, shell: { markDone: completeStep, resetFrom } });
window.qrnamapApp = {
  state: genomeStep.state,
  getWorker: genomeStep.getWorker,
  cancelGenome: genomeStep.cancel,
  nextId: genomeStep.nextId,
  loadCatalogAccession: genomeStep.loadCatalogAccession,
  getIndexBytes: genomeStep.getIndexBytes,
};
const samplesStep = createSamplesStep({ $, app: window.qrnamapApp, shell: { markDone: completeStep, resetFrom } });
window.qrnamapApp.samplesStep = samplesStep;
const resultsStep = createResultsStep({ $, app: window.qrnamapApp, shell: { completeStep, resetFrom } });
window.qrnamapApp.resultsStep = resultsStep;
const testsets = createTestsets({ $, app: window.qrnamapApp, shell: { completeStep, resetFrom } });
window.qrnamapApp.testsets = testsets;

$('btn-reset-all')?.addEventListener('click', async () => {
  const confirmed = window.confirm(
    'Reset qRNAmap and clear all memory to start over?\n\nThis unloads any in-memory genome index, samples, and cached results.'
  );
  if (!confirmed) return;
  try {
    window.qrnamapApp?.cancelGenome?.();
  } catch (_) {}
  try {
    await clearStore();
  } catch (_) {}
  try {
    for (const key of Object.keys(localStorage)) {
      if (key !== 'qrnamap-theme' && key !== 'qrnamap-layout') {
        localStorage.removeItem(key);
      }
    }
  } catch (_) {}
  location.reload();
});

document.getElementById('step-samples').dataset.ready = 'table'; // the tests wait for this
