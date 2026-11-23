// Tool shell: the three steps and the order in which they unlock. Step 1 is genome.js, step 2 is samples.js (with quant.js), step 3 is results.js.
import { createGenomeStep } from './genome.js?v=__BUILD__';
import { createSamplesStep } from './samples.js?v=__BUILD__';
import { createResultsStep } from './results.js?v=__BUILD__';
import * as wasm from '../pkg-web/qrnamap_wasm.js?v=__BUILD__';

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
window.qrnamapApp = { state: genomeStep.state, getWorker: genomeStep.getWorker, cancelGenome: genomeStep.cancel, nextId: genomeStep.nextId };
createSamplesStep({ $, app: window.qrnamapApp, shell: { markDone: completeStep, resetFrom } });
createResultsStep({ $, app: window.qrnamapApp, shell: { completeStep, resetFrom } });
document.getElementById('step-samples').dataset.ready = 'table'; // the tests wait for this
