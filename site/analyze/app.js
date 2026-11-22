// Tool shell: the three steps and the order in which they unlock. Step 1 is genome.js, step 2 is samples.js (with quant.js).
// Step 3 (results) is not built yet: a finished samples step is marked done with `markDone`, and step 3 stays locked.
import { createGenomeStep } from './genome.js?v=__BUILD__';
import { createSamplesStep } from './samples.js?v=__BUILD__';

const ORDER = ['genome', 'samples', 'results'];
const REASONS = {
  genome: 'Choose the genome first. The tool builds the transcript index from it.',
  samples: 'Locked. Finish step 1 (the genome and its index) to add samples.',
  results: 'Locked. Finish step 2 (at least one sample counted) to see results.',
};
const NOT_BUILT = { results: 'Locked. The samples are counted, but the results view is not built yet.' };

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

/** Mark `name` done and keep the next step locked, with a reason that says it is not built yet. */
export function markDone(name) {
  const next = ORDER[ORDER.indexOf(name) + 1];
  setState(name, 'done');
  if (next && NOT_BUILT[next]) document.getElementById(`reason-${next}`).textContent = NOT_BUILT[next];
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

const $ = (id) => document.getElementById(id);
const genomeStep = await createGenomeStep({ $, shell: { markDone: completeStep, resetFrom } });
// What the later steps read (PROTOCOL.md, "Page state"): genome, nTranscripts, nGenes, indexReady, tx2gene, genesTsv (step 1);
// samples, preview, samplesPhase (step 2). The index stays in the worker.
// nextId(): the id for any message to the page worker (one counter for all steps, so ids never collide).
// The event `qrnamap:worker-replaced` on window says that the worker was replaced and its index is gone.
// The event `qrnamap:samples-done` on document says that a run of all samples has ended.
window.qrnamapApp = { state: genomeStep.state, getWorker: genomeStep.getWorker, cancelGenome: genomeStep.cancel, nextId: genomeStep.nextId };
createSamplesStep({ $, app: window.qrnamapApp, shell: { markDone, resetFrom } });
document.getElementById('step-samples').dataset.ready = 'table'; // the tests wait for this
