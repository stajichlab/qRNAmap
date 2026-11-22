// Tool shell: the three steps and the order in which they unlock. No analysis logic here yet.
// A later task fills the panels and calls `completeStep` when a step finishes.
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

/** Make `name` the active step again and lock every step after it (for example when the genome changes). */
export function resetFrom(name) {
  const i = ORDER.indexOf(name);
  ORDER.forEach((n, k) => {
    if (k < i) return;
    setState(n, k === i ? 'active' : 'locked');
    document.getElementById(`reason-${n}`).textContent = REASONS[n];
  });
}

window.qrnamapShell = { completeStep, resetFrom };
