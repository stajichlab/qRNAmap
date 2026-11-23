// A short-lived module worker: reads the genome files and builds the transcripts, then the page worker ends it (page option
// `?prep_worker=1`; see prep.js). It posts { type: 'prepared', built } with the transcript FASTA transferred, or { type: 'error' }.
import { prepareTranscripts, StepError } from './prep.js?v=__BUILD__';

self.onmessage = async ({ data }) => {
  const progress = (phase, text, done = 0, total = 0) => self.postMessage({ type: 'progress', phase, text, done, total });
  try {
    const built = await prepareTranscripts(data, progress);
    self.postMessage({ type: 'prepared', built }, [built.transcripts.buffer]);
  } catch (error) {
    self.postMessage({ type: 'error', kind: error instanceof StepError ? error.kind : 'internal', message: error?.message || String(error) });
  }
};
