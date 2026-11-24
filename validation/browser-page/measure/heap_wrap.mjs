/**
 * Page-v1 measurements: wrapper for the vendored kallisto module (site/vendor/kallisto/kallisto.mjs).
 * Adapted from validation/browser-spike0/index_size_wrap.mjs (same method, other channel and file name).
 * The harness gives this file to the vendored ToolRunner as moduleUrl: heap_wrap.mjs?tool=<absolute URL of kallisto.mjs>.
 * It calls the real factory with the same options and wraps callMain to read Module.HEAPU8.length before and after the run.
 * WebAssembly linear memory does not shrink, so the value after callMain is the peak linear memory of the run.
 * The values go to the caller as /work/measure_heap.json (a requested output) and on BroadcastChannel "measure-heap".
 */
const toolUrl = new URL(import.meta.url).searchParams.get('tool');
if (!toolUrl) throw new Error('heap_wrap.mjs needs ?tool=<module URL>');
const realFactory = (await import(toolUrl)).default;

export default async function factory(options) {
  const module = await realFactory(options);
  const initialBytes = module.HEAPU8?.length ?? null;
  const call = module.callMain;
  module.callMain = args => {
    const started = performance.now();
    let threw = null;
    try {
      return call(args);
    } catch (error) {
      threw = error?.name === 'ExitStatus' ? null : String(error?.message || error).slice(0, 400);
      throw error;
    } finally {
      const record = { initialBytes, finalBytes: module.HEAPU8?.length ?? null, callMainMs: performance.now() - started, threw };
      try { module.FS.writeFile('/work/measure_heap.json', JSON.stringify(record)); } catch {}
      try { const channel = new BroadcastChannel('measure-heap'); channel.postMessage(record); channel.close(); } catch {}
    }
  };
  return module;
}
