/**
 * Browser spike 0, Task 4: measurement wrapper for one Emscripten tool module.
 * index_size.mjs gives this file to the upstream ToolRunner as moduleUrl, with the real tool module in the
 * query string: index_size_wrap.mjs?tool=<absolute URL of dist/<tool>/<tool>.mjs>.
 * The wrapper calls the real factory with the same options. It does not change the tool or its arguments.
 * It only wraps callMain to read the size of the WebAssembly linear memory (Module.HEAPU8.length) before and
 * after the run. Linear memory never shrinks, so the size after callMain is the peak linear memory of the run.
 * The values go to the caller in two ways:
 *   1. the file /work/spike0_heap.json in the tool file system (collected as a requested output on success);
 *   2. a BroadcastChannel message on channel "spike0-heap" (also sent when callMain throws, for example on OOM).
 */
const toolUrl = new URL(import.meta.url).searchParams.get('tool');
if (!toolUrl) throw new Error('index_size_wrap.mjs needs ?tool=<module URL>');
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
      try { module.FS.writeFile('/work/spike0_heap.json', JSON.stringify(record)); } catch {}
      try { const channel = new BroadcastChannel('spike0-heap'); channel.postMessage(record); channel.close(); } catch {}
    }
  };
  return module;
}
