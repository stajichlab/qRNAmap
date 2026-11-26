/** Internal worker implementation. Public callers should use ToolRunner. */
export async function execute({ moduleUrl, waitForExit = false, args, files, devices = [], outputs, outputDirectories = [], stdin }) {
  const factory = (await import(moduleUrl)).default;
  if (typeof factory !== 'function') throw new TypeError('Tool module must export an Emscripten factory as its default export');
  const stdout = new ByteCollector();
  const stderr = new ByteCollector();
  let inputOffset = 0;
  let exitCode = 0;
  let prepared = false;
  let finishExit;
  const exitEvent = new Promise(resolve => { finishExit = resolve; });
  const module = await factory({
    noInitialRun: true,
    print: (...values) => stdout.writeText(values.join(' ') + '\n'),
    printErr: (...values) => stderr.writeText(values.join(' ') + '\n'),
    onExit: status => { exitCode = status; finishExit({ status }); },
    onAbort: reason => finishExit({ error: new WebAssembly.RuntimeError(`Tool aborted: ${reason}`) }),
    // Emscripten's Node default must not terminate the worker before files are read.
    quit: (status, exception) => { exitCode = status; throw exception; },
    preRun: [instance => {
      const { FS } = instance;
      if (!FS) throw new Error('Tool build must export FS');
      FS.init(() => inputOffset < stdin.length ? stdin[inputOffset++] : null, byte => stdout.write(byte), byte => stderr.write(byte));
      FS.mkdirTree('/work');
      FS.chdir('/work');
      for (const file of files) {
        FS.mkdirTree(file.path.slice(0, file.path.lastIndexOf('/')));
        FS.writeFile(file.path, file.data);
      }
      let devMajor = 70;
      for (const dev of devices) {
        FS.mkdirTree(dev.path.slice(0, dev.path.lastIndexOf('/')));
        if (dev.blob) { mountBlobFile(FS, dev.path, dev.blob); continue; }
        const ring = new StreamRingConsumer(dev.sab);
        const devNum = FS.makedev(devMajor++, 0);
        FS.registerDevice(devNum, {
          open(stream) { stream.seekable = false; },
          close() {},
          read(stream, buffer, offset, length, pos) {
            return ring.read(buffer, offset, length);
          }
        });
        FS.mkdev(dev.path, 0o666, devNum);
      }
      prepared = true;
    }],
  });
  if (!prepared) throw new Error('Tool build did not run its preRun hooks');
  if (typeof module.callMain !== 'function') throw new Error('Tool build must export callMain');
  try {
    const result = module.callMain(args);
    if (Number.isInteger(result)) exitCode = result;
  } catch (error) {
    if (error?.name === 'ExitStatus' && Number.isInteger(error.status)) exitCode = error.status;
    else throw error;
  }
  // PROXY_TO_PTHREAD returns from callMain after starting the application.
  // The main runtime worker must keep serving filesystem proxy requests until
  // the application actually exits, before output files can be collected.
  if (waitForExit) {
    const outcome = await exitEvent;
    if (outcome.error) throw outcome.error;
    exitCode = outcome.status;
  }
  const resultFiles = [];
  const missingOutputs = [];
  const collected = new Set();
  const collectFile = path => {
    if (collected.has(path)) return;
    resultFiles.push({ path, data: module.FS.readFile(path) });
    collected.add(path);
  };
  for (const path of outputs) {
    if (!module.FS.analyzePath(path).exists) missingOutputs.push(path);
    else collectFile(path);
  }
  const collectDirectory = path => {
    const { FS } = module;
    if (!FS.isDir(FS.lstat(path).mode)) throw new Error(`Output directory is not a directory: ${path}`);
    for (const name of FS.readdir(path).filter(name => name !== '.' && name !== '..').sort()) {
      const child = `${path}/${name}`;
      const mode = FS.lstat(child).mode;
      if (FS.isDir(mode)) collectDirectory(child);
      else if (FS.isFile(mode)) collectFile(child);
      else throw new Error(`Output directory contains a non-regular file: ${child}`);
    }
  };
  for (const path of outputDirectories) {
    if (!module.FS.analyzePath(path).exists) missingOutputs.push(path);
    else collectDirectory(path);
  }
  return { exitCode, stdout: stdout.finish(), stderr: stderr.finish(), files: resultFiles, missingOutputs };
}

class ByteCollector {
  constructor() { this.chunks = []; this.chunk = new Uint8Array(65536); this.used = 0; this.length = 0; }
  writeText(value) { for (const byte of new TextEncoder().encode(value)) this.write(byte); }
  write(byte) {
    // Emscripten uses null/undefined to flush. A byte of zero is real output.
    if (byte === null || byte === undefined) return;
    if (this.used === this.chunk.length) { this.chunks.push(this.chunk); this.chunk = new Uint8Array(65536); this.used = 0; }
    this.chunk[this.used++] = byte;
    this.length++;
  }
  finish() {
    const result = new Uint8Array(this.length);
    let offset = 0;
    for (const chunk of this.chunks) { result.set(chunk, offset); offset += chunk.length; }
    result.set(this.chunk.subarray(0, this.used), offset);
    return result;
  }
}

class StreamRingConsumer {
  constructor(sab) {
    this.sab = sab;
    this.capacity = sab.byteLength - 64;
    this.header = new Int32Array(sab, 0, 16);
    this.data = new Uint8Array(sab, 64, this.capacity);
  }
  read(dest, destOffset, maxLen) {
    while (true) {
      const status = Atomics.load(this.header, 0);
      if (status === 2) throw new Error('Stream ring buffer was aborted');
      const avail = Atomics.load(this.header, 1);
      if (avail === 0) {
        if (status === 1) return 0; // EOF
        const waitVal = Atomics.load(this.header, 4);
        Atomics.wait(this.header, 4, waitVal, 50);
        continue;
      }
      const rPtr = Atomics.load(this.header, 3);
      const toRead = Math.min(maxLen, avail, this.capacity - rPtr);
      dest.set(this.data.subarray(rPtr, rPtr + toRead), destOffset);
      Atomics.store(this.header, 3, (rPtr + toRead) % this.capacity);
      Atomics.sub(this.header, 1, toRead);
      Atomics.add(this.header, 5, 1);
      Atomics.notify(this.header, 5);
      return toRead;
    }
  }
}


/**
 * Make `path` a read-only regular file whose bytes are read from `blob` (a File or Blob) only when the tool reads them,
 * like Emscripten's WORKERFS but without a rebuild: a MEMFS file node with its size taken from the blob and its `read`
 * replaced. The bytes come from `FileReaderSync` (dedicated workers have it) through one read-ahead window of
 * BLOB_WINDOW bytes, so the file is never held whole: memory is one window, whatever the size of the file.
 * No SharedArrayBuffer and no cross-origin isolation are needed. `readSync(blob, start, end)` replaces FileReaderSync
 * (tests in Node, which has no FileReaderSync).
 */
export const BLOB_WINDOW = 4 * 1024 * 1024;
export function mountBlobFile(FS, path, blob, { readSync = defaultReadSync, window = BLOB_WINDOW } = {}) {
  const slash = path.lastIndexOf('/');
  const node = FS.createFile(path.slice(0, slash) || '/', path.slice(slash + 1), {}, true, false);
  const size = blob.size;
  // MEMFS takes the file size from usedBytes (stat, and lseek with SEEK_END).
  Object.defineProperty(node, 'usedBytes', { get: () => size, configurable: true });
  let winStart = 0;
  let win = new Uint8Array(0);
  const ops = { ...node.stream_ops };
  ops.read = (stream, buffer, offset, length, position) => {
    if (position >= size || length <= 0) return 0;
    const n = Math.min(length, size - position);
    if (position < winStart || position + n > winStart + win.length) {
      winStart = position;
      win = readSync(blob, position, Math.min(size, position + Math.max(window, n)));
      if (win.length < Math.min(window, size - position)) throw new Error(`Short read of ${path} at byte ${position}`);
    }
    buffer.set(win.subarray(position - winStart, position - winStart + n), offset);
    return n;
  };
  ops.write = () => { throw new FS.ErrnoError(1); }; // EPERM: the file is read-only
  ops.mmap = () => { throw new FS.ErrnoError(19); }; // ENODEV: not mappable
  node.stream_ops = ops;
  return node;
}

function defaultReadSync(blob, start, end) {
  if (typeof FileReaderSync !== 'function') throw new Error('Reading a File during the run needs FileReaderSync (a dedicated worker)');
  return new Uint8Array(new FileReaderSync().readAsArrayBuffer(blob.slice(start, end)));
}
