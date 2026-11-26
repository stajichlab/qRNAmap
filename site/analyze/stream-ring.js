// A bounded SharedArrayBuffer ring buffer connecting a streaming reader (producer)
// to kallisto.wasm via an Emscripten block device (consumer) without unbounded memory growth.

const DEFAULT_CAPACITY = 2 * 1024 * 1024; // 2 MB buffer per mate
const HEADER_SIZE = 64; // 16 x 4-byte Int32 slots

export const STATUS_ACTIVE = 0;
export const STATUS_CLOSED = 1;
export const STATUS_ABORT = 2;

// Header offsets (Int32Array):
// [0] status: 0=active, 1=closed (EOF), 2=aborted
// [1] available: unread bytes in buffer (0..capacity)
// [2] writePtr: circular write index (0..capacity-1)
// [3] readPtr: circular read index (0..capacity-1)
// [4] readerWait: wake counter for reader
// [5] writerWait: wake counter for writer

export class StreamRingBuffer {
  static DEFAULT_CAPACITY = DEFAULT_CAPACITY;
  static HEADER_SIZE = HEADER_SIZE;

  constructor(sab = null, capacity = DEFAULT_CAPACITY) {
    this.capacity = capacity;
    if (!sab) {
      if (typeof SharedArrayBuffer === 'undefined') {
        throw new Error('SharedArrayBuffer is not supported in this environment');
      }
      this.sab = new SharedArrayBuffer(HEADER_SIZE + capacity);
    } else {
      this.sab = sab;
      this.capacity = sab.byteLength - HEADER_SIZE;
    }
    this.header = new Int32Array(this.sab, 0, 16);
    this.data = new Uint8Array(this.sab, HEADER_SIZE, this.capacity);
  }

  async waitForSpace() {
    while (true) {
      if (Atomics.load(this.header, 0) === STATUS_ABORT) {
        throw new Error('Stream ring buffer was aborted');
      }
      const avail = Atomics.load(this.header, 1);
      const free = this.capacity - avail;
      if (free > 0) return free;
      const waitVal = Atomics.load(this.header, 5);
      const waitRes = typeof Atomics.waitAsync === 'function'
        ? Atomics.waitAsync(this.header, 5, waitVal, 50)
        : null;
      if (waitRes && waitRes.async) {
        await waitRes.value;
      } else {
        await new Promise((resolve) => setTimeout(resolve, 1));
      }
    }
  }

  /**
   * Producer: Write a Uint8Array chunk into the ring buffer.
   * Yields asynchronously when the buffer is full, preventing event loop blocking.
   */
  async write(chunk) {
    let offset = 0;
    while (offset < chunk.length) {
      let free = this.capacity - Atomics.load(this.header, 1);
      if (free === 0) {
        free = await this.waitForSpace();
      }
      const wPtr = Atomics.load(this.header, 2);
      const toWrite = Math.min(chunk.length - offset, free, this.capacity - wPtr);
      this.data.set(chunk.subarray(offset, offset + toWrite), wPtr);
      Atomics.store(this.header, 2, (wPtr + toWrite) % this.capacity);
      Atomics.add(this.header, 1, toWrite);
      // Wake reader
      Atomics.add(this.header, 4, 1);
      Atomics.notify(this.header, 4);
      offset += toWrite;
    }
  }

  /** Producer: signal EOF. */
  close() {
    Atomics.store(this.header, 0, STATUS_CLOSED);
    Atomics.add(this.header, 4, 1);
    Atomics.notify(this.header, 4);
  }

  /** Producer or Consumer: signal error or cancellation. */
  abort() {
    Atomics.store(this.header, 0, STATUS_ABORT);
    Atomics.add(this.header, 4, 1);
    Atomics.add(this.header, 5, 1);
    Atomics.notify(this.header, 4);
    Atomics.notify(this.header, 5);
  }

  /**
   * Consumer: Synchronous read from ring buffer into `dest` array at `destOffset`.
   * Called inside Emscripten device driver read callback.
   * Blocks with Atomics.wait until data arrives, EOF is signaled, or abort.
   * Returns bytes read, or 0 on EOF.
   */
  read(dest, destOffset, maxLen) {
    while (true) {
      const status = Atomics.load(this.header, 0);
      if (status === STATUS_ABORT) throw new Error('Stream ring buffer was aborted');
      const avail = Atomics.load(this.header, 1);
      if (avail === 0) {
        if (status === STATUS_CLOSED) return 0; // EOF
        const waitVal = Atomics.load(this.header, 4);
        Atomics.wait(this.header, 4, waitVal, 50);
        continue;
      }
      const rPtr = Atomics.load(this.header, 3);
      const toRead = Math.min(maxLen, avail, this.capacity - rPtr);
      dest.set(this.data.subarray(rPtr, rPtr + toRead), destOffset);
      Atomics.store(this.header, 3, (rPtr + toRead) % this.capacity);
      Atomics.sub(this.header, 1, toRead);
      // Wake writer
      Atomics.add(this.header, 5, 1);
      Atomics.notify(this.header, 5);
      return toRead;
    }
  }
}
