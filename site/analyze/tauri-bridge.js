// Tauri desktop IPC bridge for native capabilities and path execution.
// Safely detects whether running inside Tauri desktop app and provides typed wrappers.

export function isTauri() {
  return typeof window !== 'undefined' && Boolean(window.__TAURI__?.core?.invoke || window.__TAURI_INTERNALS__?.invoke);
}

export async function tauriInvoke(cmd, args = {}) {
  if (window.__TAURI__?.core?.invoke) {
    return window.__TAURI__.core.invoke(cmd, args);
  }
  if (window.__TAURI_INTERNALS__?.invoke) {
    return window.__TAURI_INTERNALS__.invoke(cmd, args);
  }
  throw new Error('Tauri IPC is not available in standard browser environment');
}

let cachedCaps = null;

export async function checkNative() {
  if (!isTauri()) return { available: false, reason: 'browser' };
  if (cachedCaps) return cachedCaps;
  try {
    cachedCaps = await tauriInvoke('check_native_capabilities');
    return cachedCaps;
  } catch (err) {
    console.warn('Failed checking native capabilities:', err);
    return { available: false, error: err?.message || String(err) };
  }
}

export async function pickFile(title, filterName, filterExts) {
  if (!isTauri()) return null;
  return tauriInvoke('pick_file_dialog', { title, filterName, filterExts });
}

export async function pickDirectory(title) {
  if (!isTauri()) return null;
  return tauriInvoke('pick_directory_dialog', { title });
}

export async function readLocalSheet(path, baseDir = null) {
  return tauriInvoke('read_local_sample_sheet', { req: { path, base_dir: baseDir } });
}

export async function prepareNativeIndex(params) {
  return tauriInvoke('prepare_native_index', { req: params });
}

export async function runNativeQuant(params) {
  return tauriInvoke('run_native_quant', { req: params });
}

export async function resolvePath(path, baseDir = null) {
  return tauriInvoke('resolve_path', { path, baseDir });
}

/** Pick several files with the native dialog. Returns absolute paths (empty when cancelled). */
export async function pickFiles(title, filterName, filterExts) {
  if (!isTauri()) return [];
  return tauriInvoke('pick_files_dialog', { title, filterName, filterExts });
}

/**
 * Write the kallisto index (page-built bytes, or the gzip of a prebuilt index with `gzip: true`) to disk for native kallisto.
 * The bytes go as the raw IPC body: a JSON number array (Array.from) of a 100 to 300 MB index is several GB of text.
 */
export async function writeNativeIndex(bytes, { gzip = false } = {}) {
  const body = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
  const options = { headers: { 'x-qrnamap-gzip': gzip ? '1' : '0' } };
  if (window.__TAURI__?.core?.invoke) return window.__TAURI__.core.invoke('write_native_index', body, options);
  return window.__TAURI_INTERNALS__.invoke('write_native_index', body, options);
}

/** Can the native kallisto load the index at `indexPath`? Returns { ok, detail }. */
export async function checkNativeIndex(indexPath) {
  return tauriInvoke('check_native_index', { indexPath });
}

/** The native index this machine built before for genome `key` (an accession) with the current kallisto, or null. */
export async function cachedNativeIndex(key) {
  return tauriInvoke('cached_native_index', { key });
}

/**
 * Build a native kallisto index from the transcript FASTA bytes (the raw IPC body), kept under genome `key`. Used when the
 * native kallisto cannot read the page's kallisto.wasm index (index files are not portable between kallisto builds).
 */
export async function buildNativeIndex(fastaBytes, key) {
  const body = fastaBytes instanceof Uint8Array ? fastaBytes : new Uint8Array(fastaBytes);
  const options = { headers: { 'x-qrnamap-key': key || 'upload' } };
  if (window.__TAURI__?.core?.invoke) return window.__TAURI__.core.invoke('build_native_index', body, options);
  return window.__TAURI_INTERNALS__.invoke('build_native_index', body, options);
}

/** Stop the native kallisto run that is going now. */
export async function cancelNativeQuant() {
  if (!isTauri()) return false;
  return tauriInvoke('cancel_native_quant');
}

/** File extensions offered by the native FASTQ pickers. */
export const FASTQ_EXTS = ['fastq', 'fq', 'gz'];
