# Parallel Development Branches & Agent Architecture

This document describes the two parallel development tracks and agent roles established for `qRNAmap` to solve large-scale RNA-seq quantification without browser memory exhaustion.

---

## 1. Overview of the Two Tracks

| Attribute | Track 1: Wasm Ring-Buffer Streaming | Track 2: Desktop Native Execution |
|---|---|---|
| **Branch** | `stream-quant` | `desktop-native-exec` |
| **Working Directory** | `/Users/jstajich/projects/qRNAmap` (repo root) | `/Users/jstajich/.local/share/opencode/worktree/4f6acb/desktop-native-exec` |
| **Target Environment** | Web browsers (Chrome, Firefox, Safari) & WebAssembly | Tauri Desktop application (macOS, Linux, Windows) |
| **Execution Engine** | `kallisto.wasm` via Emscripten in a Web Worker | System native `kallisto` 64-bit binary via Tauri IPC |
| **Read Depth** | Preview or configurable depth with constant small memory footprint | Full depth (100% of all reads directly from disk) |
| **Memory Model** | SharedArrayBuffer / Atomics bounded ring buffer (~16 MB) | Native process address space, direct file streaming |

---

## 2. Track 1: Agent 1 (`stream-quant`)

### Goal
Provide a WebAssembly-compatible quantification pipeline where FASTQ files are decompressed and streamed **on-the-fly** into `kallisto.wasm` without holding uncompressed FASTQ records or entire files in memory.

### Architecture
- **In-Memory Ring Buffer (`site/analyze/stream-ring.js`)**:
  - Implements `StreamRingBuffer` using `SharedArrayBuffer` (or typed array fallbacks) with `Atomics.wait` and `Atomics.notify`.
  - Coordinates a producer (JavaScript decompression worker) and a consumer (`kallisto.wasm` Emscripten virtual file system).
- **Single-Pass Streaming (`site/analyze/fastq.js`)**:
  - Implements `streamFastqRecords()` using `fflate` streaming chunks.
  - Decompresses block-by-block and writes directly to the ring buffer, yielding to avoid worker starvation.
- **Emscripten Device Registration (`site/vendor/kallisto/runtime/execute.mjs`)**:
  - Mounts custom device nodes (`/dev/stream1`, `/dev/stream2`) via `FS.registerDevice` / `FS.mkdev`.
  - Maps synchronous C/C++ `read()` syscalls inside `kallisto` directly to `StreamRingBuffer.readSync()`.
- **Worker Pipeline (`site/analyze/worker.js`)**:
  - Manages the quantification worker lifecycle, memory reset, and streaming cancellation.

### How to Run & Test Track 1
```bash
cd /Users/jstajich/projects/qRNAmap
# Run all unit tests
node --test validation/browser-page/tests/units.test.mjs
# Run cargo workspace tests
cargo test --workspace
# Build static web distribution
bash scripts/build-site.sh
```

---

## 3. Track 2: Agent 2 (`desktop-native-exec`)

### Goal
Allow the desktop app to bypass WebAssembly memory limits and single-thread constraints entirely. In desktop mode, `qRNAmap` can quantify full-depth, gigabyte-sized FASTQ datasets directly from local disk with multi-core performance.

### Architecture
- **Tauri IPC Commands (`apps/desktop/src/lib.rs`)**:
  - `check_native_capabilities`: Detects system `kallisto` binary version, available CPU threads, and temp directories.
  - `read_local_sample_sheet`: Parses local sample sheets, resolves relative FASTQ paths against the sheet directory or project root, and verifies file existence on disk.
  - `prepare_native_index`: Prepares native kallisto indexes (from prebuilt `.kidx`, compressed files, or builds from FASTA).
  - `run_native_quant`: Spawns native `kallisto quant` with multi-threading (`-t <CPUs>`), streaming directly from disk FASTQs at full depth.
  - `pick_file_dialog` / `pick_directory_dialog`: Native OS file pickers.
- **Tauri Bridge (`site/analyze/tauri-bridge.js`)**:
  - Safe bridge layer detecting if `window.__TAURI_INTERNALS__` exists before making IPC invocations.
- **UI Mode Toggle & Limit Relaxation (`site/analyze/samples.js`, `site/analyze/index.html`)**:
  - When running in Tauri, displays execution mode options:
    - **Native Execution (Full depth)**: Default in desktop. Bypasses subsampling; all reads from disk are quantified (no 200,000 read limit).
    - **Browser Preview**: Subsamples into memory; limit relaxed up to 100,000,000 pairs.
  - Automatically loads and verifies sample paths from local sheets (e.g., `/Users/jstajich/projects/Dhingra/samples.txt`).

### How to Run & Test Track 2
```bash
cd /Users/jstajich/.local/share/opencode/worktree/4f6acb/desktop-native-exec
# Run desktop tests (includes real-data quant test)
cd apps/desktop && cargo test
# Launch desktop app in debug mode
cargo run
# Launch desktop app in optimized release mode
cargo run --release
```

---

## 4. Shared Improvements in Both Branches

Both branches incorporate these recent enhancements:

1. **Volcano Plot Collision Resolution & Pointer Lines (`site/analyze/results.js`)**:
   - **Vertical Repulsion**: Iterative layout eliminates overlapping text for top-ranked genes with close $p$-values and fold changes (maintains $\ge 14\text{ px}$ baseline gap).
   - **Horizontal Jitter / Stagger**: Alternating horizontal offsets for adjacent labels.
   - **Leader Lines**: Dashed pointer lines (`stroke-dasharray="2,2"`) connect displaced labels directly to the data point circles.
   - **SVG Text Halo**: Added `paint-order: stroke fill` with a background stroke to prevent visual interference from background grid lines.
2. **Sample Sheet Name Normalization**:
   - Fixed `s.sample` vs `s.name` resolution in the sheet loader to prevent empty sample names and avoid triggering `SheetError::BadName` ("contains forbidden characters").

---

## 5. Working Guidelines

- **Two Worktrees**: Keep changes for the browser streaming path on `stream-quant` (repo root) and desktop native changes on `desktop-native-exec` (worktree).
- **Git Push Protection**: Do **not** push either branch to `origin` without explicit user permission.
