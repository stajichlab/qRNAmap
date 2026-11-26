# Provenance of the vendored kallisto build

This directory holds a WebAssembly build of kallisto 0.52.0 and the runtime files that load it in a browser worker.
Vendored on 2026-11-21.

## Sources

- Upstream repository: https://github.com/comparativechrono/genomics-wasm
- Upstream pin: commit `0b015f63e8048eb57b72ffd2706b79b803862b7c` (release `v0.1.0-rc.1`).
- kallisto version: 0.52.0.
- kallisto source revision: `4e9f29cf3b021260415430c057a22469ca081391`.
- kallisto source URL: https://codeload.github.com/pachterlab/kallisto/tar.gz/4e9f29cf3b021260415430c057a22469ca081391
  (SHA-256 of the tarball, from `sources.lock.json`: `e742c36bf8b48b81c18444e212c93ac7ba6e315363c6930e0fb2f0d27923ff56`).
- Emscripten version: 3.1.64 (from `sources.lock.json`).
- zlib version: 1.3.1.
- Our build run: https://github.com/stajichlab/genomics-wasm/actions/runs/44939577353
  (fork `stajichlab/genomics-wasm`, branch `qrnamap-spike0`, built from the upstream pin).
- Artifact name: `quant-compiled-kallisto`.

## How the files were obtained

- `kallisto.mjs` and `kallisto.wasm` come from `dist/kallisto/` in the artifact above.
  Their SHA-256 values equal the values in `build-manifest.json` (`assets`) of the same artifact.
- `browser.patch`, `sources.lock.json`, `build-manifest.json` and the licence files `LICENSE`, `LICENSE.bifrost`,
  `LICENSE.zlib` and `THIRD_PARTY_LICENSES.txt` also come from that artifact directory.
  `LICENSE`, `LICENSE.bifrost` and `THIRD_PARTY_LICENSES.txt` are byte-identical to `tools/kallisto/` at the pin.
- `LICENSE.BooPHF`, `LICENSE.CRoaring` and `LICENSE.libdivide` come from `tools/kallisto/` at the pin.
  They were not in the artifact directory.
- `LICENSE.genomics-wasm` is the root `LICENSE` (MIT) of the upstream repository at the pin.
- `runtime/index.mjs`, `runtime/worker.mjs` and `runtime/execute.mjs` come from `runtime/` at the pin.
  The spike harness (`validation/browser-spike0/parity.mjs`) imports `runtime/index.mjs` (`ToolRunner`).
  `index.mjs` loads `worker.mjs`, and `worker.mjs` imports `execute.mjs`. No other runtime file is needed for kallisto.
- Files were read and copied. None of them was run on the cluster.
- The kallisto source tarball and the zlib tarball from the artifact are not vendored (about 7.5 MB).

## Local changes

- `runtime/index.mjs` has one local change: the list of protocols a module of the tool may be loaded
  from now also holds the protocol of the file itself (`ALLOWED_PROTOCOLS`). The upstream list is
  `file`, `http` and `https` only, which refused the desktop app: its page is served over the
  app's own protocol (`tauri:`), and the tool module and the tool worker are its siblings, so they
  have that protocol too.
  It also supports passing `devices` (virtual block devices for streaming ring buffers) into the worker.
- `runtime/execute.mjs` registers virtual block devices via `FS.registerDevice` and `FS.mkdev` when `devices` are passed,
  allowing kallisto to stream reads on-the-fly through `StreamRingConsumer`.
  The changes are in the SHA-256 list below.

## What is not verified

- UNVERIFIED: whether the artifact equals a fresh upstream build. We did not rebuild it.
  Our build ran in a fork, from the upstream pin. We did not compare it with an upstream release asset.
- UNVERIFIED: that `browser.patch` and the kallisto source in the artifact match the licence texts in this directory.
  We copied the licence files that upstream ships. We did not audit the sources.
- UNVERIFIED: the runtime files were copied from the GitHub contents API at the pin. We did not compare them with the copy in the spike run.

- Fact: no Emscripten runtime licence file is included here. `kallisto.mjs` is Emscripten output. `THIRD_PARTY_LICENSES.txt` has no Emscripten entry.
- The check script trusts `PROVENANCE.md` itself. A change to `PROVENANCE.md` together with a change to a file is not detected by the script. Review of the diff is the control.

## Licences

See `THIRD_PARTY_LICENSES.txt` and the `LICENSE.*` files here.
kallisto and Bifrost are BSD-2-Clause. CRoaring is Apache-2.0. BooPHF is MIT. libdivide, kseq (MIT), wyhash (Unlicense) and zlib 1.3.1 are listed in the notices.
The MIT licence of the upstream genomics-wasm repository (`LICENSE.genomics-wasm`) covers its build, runtime and validation code.

## SHA-256 of every vendored file (except this file)

`scripts/check-vendor.sh <absolute dir>` recomputes these values and compares them.
It exits 1 on a mismatch, a missing file or an extra file.

```
efd8a2f8f2b64a250a9ecebe8a2b09d5c1bf0aff21f51e3f24d1a31ad1221310  LICENSE
354ef1d76ae18655b0f98d11606075af248b25b026d02f8ef0ad3f73606f9ae9  LICENSE.BooPHF
5b223235c64f595576d15ba2d983ef14f8e9e3ec45aafcad689dc40008e1b29e  LICENSE.CRoaring
3cd11ff854bb520a4dec0b2d47f26334c8375943018b867db75659c7e2cb59c2  LICENSE.bifrost
de4f8aee5c06ddddd7f1e9168cea65b665fe1239d9d9114cfe1975431b3a2514  LICENSE.genomics-wasm
524ffe6eed9d379255605a01076282517f8a41657c137fe00073e2c5659f4584  LICENSE.libdivide
845efc77857d485d91fb3e0b884aaa929368c717ae8186b66fe1ed2495753243  LICENSE.zlib
442fb228e4af03fe0a5b251a969f1c1406247f246e66b1e4908962060a22b6bb  THIRD_PARTY_LICENSES.txt
bc48df6b597a8d378c0674d4771b9cd02ab2626c1297e45d6b9a91cd74582249  browser.patch
46c819c4cbf39bb58c6da337ef1cdb0ad3b712cd9349d4b800a292d598bf9469  build-manifest.json
9836e42fc304d7616ec73d89262498b31cc96e2b35d2d7460b64f358b823b9f1  kallisto.mjs
01ef635d8438b7e6ac9c74c9d04a89f83f6801959f799f9f657523c90d37adfc  kallisto.wasm
e727ce1205fc8b6b878d1603e4c0c6056b272016dc7951627de3b8cc9dfb6ac7  runtime/execute.mjs
df29f87184a7c4e9096b7518d7321279ab040916fbf2701a3cd996d277dc004f  runtime/index.mjs
19ac7d452ace47900eb93b7a65a55ff47d3cd5df33ddee35664b93cf7d9ddc4a  runtime/worker.mjs
d3fb3b68d9517d4b982786a3825219d0845a7896772126c4286f1e8d762dc38a  sources.lock.json
```
