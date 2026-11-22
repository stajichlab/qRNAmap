# Browser spike 0 results

This file holds the spike 0 results. This section covers licences only.

## Licences (spike 0e)

This section lists facts from licence texts. It gives no legal conclusion. The decisions are in the last list.

### Scope and sources read

- Repository: `comparativechrono/genomics-wasm`, pinned commit `0b015f63e8048eb57b72ffd2706b79b803862b7c`.
- Files read through `gh api` (read-only) at the pin: `LICENSE`, `THIRD_PARTY.md`, `README.md`, and for kallisto, HISAT2, featureCounts and DESeq2 the files in `tools/<tool>/`.
- Bundles that the spike used, listed with `unzip -l` or `find`, under `/bigdata/stajichlab/jstajich/projects/qRNAmap_data/spike0/baseline/`:
  - `quant/quant-compiled-kallisto/dist/kallisto/` (unpacked build output).
  - `ci/tool-packages/hisat2-2.2.2-wasm.zip`.
  - `rna-analysis/rna-analysis-packages/featurecounts-2.1.1-wasm.zip`.
  - `release/deseq2-package.zip` (contains one inner file, `deseq2-webr-1.52.0.zip`, about 438 MB).
- I did not run any file from the repository or the bundles.
- I did not check that the spike's bundles were built from the pinned commit. The kallisto `sources.lock.json` in the bundle has the same `revision` as the pin.

### Repository licence (as written)

- `LICENSE` (repo root), first line: "MIT License". Copyright line: "Copyright (c) 2026 genomics-wasm contributors".
- `LICENSE`, last paragraph: "This licence covers original build infrastructure, JavaScript runtime and validation code only. Upstream programs, libraries, and code derived from them retain their own licences."
- `README.md`, heading "Licence": "compiled packages are not all MIT-licensed."
- `THIRD_PARTY.md`, last paragraph: "The repository-level MIT licence applies only to original infrastructure, runtime and validation code."
- `THIRD_PARTY.md` heading "Upstream software", first sentence: the WebAssembly packages "contain upstream software, not reimplementations of the algorithms."
- `THIRD_PARTY.md`, last paragraph: "Each binary package includes the relevant notices, build recipes, portability changes and the corresponding upstream sources."
- GitHub API field `license` for the repository: key `other`, name "Other", `spdx_id` `NOASSERTION`.
- The texts say the repository MIT licence covers the repository's own code. They say it does not cover the tool programs. They do not say that the compiled bundles are MIT.
- The earlier note in the qRNAmap design that "genomics-wasm's own LICENSE is MIT" matches the root `LICENSE` text. The texts do not say why GitHub shows "Other".

### Per tool

#### kallisto 0.52.0

- Tool licence as written:
  - `tools/kallisto/README.md`, last paragraph: "kallisto and Bifrost use the BSD 2-Clause licence."
  - `tools/kallisto/LICENSE`, first line: "BSD 2-Clause License". Copyright line: "Copyright (c) 2017, Nicolas Bray, Harold Pimentel, Páll Melsted and Lior Pachter".
  - `THIRD_PARTY.md`, entry "kallisto": "BSD-2-Clause".
- Bundled library licences (from `tools/kallisto/THIRD_PARTY_LICENSES.txt`; headings are file names in that file):
  - `LICENSE.bifrost`: "BSD 2-Clause License", copyright "[2018], [Pall Melsted, Guillaume Holley]".
  - `LICENSE.CRoaring`: "Apache License Version 2.0, January 2004".
  - `LICENSE.BooPHF`: "The MIT License (MIT)", copyright 2015 Guillaume Rizk.
  - `LICENSE.libdivide`: "libdivide is made available under two licenses" (zlib and Boost).
  - `src/kseq.h`: "The MIT License", copyright Attractive Chaos.
  - `src/unordered_dense.h`: header says "ankerl::unordered_dense" version 3.1.1 and `SPDX-License-Identifier: MIT`.
  - `ext/bifrost/src/libpopcnt.h`: notice text, copyright 2016 - 2017 Kim Walisch.
  - `ext/bifrost/src/wyhash.h`: "released under the Unlicense".
  - File text, line 3: "HDF5, HTSlib and zlib-ng are not linked into this build."
  - zlib 1.3.1: pinned in `sources.lock.json`. The text "Zlib license is supplied separately as LICENSE.zlib" is in `THIRD_PARTY_LICENSES.txt`. `LICENSE.zlib` is not in `tools/kallisto/` at the pin. It is in the built `dist/kallisto/`.
  - Emscripten: `sources.lock.json` pins `"emscripten": "3.1.64"`. `tools/kallisto/` has no `LICENSE.emscripten` at the pin. The built `dist/kallisto/` has no `LICENSE.emscripten` either.
- Source versions (`tools/kallisto/sources.lock.json`): kallisto `0.52.0`, revision `4e9f29cf3b021260415430c057a22469ca081391`; zlib `1.3.1`; Emscripten `3.1.64`.
- Redistribution text: `THIRD_PARTY_LICENSES.txt` line 3: "The exact upstream source archive accompanies this package." `README.md`: "the exact source archives accompany each package."
- Files in the spike's built `dist/kallisto/`: `LICENSE`, `LICENSE.bifrost`, `LICENSE.zlib`, `THIRD_PARTY_LICENSES.txt`, `sources.lock.json`, `browser.patch`, `build-manifest.json`, `kallisto-4e9f29cf...tar.gz`, `zlib-1.3.1.tar.gz`, `kallisto.mjs`, `kallisto.wasm`, `serial-thread.h`. Separate `LICENSE.BooPHF`, `LICENSE.CRoaring` and `LICENSE.libdivide` files are not in that folder. Their texts are inside `THIRD_PARTY_LICENSES.txt`.
- No kallisto release zip was in the baseline folder. I only listed the unpacked build output.

#### HISAT2 2.2.2

- Tool licence as written:
  - `tools/hisat2/README.md`, last paragraph: "HISAT2 and the port are distributed under GPL-3.0-or-later."
  - `tools/hisat2/LICENSE`: "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007" (full text, 674 lines).
  - `THIRD_PARTY.md`, entry "HISAT2": "including code from Bowtie 2; GPL-3.0-or-later."
- Bundled component notices (`tools/hisat2/THIRD_PARTY_LICENSES.txt`):
  - Line 3: "Full upstream sources and GPL-3.0 license are included in this package."
  - `tinythread.h`: copyright 2010-2012 Marcus Geelnard, with a zlib-style permission text; modified by Nigel Dyer 2016.
  - `third_party/cpuid.h`: Free Software Foundation, with the "GCC Runtime Library Exception, version 3.1" (full exception text in the file).
  - Bowtie 2 code: "Copyright 2011, Ben Langmead", under "GNU General Public License ... version 3".
  - `qsufsort.c`: "Copyright 1999, N. Jesper Larsson"; "This software may be used freely for any purpose. However, when distributed, the original source must be clearly stated".
  - The file has no zlib, htslib or Emscripten section. `sources.lock.json` lists no zlib dependency.
- Source versions (`tools/hisat2/sources.lock.json`): HISAT2 `2.2.2`, revision `99583d7536b9ee017ac07de8834017a3bf99a2fe`; Emscripten `3.1.64`.
- Redistribution text: `README.md`: "The package includes that archive, the patch, dispatcher, license notices and build script."
- In `hisat2-2.2.2-wasm.zip`: the root `LICENSE` starts with "MIT License". `dist/hisat2/LICENSE` and `tools/hisat2/LICENSE` start with "GNU GENERAL PUBLIC LICENSE". The zip also holds `THIRD_PARTY.md`, `source-manifest.json`, `dist/hisat2/THIRD_PARTY_LICENSES.txt`, and `build/downloads/hisat2-99583d75...tar.gz`.

#### featureCounts 2.1.1 (Subread)

- Tool licence as written:
  - `tools/featurecounts/README.md`, heading "Sources and licences": "Subread is distributed under GPL version 3 or later."
  - Same section: "The build scripts and small port helpers use the repository's MIT licence; this does not change the licence of the compiled program."
  - `tools/featurecounts/LICENSE`: "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007".
  - `tools/featurecounts/THIRD_PARTY_LICENSES.txt`: "Licence: GNU General Public License, version 3 or later."
  - `THIRD_PARTY.md`, entry "featureCounts": "GPL-3.0".
- Bundled libraries:
  - zlib 1.3.1: "Licence: zlib licence. LICENSE.zlib accompanies the compiled distribution". `LICENSE.zlib` is not in `tools/featurecounts/` at the pin. It is in `dist/featurecounts/` in the bundle.
  - Emscripten 3.1.64: `LICENSE.emscripten` first line: "Emscripten is available under 2 licenses, the MIT license and the University of Illinois/NCSA Open Source License." `THIRD_PARTY_LICENSES.txt`: the `mailbox.js` override "derives from src/library_pthread.js" in Emscripten; "their notices remain in that SDK's source distribution."
- Source versions (`tools/featurecounts/sources.lock.json`): Subread `2.1.1` (`subread-2.1.1-source.tar.gz`); zlib `1.3.1`; Emscripten `3.1.64`.
- Redistribution text: `README.md`: "The distribution includes its source lock, build manifest, patches, complete Subread and zlib source archives, and licence notices."
- In `featurecounts-2.1.1-wasm.zip`: root `LICENSE` starts with "MIT License". `dist/featurecounts/LICENSE` starts with "GNU GENERAL PUBLIC LICENSE". The zip also has `dist/featurecounts/LICENSE.emscripten`, `LICENSE.zlib`, `THIRD_PARTY_LICENSES.txt`, `sources.lock.json`, and three source tarballs: `build/downloads/subread-2.1.1-source.tar.gz`, `build/downloads/zlib-1.3.1.tar.gz` and `build/rna-analysis-sources/hisat2-source.tar.gz`. It also has test-data licence files `datasets/rna-analysis/provenance/HISAT2_LICENSE` and `datasets/rna-analysis/provenance/yeast-UPSTREAM_LICENSE.txt`.

#### DESeq2 1.52.0 (webR, R 4.6.0)

- Files at the pin: `tools/deseq2/` has `README.md`, `sources.lock.json`, `licenses/Noto-OFL.txt`, `licenses/Noto-font-notices.txt`. It has no `LICENSE`, no `LICENSE*` file and no `THIRD_PARTY_LICENSES.txt` at the pin (not found at the pin).
- Tool licence as written:
  - `tools/deseq2/README.md`, last paragraph: "DESeq2 is LGPL-3.0-or-later. The webR binary distribution is GPL-3.0 and dependencies retain their individual licences."
  - Same paragraph: "The archive includes upstream licence notices, corresponding sources and build information. Original JavaScript and validation infrastructure in this repository are MIT-licensed."
  - `sources.lock.json`, package `DESeq2` 1.52.0: `"license": "LGPL (>= 3)"`.
  - `THIRD_PARTY.md`, entry "DESeq2": "LGPL-3.0-or-later."
- GPL parts, as written:
  - `THIRD_PARTY.md`, entry "webR": "distribution binaries are GPL-3.0. Build scripts and additional application software are MIT unless stated otherwise upstream. R and bundled libraries retain their respective licences."
  - `THIRD_PARTY.md`, entry "R": "GPL-2.0-or-later."
  - `sources.lock.json`, `r` section: `"license": "GPL (>= 2)"`; R source `R-4.6.0.tar.gz`.
  - `sources.lock.json`, package `locfit` 1.5-9.12: `"license": "GPL (>= 2)"`. `THIRD_PARTY.md`: "`locfit` is GPL-2.0-or-later; the build changes its compiler configuration without changing the statistical implementation."
  - `sources.lock.json` lists 43 packages, each with a `license` field. I read only DESeq2, locfit and R. I did not tabulate the other 41 licences.
- Other bundled licences: `tools/deseq2/licenses/Noto-font-notices.txt`: Noto fonts "licensed under the SIL Open Font License, Version 1.1". The inner zip has `dist/deseq2/webr/LICENSE.md`, a webR "Software Licensing" list (R GPL v2 and v3, PCRE2 BSD, XZ public domain, jsTree MIT, libgfortran GPL v2 and v3, Dragonegg GPL v2, LLVM Flang Apache 2.0, xterm.js MIT, and more; I read the first 8 entries).
- Source versions (`tools/deseq2/sources.lock.json`): DESeq2 `1.52.0`; R `4.6.0`; webR `0.6.0` (npm `webr-0.6.0.tgz` and `webr-0.6.0-source.tar.gz`); Emscripten source `5.0.7` appears under `build/deseq2/sources/` in the bundle.
- Redistribution text: see the README paragraph above.
- In `deseq2-package.zip`: one file, `deseq2-webr-1.52.0.zip`. That inner zip has 313 entries. Its root `LICENSE` starts with "MIT License" and matches the repository `LICENSE`. It has `THIRD_PARTY.md` (119 lines; it differs from the repository file at the pin, I did not diff the content), `dist/deseq2/webr/LICENSE.md`, `tools/deseq2/licenses/Noto-*.txt`, and 75 files under `build/deseq2/sources/` (including `R-4.6.0.tar.gz` and `webr-0.6.0-source.tar.gz`). It has no `LICENSE` file for DESeq2 itself and no `THIRD_PARTY_LICENSES.txt`.

### Summary table (texts as written)

| Tool | Version | Tool licence as written | Source archive in bundle |
|---|---|---|---|
| kallisto | 0.52.0 | BSD 2-Clause | yes (kallisto, zlib) |
| HISAT2 | 2.2.2 | GPL-3.0-or-later | yes (HISAT2) |
| featureCounts | 2.1.1 | GPL-3.0 or later | yes (Subread, zlib) |
| DESeq2 | 1.52.0 | LGPL-3.0-or-later; webR binaries GPL-3.0; R GPL-2.0-or-later | yes (75 source files, including R and webR) |

### What the texts say and do not say

- All four READMEs and `THIRD_PARTY.md` say each package ships notices and the upstream sources. No text states a rule about hosting the bundles on a web page or about fetching them from another site.
- No text in the files I read says how a user or a third party may redistribute the built bundles. The repository `LICENSE` says it does not cover upstream code.
- The texts give MIT for the repository code and the tool licences above for the tools. They do not say which licence applies to a zip that holds both. What I found in the bundles: HISAT2 zip and featureCounts zip each have an MIT root `LICENSE` and a GPL `LICENSE` in `dist/`. The DESeq2 inner zip has an MIT root `LICENSE` and `dist/deseq2/webr/LICENSE.md` (a webR licence list); it has no GPL `LICENSE` file. For kallisto, no zip exists in the baseline folder; the unpacked build output has `LICENSE` (BSD 2-Clause) and no MIT root `LICENSE`.
- `THIRD_PARTY.md` does not list the Emscripten runtime as its own entry. The featureCounts files name it (MIT and University of Illinois/NCSA). The kallisto and HISAT2 files at the pin do not name it.
- GitHub shows "Other" (`NOASSERTION`). The root `LICENSE` text says MIT.
- The DESeq2 directory at the pin has no tool `LICENSE` file. The README gives the licence statement.

### qRNAmap side (facts)

- qRNAmap has a `LICENSE` file at the repository root. Its text is "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007".
- `Cargo.toml` (workspace) has `license = "GPL-3.0-only"`. `README.md` line 4 says "Licence: GPL-3.0-only."
- `docs/superpowers/specs/2026-11-19-qrnamap-design.md` (lines 174 and 178) says: qRNAmap is GPL-3.0-only (user, 2026-11-19); qmutfind is GPL-3.0-only, verified in its `Cargo.toml`.
- The same design (line 210) says "A page that serves GPL binaries also offers the matching source archives, which the bundles include." It marks this as proposed.
- Size of the DESeq2 package zip: about 438 MB (inner zip 437,694,988 bytes). GitHub Pages size limits are not checked in this section.

### Questions for the user (decision)

1. Which licence does qRNAmap use for the page code that loads the tool bundles? The repo file says GPL-3.0-only now. Do you keep that?
2. May the qRNAmap page host (bundle) the tool bundles, or must the page fetch them from the genomics-wasm releases at run time?
3. If the page hosts them: do you want to ship the source archives from each bundle with the page, and where (same site, or a link)?
4. Which attribution text and notice files must the page show for kallisto, HISAT2, featureCounts, DESeq2, webR, R and the Emscripten runtime? Where in the page does it go?
5. Do you want to ask the genomics-wasm author to clarify the repository licence (GitHub shows "Other") and the status of the compiled bundles?
6. For DESeq2: the bundle has no tool `LICENSE` file at the pin, and 41 package licences are not tabulated here. Do you want a full licence table for the 43 R packages before release?
7. Should the project pin the bundles by release asset or by commit, given that the licence texts differ between files at the pin and in the zips (for example `THIRD_PARTY.md` in the DESeq2 zip)?
