#!/usr/bin/env bash
# Assemble the static site in dist-site/: the landing page, the tool page, the copied modules, the vendored
# kallisto and gzip reader, and the web package. Run from the repository root. No BASH_SOURCE is used.
# Build the web package first: pixi run build-wasm-web (writes pkg-web/).
# SITE_CNAME (optional): write it to dist-site/CNAME.
set -euo pipefail
test -f pkg-web/qrnamap_wasm.js && test -f pkg-web/qrnamap_wasm_bg.wasm || { echo "pkg-web is missing: run pixi run build-wasm-web first" >&2; exit 1; }
test -f site/vendor/kallisto/kallisto.wasm || { echo "site/vendor/kallisto is missing" >&2; exit 1; }
rm -rf dist-site
mkdir -p dist-site/analyze dist-site/sheet dist-site/vendor dist-site/pkg-web
cp site/index.html site/style.css site/theme.js site/testsets.json dist-site/
# The test-set reads are not part of the site: the page fetches them from the S3 mirror of release testsets-v1
# (site/analyze/testsets.js). A local site/testsets/ left from an older checkout is not copied.
cp site/analyze/*.html site/analyze/*.js site/analyze/*.json dist-site/analyze/
cp site/sheet/*.html site/sheet/*.js dist-site/sheet/
cp site/gz.js site/persist.js site/reference-panel.js dist-site/
cp site/vendor/fflate.js site/vendor/fflate.LICENSE dist-site/vendor/
cp -r site/vendor/kallisto dist-site/vendor/kallisto
cp pkg-web/qrnamap_wasm.js pkg-web/qrnamap_wasm_bg.wasm dist-site/pkg-web/

# reads.js is not used by the page (see site/COPIED.md), so it is not copied to dist-site.

# Name every file the page loads with one stamp (a hash of the source content). GitHub Pages lets a browser keep a file for
# 10 minutes, so an old script could meet a new page. The source files hold the token __BUILD__ for the stamp.
STAMP=$(
  { find site -type f ! -name COPIED.md ! -path 'site/testsets/*' | LC_ALL=C sort | xargs cat
    cat pkg-web/qrnamap_wasm.js pkg-web/qrnamap_wasm_bg.wasm; } | sha1sum | cut -c1-10
)
# Stamp the page and script files (not the vendored kallisto files: they are hash-checked and hold no token).
# `sed -i.bak` is the form that works with GNU sed and with BSD sed (macOS), which needs an argument
# after -i. The copies are removed before the check for a leftover token.
find dist-site -type f \( -name '*.html' -o -name '*.js' -o -name '*.css' \) ! -path 'dist-site/vendor/kallisto/*' ! -path 'dist-site/pkg-web/*' \
  -exec sed -i.bak "s/__BUILD__/$STAMP/g" {} +
find dist-site -type f -name '*.bak' -delete
if grep -rIl "__BUILD__" dist-site; then echo "a __BUILD__ token is left in dist-site" >&2; exit 1; fi
echo "build stamp: $STAMP"
if [ -n "${SITE_CNAME:-}" ]; then echo "$SITE_CNAME" > dist-site/CNAME; fi
touch dist-site/.nojekyll
bash scripts/check-vendor.sh "$PWD/dist-site/vendor/kallisto"
ls -l dist-site
