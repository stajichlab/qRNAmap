#!/usr/bin/env bash
# Build the WebAssembly package: scripts/build-wasm-pkg.sh <nodejs|web> <output dir>
# Run from the repository root (the pixi tasks do). No BASH_SOURCE is used.
# The wasm-bindgen CLI version must equal the wasm-bindgen version in Cargo.lock.
set -euo pipefail
target_kind="${1:?usage: build-wasm-pkg.sh <nodejs|web> <output dir>}"
out_dir="${2:?usage: build-wasm-pkg.sh <nodejs|web> <output dir>}"
target="${CARGO_TARGET_DIR:-$PWD/target}"

lock_version="$(awk '/^name = "wasm-bindgen"$/ {getline; gsub(/[^0-9.]/, "", $0); print; exit}' Cargo.lock)"
cli_version="$(wasm-bindgen --version | awk '{print $2}')"
if [ "$lock_version" != "$cli_version" ]; then
    echo "wasm-bindgen CLI $cli_version does not match Cargo.lock $lock_version" >&2
    exit 1
fi

cargo build --release --target wasm32-unknown-unknown -p qrnamap-wasm
rm -rf "$out_dir"
wasm-bindgen --target "$target_kind" --out-dir "$out_dir" "$target/wasm32-unknown-unknown/release/qrnamap_wasm.wasm"
ls -l "$out_dir"
