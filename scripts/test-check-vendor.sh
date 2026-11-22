#!/usr/bin/env bash
# Test scripts/check-vendor.sh. Run from the repository root. No BASH_SOURCE is used.
set -euo pipefail
root="$PWD"
real="$root/site/vendor/kallisto"
check="$root/scripts/check-vendor.sh"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
fail() { echo "FAIL: $1" >&2; exit 1; }

# 1. Real directory passes.
bash "$check" "$real" >/dev/null || fail "real directory should pass"

# 2. Usage errors give exit 2.
rc=0; bash "$check" >/dev/null 2>&1 || rc=$?; [ "$rc" -eq 2 ] || fail "no argument should exit 2 (got $rc)"
rc=0; bash "$check" relative/dir >/dev/null 2>&1 || rc=$?; [ "$rc" -eq 2 ] || fail "relative path should exit 2 (got $rc)"
rc=0; bash "$check" "$tmp/none" >/dev/null 2>&1 || rc=$?; [ "$rc" -eq 2 ] || fail "missing directory should exit 2 (got $rc)"

expect_fail() { # name dir
  local rc=0
  bash "$check" "$2" >/dev/null 2>&1 || rc=$?
  [ "$rc" -eq 1 ] || fail "$1 should exit 1 (got $rc)"
}

# 3. One changed byte.
cp -r "$real" "$tmp/changed"
printf 'X' | dd of="$tmp/changed/kallisto.wasm" bs=1 seek=100 conv=notrunc status=none
expect_fail "changed byte" "$tmp/changed"

# 4. Extra file.
cp -r "$real" "$tmp/extra"
echo extra > "$tmp/extra/extra.txt"
expect_fail "extra file" "$tmp/extra"

# 5. Missing file.
cp -r "$real" "$tmp/missing"
rm "$tmp/missing/kallisto.mjs"
expect_fail "missing file" "$tmp/missing"

echo "ok: check-vendor tests passed"
