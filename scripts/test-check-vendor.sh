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

# 6. Symlink as extra entry.
cp -r "$real" "$tmp/symextra"
ln -s kallisto.mjs "$tmp/symextra/link.mjs"
expect_fail "extra symlink" "$tmp/symextra"

# 7. Listed file replaced by a symlink to identical content.
cp -r "$real" "$tmp/symrepl"
cp "$real/LICENSE" "$tmp/target_license"
rm "$tmp/symrepl/LICENSE"
ln -s "$tmp/target_license" "$tmp/symrepl/LICENSE"
expect_fail "symlink replacing a listed file" "$tmp/symrepl"

# Helper: copy the real directory and edit PROVENANCE.md.
mod_prov() { # name  awk-or-sed program applied via function body in $2
  cp -r "$real" "$tmp/$1"
  eval "$2" < "$real/PROVENANCE.md" > "$tmp/$1/PROVENANCE.md"
}
row_index="$(grep -E '^[0-9a-f]{64}  runtime/index.mjs$' "$real/PROVENANCE.md")"
hash_index="${row_index:0:64}"

# 8. Identical duplicate row.
mod_prov dup "sed 's|^$row_index\$|&\n$row_index|'"
expect_fail "duplicate row" "$tmp/dup"

# 9. Alias row with ./ prefix.
mod_prov alias "sed 's|^$row_index\$|$hash_index  ./runtime/index.mjs|'"
expect_fail "./ alias row" "$tmp/alias"

# 10. Row with ..
mod_prov dotdot "sed 's|^$row_index\$|$hash_index  runtime/../runtime/index.mjs|'"
expect_fail ".. row" "$tmp/dotdot"

# 11. Malformed row inside the block.
mod_prov bad "sed 's|^$row_index\$|&\nnot-a-hash  runtime/x.mjs|'"
expect_fail "malformed row" "$tmp/bad"

# 12. Missing fences.
mod_prov nofence "grep -v '^\`\`\`\$'"
expect_fail "missing fences" "$tmp/nofence"

echo "ok: check-vendor tests passed"
