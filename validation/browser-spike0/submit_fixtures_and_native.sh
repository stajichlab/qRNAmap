#!/usr/bin/env bash
# Records the two Task 2 job submissions of 2026-11-21 (job IDs 32356703 and 32356735).
# Written after the jobs ran. Not executed as a script. Its parameters match sacct for both jobs:
# partition short, 2 CPUs, 8G, 2 h.
# Usage: submit_fixtures_and_native.sh WORKTREE DATA
#   WORKTREE  absolute path of the worktree (pixi env at WORKTREE/.pixi/envs/default)
#   DATA      absolute path of /bigdata/stajichlab/jstajich/projects/qRNAmap_data
# Run `fixtures.sh fetch` on the login node first (it needs internet). No BASH_SOURCE.
# The 8G request was not a measured need: see the README (peak 111040K and 1231516K).
set -euo pipefail
W=${1:?WORKTREE}; D=${2:?DATA}
mkdir -p "$D/spike0/logs"
sbatch --parsable -p short -c 2 --mem=8G -t 02:00:00 -J s0fix -o "$D/spike0/logs/fixtures-%j.out" \
 --export=ALL,FIX=$D/spike0/fixtures,QRNAMAP=$D/target-spike0/release/qrnamap,PIXI_BIN=$W/.pixi/envs/default/bin,SCRIPTS=$W/validation/browser-spike0,AMB=$D/amb \
 --wrap "$W/validation/browser-spike0/fixtures.sh build"
# The native job was submitted after the fixtures job had finished (add --dependency=afterok:<id> to chain them).
sbatch --parsable -p short -c 2 --mem=8G -t 02:00:00 -J s0native -o "$D/spike0/logs/native-%j.out" \
 --export=ALL,FIX=$D/spike0/fixtures,OUT=$D/spike0/fixtures/native,PIXI_BIN=$W/.pixi/envs/default/bin,SCRIPTS=$W/validation/browser-spike0 \
 --wrap "$W/validation/browser-spike0/native_refs.sh"
