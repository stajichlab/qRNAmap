# qRNAmap CLI hardening: execution ledger

- Plan: `docs/superpowers/plans/2026-11-20-qrnamap-cli-hardening.md`
- Branch: `harden-cli`
- Date: 2026-11-20

## Rulings

- R30: Ruling: the gene flags and the mapping summary are sidecar files; the count matrix is unchanged (user chose option A on 2026-11-20) — cost if wrong: a different layout later means a format change.
- R31: Ruling: the `--mem-gb` default is 8 (the user's option A was taken as the default when no answer had arrived; the user may change it) — cost if wrong: jobs for large inputs may need `--mem-gb` raised or the default lowered.
- R32: Ruling: defect D3 (no `$SCRATCH` in jobs submitted from a bare environment) is handled by `--tmpdir` for `run --slurm` plus a clearer message, not by a site-specific `module load` in the scripts (the module approach worked in one test job on one node r47 but is site-specific and prints a long Tcl trace) — cost if wrong: bare-environment submitters must pass `--tmpdir`.
- R33: Ruling: `strand-check` removes an old `strand_check.tsv` first, after the input checks (R24 order) — cost if wrong: none.
- R34: Ruling: `run --slurm --tmpdir` under `$SCRATCH` or `/scratch` warns (not an error); the help says the job creates the directory when its parent is writable — cost if wrong: none.

The full SDD progress ledger of this plan is kept in the git-ignored workspace and is appended below by the controller at the end.
