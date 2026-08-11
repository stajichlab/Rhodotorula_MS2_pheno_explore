## What was worked on

Built a corrected re-analysis of the metabolite-color phenotype correlations reported in
`docs/FEATURE_ANALYSIS.md`, at the user's request to establish confidence in the top
hits before further annotation work.

- Found the legacy `scripts/phase2_correlation_analysis.py` never actually regressed
  out Species (its `'pooled' in decision['strategy']` branch was never taken) despite
  Phase 0 flagging a significant species effect, plus a residualization bug and a
  mis-encoded Plate covariate.
- New analysis `analysis/phenotype_metabolite_association/` (scripts 01-05): joint
  rank-space partial correlation (Species + Plate + C/SUP sample_type), strain-block
  permutation null, species-stratified holdout replication, and a single-species
  (*R. mucilaginosa*) re-run per the user's follow-up question about genetic control.
- Result: **0/12,269** original Tier-1 hits survive; all 6 headline features fail
  permutation testing; holdout replication is at chance (2/75); the mucilaginosa-only
  model finds no stronger within-species signal either. The original ρ up to 0.735 were
  almost entirely a species confound (Simpson's paradox).
- Logged 2 findings (`.living/findings/color-phenotype-metabolomics.md`), 1 decision, 2
  learnings; `docs/FEATURE_ANALYSIS.md` now carries an update section pointing to the
  corrected results; `analysis/ANALYSIS_MANIFEST.md` updated.
- Also diagnosed the SIRIUS annotation pipeline (`analysis/secreted_products/sirius_annotation/`):
  both prior runs failed at CSI:FingerID/CANOPUS with a missing `sirius login`, and the
  latest run was Slurm-cancelled mid-job; gave the user setup instructions.
- Committed (`5c9f4d9`) but **push failed** — no SSH agent/key available in this shell
  for `git@github.com`. User needs to push manually from a shell with their key loaded.

## Current state

- Branch: `main`, 1 local commit ahead of `origin/main` (not pushed).
- `analysis/secreted_products/sirius_annotation/outputs/sirius_project/` intentionally
  left untracked (68M SIRIUS intermediate output from an incomplete run).

## Next steps (see analysis doc for full list)

1. User to push the local commit (`git push`) once SSH auth is available.
2. Fix SIRIUS login, rerun with more walltime.
3. Consider multivariate/module-level test (PLS-DA) since no single feature clears
   univariate testing even within R. mucilaginosa.
4. Holdout replication specifically within R. mucilaginosa (script 05 has permutation
   but no train/test split yet).
