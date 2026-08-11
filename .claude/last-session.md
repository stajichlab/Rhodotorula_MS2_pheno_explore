## What was worked on

Built and completed a corrected re-analysis of the metabolite-color phenotype
correlations reported in `docs/FEATURE_ANALYSIS.md`, at the user's request to establish
confidence in the top hits before further annotation work, then extended it with the
user's own follow-up questions (single-species genetics test, multivariate/module-level
tests, non-linear model).

- Found the legacy `scripts/phase2_correlation_analysis.py` never actually regressed
  out Species (its `'pooled' in decision['strategy']` branch was never taken) despite
  Phase 0 flagging a significant species effect, plus a residualization bug and a
  mis-encoded Plate covariate. Also found `Species` was NaN for ~55% of samples
  (almost all `SUP_*` rows) and one strain-ID collision.
- New analysis `analysis/phenotype_metabolite_association/` (scripts 01-08): joint
  rank-space partial correlation (Species + Plate + C/SUP sample_type), strain-block
  permutation null, species-stratified holdout replication, a single-species
  (*R. mucilaginosa*) re-run + its own dedicated holdout, a linear (PLS) multivariate
  module-level test, and a non-linear (Random Forest) module-level test — all with the
  same strain-block permutation design and, for the multivariate tests, with
  hyperparameter selection nested inside the permutation loop to avoid winner's-curse.
- **Final result: five independent checks all agree on a null.** 0/12,269 original
  Tier-1 hits survive; all 6 headline features fail permutation testing; pooled and
  mucilaginosa-only holdout replication are both at chance (2/75, 3/75); neither PLS
  (p=0.56) nor Random Forest (p=0.25) find a joint multivariate signal. The original
  ρ up to 0.735 were almost entirely a species confound (Simpson's paradox).
- Logged 4 findings (F-001 through F-004 in
  `.living/findings/color-phenotype-metabolomics.md`), 2 decisions, 2 learnings;
  `docs/FEATURE_ANALYSIS.md` now carries an update section pointing to the corrected
  results; `analysis/ANALYSIS_MANIFEST.md` updated.
- Regenerated `analysis/phase2_summary.json`, which was found truncated/corrupted
  (invalid JSON from an interrupted write).
- Also diagnosed the SIRIUS annotation pipeline (`analysis/secreted_products/sirius_annotation/`):
  both prior runs failed at CSI:FingerID/CANOPUS with a missing `sirius login`, and the
  latest run was Slurm-cancelled mid-job; gave the user setup instructions.
- Committed across several commits (latest local: see `git log`) but **push has
  repeatedly failed** — no SSH agent/key available in this shell for `git@github.com`.
  User needs to push manually from a shell with their key loaded.

## Current state

- Branch: `main`, several local commits ahead of `origin/main` (not pushed).
- `analysis/secreted_products/sirius_annotation/outputs/sirius_project/` intentionally
  left untracked (68M SIRIUS intermediate output from an incomplete run).
- `analysis/phenotype_metabolite_association/` is feature-complete for this round:
  scripts 01-08, `run.sh` reproduces all of it, `scilintr` clean.

## Next steps (see analysis doc for full list)

1. User to push the local commits (`git push`) once SSH auth is available.
2. Fix SIRIUS login, rerun with more walltime — but note there's currently no
   defensible candidate feature list from this analysis to annotate.
3. A formal power analysis (what effect size is even detectable at n≈210-550) is the
   main open question after five converging negative checks.
4. Longer-term: bring in phenotype replicates (mixed-model framework already documented
   in the analysis doc) or a different data type (genotype/transcriptomics) rather than
   further re-slicing this same MS2 dataset.
