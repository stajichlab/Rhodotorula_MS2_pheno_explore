# Phenotype-Metabolite Association (Corrected Re-Analysis)

**Date:** 2026-08-11
**Status:** complete (v1) — supersedes the confidence claims in `docs/FEATURE_ANALYSIS.md`
**Parent analysis:** `scripts/phase0_batch_assessment.py` → `phase1_feature_filtering.py` →
`phase2_correlation_analysis.py` (legacy pipeline, pre-dates mycelium; see
`analysis/README.md`). This analysis reuses Phase 0/1's cleaned feature matrix and
phenotype table but replaces Phase 2's correlation method.

## Question

`docs/FEATURE_ANALYSIS.md` reported very strong pooled Spearman correlations (ρ up to
0.735) between individual MS2 features and CIELab color phenotypes (L*, a*, b*), and
named these as candidate brightness/pigment-controlling metabolites. Can we trust those
numbers, and if not, is there a real (smaller, harder-to-find) signal underneath?

## Why the original numbers are not trustworthy

Two problems, both root-caused by reading the actual code rather than the docs:

1. **Species was never regressed out.** `phase0_decision.json` set
   `"strategy": "stratified_with_plate"` because Phase 0 detected a highly significant
   species effect (F=32.65, p=1.1e-16). `phase2_correlation_analysis.py` only adds
   Species as a covariate when `'pooled' in decision['strategy']`, a branch that is
   never taken for this run — so every ρ in `FEATURE_ANALYSIS.md` is a **pooled,
   species-confounded** correlation across 16-17 species. `docs/PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md`
   already found 0% significant within-species hits in *R. mucilaginosa* (n=205, the
   dominant species) using a smaller top-200-feature scan — consistent with this being
   a between-species (Simpson's-paradox-style) artifact, not confirmed until now.
2. **The `spearman_partial_corr()` helper had a residualization bug**: it looped over
   covariates but reassigned (not accumulated) the residual each iteration, so with
   >1 covariate only the *last* one was actually regressed out. This never surfaced
   because the run only ever passed one covariate (Library Plate, itself mis-encoded
   as a continuous variable via `pd.factorize` rather than one-hot). It would have
   silently under-corrected the moment Species was added.

A third, independent data-quality problem surfaced while building the fix (see
`scripts/01_prepare_data.py` docstring): the `Species` column in
`phase1_phenotype_data.csv.gz` is NaN for 321/590 samples — almost every `SUP_*` row —
so even a naive "just add Species as a covariate" patch to Phase 2 would have silently
dropped or miscoded most supernatant samples. `ATTRIBUTE_species` from the richer
`input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz` file is
used instead (only 30/590 NaN). One strain ID collision (`17-332Y-1`, 4 rows) was also
found and dropped (documented in that script).

## Method

`scripts/01_prepare_data.py` — merges strain/species/sample-type identity onto the
Phase-1 feature and phenotype tables; drops 40/590 samples for missing covariates or
the strain-ID collision (550 retained, 281 strains, 17 species).

`scripts/02_corrected_correlation.py` — partial Spearman correlation via joint
rank-space OLS residualization against **Species + Library Plate + sample_type (C/SUP)**
all at once (fixing both bugs above), vectorized across all 7,341 features. Same
two-stage BH-FDR and tiering thresholds as the original Phase 2, for direct comparison.

`scripts/03_permutation_null.py` — strain-block permutation null (5,000 permutations)
that respects the C/SUP within-strain pairing: paired strains' (C, SUP) phenotype-residual
pairs are permuted as a unit; unpaired strains are permuted within their own singleton
stratum. Tested on the six headline `FEATURE_ANALYSIS.md` features plus the top-25
nominal hits per phenotype from the corrected scan.

`scripts/04_replication_holdout.py` — species-stratified, strain-level 80/20 holdout
(species with <5 strains kept entirely in train). Reports (a) Spearman correlation of
ρ_train vs ρ_test across all features (global calibration) and (b) whether the
train-selected top-25-per-phenotype hits replicate (same sign + test p<0.05) in the
held-out strains.

`scripts/05_within_species_mucilaginosa.py` — the corrected pipeline re-run restricted
to *R. mucilaginosa* alone (n=415 rows, 210 strains — the only species in this dataset
with enough strains for this), dropping Species from the covariate set. Directly tests
whether a real within-species (candidate-genetic) signal exists once between-species
variation can't contribute to it.

`scripts/06_multivariate_module_test.py` — a module-level test: PLS regression of the
same rank-residualized 7,341-feature matrix against residualized [L*, a*, b*], scored by
out-of-fold R² under `GroupKFold` (grouped by `strain_id`, same non-independence
handling as elsewhere). The number of PLS components is chosen by a small grid search
inside the CV, and — to avoid a repeat of the winner's-curse problem seen with the
"top nominal" univariate candidates — the *entire* grid-search procedure (not a single
fixed model) is re-run inside a strain-block permutation null (300 permutations). This
tests whether a joint, many-feature signal exists even if no single feature is
individually significant.

`scripts/07_mucilaginosa_holdout.py` — the direct within-species analogue of
`04_replication_holdout.py`: a plain strain-level 80/20 split restricted to
*R. mucilaginosa* (species-stratification isn't needed with only one species), same
global-calibration and top-hit-replication checks as script 04.

## Results

| Check | Result |
|---|---|
| Corrected pooled model, Tier1 (q<0.05) hits | **0** (vs. 12,269 in the uncorrected original) |
| Feature 2755 (headline "master brightness metabolite"), corrected ρ | **0.024** (was 0.735); permutation p=0.69 |
| Feature 6926, corrected ρ | 0.083; permutation p=0.14 |
| Feature 6188, corrected ρ | 0.088; permutation p=0.12 |
| All 6 headline features (2755/6926/6188/1560/5740/2308) × 3 phenotypes, permutation p<0.05 | **0/17 pairs** |
| Max \|ρ\| anywhere in the corrected 22,023-test scan | 0.19 |
| Holdout: Spearman(ρ_train, ρ_test) across all features | 0.10 (L*), **-0.43** (a*), 0.20 (b*) — near zero to negative |
| Holdout: train-selected top-25/phenotype hits replicating in test | **2/75 (2.7%)** — chance level |
| *R. mucilaginosa*-only (n=415, 210 strains) corrected model, Tier1 hits | **0**; max \|ρ\|=0.196 |
| *R. mucilaginosa* phenotype spread | L\* CV=2.3% (range 67.4-78.8); a\* CV=13.3% (0.6-14.4); b\* CV=34.7% (-0.3-8.7) — a\*/b\* have real spread, not a floor effect |
| Multivariate/module-level test (PLS, best of 2/5/10 components, 4-fold GroupKFold by strain) | observed CV R²=**-0.09**; null 95% range=[-0.14, -0.03]; **permutation p=0.56** — indistinguishable from chance |
| *R. mucilaginosa*-only holdout (strain-level 80/20, n=210 strains) | Spearman(ρ_train,test): L\*=0.19, **a\*=-0.46**, b\*=0.15; top-25/phenotype hit replication **3/75 (4.0%)** — chance level, same pattern as the pooled holdout |

**Bottom line:** once Species, Library Plate, and C/SUP sample type are jointly and
correctly controlled for, **none of the original "high-confidence" metabolite-color
associations survive** — not FDR correction, not a permutation null, not out-of-sample
replication. Restricting to *R. mucilaginosa* alone (which has genuine phenotype
variance to work with, so this isn't a floor-effect issue) doesn't recover a
stronger within-species signal either; effect sizes stay in the same ρ≈0.2 ceiling
seen pooled, and its own dedicated strain-level holdout replicates at the same
chance rate (3/75, 4.0%) as the pooled holdout. The multivariate/module-level test
(script 06) rules out the remaining obvious alternative — that a joint, many-feature
signal exists even though no single feature is significant: observed CV R²=-0.09 sits
inside the permutation null's 95% range (p=0.56). The pooled ρ=0.7+ correlations in
`FEATURE_ANALYSIS.md` were almost entirely a between-species confound (Simpson's
paradox), compounded by a residualization bug and a missing/miscoded Species column.

This is a **negative result for this dataset/method as analyzed**, not proof that no
metabolite controls pigmentation — see Next Steps. Four independent checks (FDR,
permutation, holdout, and multivariate) now agree, which is a reasonably thorough case
that there's no detectable single- or joint-feature MS2 signal for color phenotype in
this dataset at this sample size, rather than a fragile null from any one test.

## Caveats of this re-analysis itself

- Univariate per-feature tests only; a combination of features (module/pathway-level
  signal) could still exist and wouldn't show up here. PLS-DA / sparse multivariate
  methods were flagged as future work in the original `analysis/README.md` and remain
  untried.
- The C/SUP `sample_type` covariate is additive only; feature×compartment interactions
  (a metabolite that matters only in the secreted fraction) aren't modeled.
- The "top nominal" permutation-test candidates in step 3 exhibit obvious winner's-curse
  behavior (75/92 candidates "significant" at perm_pval=0.0002, the permutation floor)
  because they were selected for having the smallest p-value out of 7,341 features to
  begin with — that is exactly why step 4's holdout replication (not step 3's
  permutation test alone) is the credible check, and it shows near-chance replication.
- Within-*R. mucilaginosa*, only 8/17 species had ≥5 strains for the holdout split in
  script 04 (species-level, not applicable to 05); the mucilaginosa-only model in
  script 05 has no dedicated holdout yet (see Next Steps).

## Framework for incorporating future phenotype replicates

More color-phenotype measurement replicates (additional imaging batches, repeat
plate-scans, or new strains) are expected later. To slot them in without re-deriving
this pipeline:

1. **Keep replicates in long format, never pre-averaged.** Extend `sample_design.csv`'s
   schema with two new columns: `replicate_batch` (e.g. an imaging-run ID or date) and
   keep `strain_id` as the grouping key — do not collapse multiple phenotype
   measurements of the same strain into a single mean before this pipeline sees them.
   Averaging early throws away the within-strain variance needed to (a) estimate
   measurement noise and (b) keep the permutation/holdout designs valid.
2. **The strain is the unit of non-independence, not the sample-type row.** The
   strain-block permutation in `03_permutation_null.py` and `05_within_species_mucilaginosa.py`
   already group rows by `strain_id`; adding a `replicate_batch` dimension only means a
   strain's block grows from ≤2 rows (C, SUP) to ≤2×n_replicates rows. The block-permute
   logic (`strain_rows[strain][sample_type]` → extend to
   `strain_rows[strain][(sample_type, replicate_batch)]`) generalizes directly as long
   as every strain has the same replicate structure, or the "paired/singleton stratum"
   split in the current code is generalized to "stratify by the strain's replicate
   signature" (a small code change, not a redesign).
3. **Once real replicates exist, prefer a mixed-effects model over the permutation
   workaround.** The permutation-block design here is a defensible substitute for a
   proper mixed model *because* there wasn't enough replication to fit one reliably.
   With repeated measurements per strain, `statsmodels.regression.mixed_linear_model.MixedLM`
   (already available in this environment — confirmed via `python3 -c "import
   statsmodels.regression.mixed_linear_model"`) with `groups=strain_id` and a
   random intercept (or slope, if replicate count per strain is large enough) becomes
   preferable: it directly estimates a strain-level variance component and gives proper
   standard errors, rather than relying on permutation to launder the same problem.
   The covariate set (Species, Library Plate, sample_type) carries over unchanged as
   fixed effects.
4. **Register a `phenotype_source` version tag.** When new replicate data lands, treat
   it as a new value of `replicate_batch`, not a new file to swap in — keep
   `input_data/` immutable (per `CLAUDE.md`) and append. Re-run `01_prepare_data.py`
   with the extended metadata; everything downstream is parametrized off
   `sample_design.csv` and needs no other changes except the mixed-model upgrade in
   point 3 once there's enough replication to justify it (rule of thumb: ≥3 strains
   with ≥3 replicates each before a random-slope model is worth attempting over a
   random-intercept one).

## Next Steps

1. ~~**Multivariate/module-level test**~~ — done (script 06): PLS regression against
   L*/a*/b*, CV R² permutation-tested including the component-count selection step;
   p=0.56, no evidence of a joint many-feature signal either.
2. **SIRIUS structural annotation** — still not attempted for any feature from this
   analysis, and arguably no longer well-motivated: after FDR, permutation, holdout, and
   multivariate tests all agree on a null result, there currently isn't a defensible
   candidate feature list to annotate. Fixing the SIRIUS login/walltime issues in
   `analysis/secreted_products/sirius_annotation/` remains useful for that analysis's
   own (secretion-based, not phenotype-correlation-based) target list.
3. ~~**Holdout replication within *R. mucilaginosa***~~ — done (script 07): strain-level
   80/20 split, 3/75 (4.0%) top-hit replication, same chance-level pattern as the pooled
   holdout in script 04.
4. Once additional phenotype replicates arrive, revisit with the mixed-model upgrade
   described above — measurement noise may currently be swamping small true effects,
   and a proper strain-level variance component would help separate the two. Given four
   independent negative checks, this (or a fundamentally different data type, e.g.
   transcriptomics/genotype rather than a second cut of the same MS2 data) is probably
   the more promising direction than further re-slicing this dataset.

## Key outputs

- `outputs/sample_design.csv`, `outputs/features_cleaned.csv.gz` — cleaned, merged
  sample table used by every downstream step.
- `outputs/corrected_all_correlations.csv.gz`, `outputs/corrected_tier1_hits.csv.gz`,
  `outputs/corrected_correlation_summary.json` — corrected pooled model.
- `outputs/permutation_null_results.csv`, `outputs/permutation_null_summary.json`.
- `outputs/holdout_calibration.csv`, `outputs/holdout_hit_replication.csv`,
  `outputs/holdout_summary.json`.
- `outputs/mucilaginosa_all_correlations.csv.gz`, `outputs/mucilaginosa_permutation_results.csv`,
  `outputs/mucilaginosa_phenotype_spread.csv`, `outputs/mucilaginosa_summary.json`.
- `outputs/multivariate_permutation_null.csv`, `outputs/multivariate_module_test_summary.json`.
- `outputs/mucilaginosa_holdout_calibration.csv`, `outputs/mucilaginosa_holdout_hit_replication.csv`,
  `outputs/mucilaginosa_holdout_summary.json`.
- `outputs/numbers.json` — registered reportable values (via `register_value`).
