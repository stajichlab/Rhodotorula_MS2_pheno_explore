# Copper-AUC Metabolite Association

**Date:** 2026-08-11
**Status:** complete (v1)
**Parent data:** `data/metadata/rhodotorula_auc_copper/` (copper stress-response AUC phenotype,
ingested 2026-08-11), joined onto the same Phase 1 QC'd feature matrix
(`analysis/phase1_features_filtered.csv.gz`) used by `analysis/phenotype_metabolite_association/`.
**Methodology:** directly reuses and extends the confound-correction / permutation / holdout /
power-analysis / multivariate-module-test methodology validated in
`analysis/phenotype_metabolite_association/` for the color-phenotype re-analysis.

## Question

Does metabolite composition (cell or supernatant) predict copper stress-response AUC
(`mean_auc_rate`), once species and batch confounds are properly controlled for?

## Design: cell and supernatant as parallel tracks

Unlike the color-phenotype analysis (which pooled `C_*`/`SUP_*` rows with `sample_type`
as a covariate), this analysis runs **two fully independent tracks** — cell (`C_*`) and
supernatant (`SUP_*`) — through the entire pipeline from step 1, rather than pooling them
or picking one fraction first. Rationale: copper resistance could plausibly be driven by
either intracellular defense metabolites or secreted compounds, and there was no strong
enough prior to justify testing one fraction before the other (user decision, logged in
`.living/decisions.md`).

## Data preparation (`scripts/01_prepare_data.py`)

- Copper AUC source: `input_data/Rhodotorula_AUC_copper.20260811.csv.gz`, 275 rows.
- Excluded permanently (user decision — no reconciliation attempted): the 2 rows for
  `SAMPLE_NAME = "TFCN_17-332M-1"` (conflicting `mean_auc_rate`/`Strain ID` values both
  mapped to `C_190`/`SUP_190`).
- Excluded: 2 rows with `MS2_SAMPLE_Cell`/`MS2_SAMPLE_Supernatant` = literal `"No MS2 Data"`.
- Species comes from `ATTRIBUTE_species` in the extended strain-traits metadata (not the
  copper file's own `SPECIES` column, not phase1's sparse `Species` column) — same fix
  as the color-phenotype analysis, for the same reason (recurring incomplete/inconsistent
  species columns in this project's legacy files).
- 7 more rows per track dropped because that sample was already excluded upstream in
  Phase 1 QC (blanks/failed-QC samples never entered `phase1_features_filtered.csv.gz`).
- **Result: 264 strains, 17 species, in each of the cell and supernatant tracks.**

## Confound check (`scripts/02_confound_check.py`)

Is `mean_auc_rate` itself species/plate-structured, before any feature-correlation work?

| Track | Species (Kruskal-Wallis) | Library Plate |
|---|---|---|
| cell | H=61.88, **p=5.8e-10** | H=21.19, **p=9.6e-5** |
| supernatant | H=61.88, **p=5.8e-10** | constant (=1.0 for all 264 rows — an upstream data property, not a bug; excluded as a covariate for this track) |

Copper AUC is strongly species-structured in both tracks — the species correction below is
doing real, necessary work, not a formality. Library Plate is also a real confound in the
cell track, but happens to be entirely constant within the supernatant track's samples (all
recorded as plate 1.0), so it is dropped from that track's design matrix rather than
included as a degenerate covariate.

## Corrected correlation (`scripts/03_corrected_correlation.py`)

Joint rank-space OLS partial-Spearman correlation (species [+ plate for cell] regressed out
of both feature ranks and phenotype ranks in one step, not sequentially), two-stage BH-FDR,
same Tier1/2/3 thresholds as the color-phenotype analysis (|ρ|>0.30 & q<0.05 / |ρ|>0.25 &
q<0.05 / |ρ|>0.20 & q<0.10).

| Track | n features tested | max \|ρ_corrected\| | Tier1 | Tier2 | Tier3 |
|---|---|---|---|---|---|
| cell | 7,341 | 0.224 | 0 | 0 | 0 |
| supernatant | 7,333 (8 dropped, zero residual variance) | 0.214 | 0 | 0 | 0 |

**No feature clears FDR<0.05 in either track.**

## Power analysis (`scripts/04_power_analysis.py`) — run early, not as an afterthought

Explicit lesson carried over from the color-phenotype analysis: its own power analysis,
written last, revealed the "five converging null checks" result was only informative
above ρ≈0.34 — a fact that should have shaped how the earlier null results were reported.
Run here immediately after the primary correlation scan instead.

Simulation-based: synthesize a feature correlated with the real phenotype residual at a
grid of true ρ, re-run the two-stage BH-FDR against the real background null p-values from
that track's own 7,341-feature scan, and measure the fraction of trials clearing Tier1.

| Track | Minimum detectable ρ at 80% power |
|---|---|
| cell (n=264) | **0.340** |
| supernatant (n=264) | **0.345** |

The largest ρ actually observed in either track (0.22) is well below this floor — power at
ρ=0.22 is low (consistent with the ~5% power at ρ=0.20 seen in the color-phenotype analysis
at a similar sample size / testing burden). **The null correlation result is uninformative
about a true ρ in the 0.15–0.30 range; it only confidently rules out anything ≥~0.34.**

## Permutation null (`scripts/05_permutation_null.py`)

Simpler than the color-phenotype analysis's strain-block permutation: each track already
has exactly one row per strain (no C/SUP pairing within a track), so a plain row permutation
is already strain-safe. Top-25 nominal candidates per track, 5,000 permutations each.

All 25/25 candidates in both tracks individually clear perm_p<0.05 — expected for the
smallest nominal p-values out of 7,341 tests (this confirms the individual p-values aren't
degenerate, it does **not** imply they survive multiple-testing correction; they don't, per
the FDR result above).

## Holdout replication (`scripts/06_holdout_replication.py`)

Species-stratified, strain-level 80/20 holdout (species with <5 strains kept entirely in
train). Two checks: global rank calibration (train ρ vs. test ρ across all features) and
replication of the top-25 train-nominal hits in the held-out test set (same sign AND test
p<0.05).

| Track | Calibration (Spearman, train ρ vs test ρ) | Top-25 hit replication |
|---|---|---|
| cell | **-0.226** (negative) | **0/25** |
| supernatant | 0.031 (~zero) | **0/25** |

Both tracks replicate at chance or worse — no evidence the nominal top hits carry real,
out-of-sample signal.

## Multivariate module tests (`scripts/07_multivariate_module_test.py`, PLS; `scripts/08_random_forest_module_test.py`, RF)

Is there a joint many-feature signal even though no single feature clears FDR? Both use
GroupKFold CV (grouped by strain_id) with the hyperparameter grid search (PLS n_components;
RF max_depth) repeated inside the permutation null itself, so the reported p-value accounts
for the selection procedure (nested-design fix validated in the color-phenotype analysis).

**PLS (linear), completed:**

| Track | Observed best CV R² | n_components | Null R² (mean, 95% range) | Permutation p |
|---|---|---|---|---|
| cell | -0.0496 | 10 | -0.217 [-0.379, -0.058] | **0.017** |
| supernatant | -0.0673 | 10 | -0.226 [-0.410, -0.074] | **0.020** |

Subtle result, worth reading carefully: the observed CV R² is **negative in both tracks**
(the model does not predict copper AUC well in absolute out-of-fold terms — a negative R²
means worse than predicting the mean). But it is permutation-significant in both tracks
(p<0.05): the real data's R² is *less negative* than R² from label-shuffled data, meaning
a PLS model fit to the real feature-phenotype relationship overfits/generalizes worse than
a null baseline by *less* than a model fit to pure noise does. This is consistent with a
weak joint multivariate signal across many features — smaller than any single feature's
own detectable effect, and far too weak to actually predict copper AUC usefully, but
statistically distinguishable from no relationship at all. This is a genuinely different
result from the color-phenotype analysis's own PLS test (CV R²=-0.09, permutation p=0.56,
not significant) — do not equate the two.

**Random Forest (non-linear), completed:**

| Track | Observed best CV R² | max_depth | Null R² (mean, 95% range) | Permutation p |
|---|---|---|---|---|
| cell | 0.0005 | None | -0.062 [-0.128, 0.0003] | **0.030** |
| supernatant | -0.0288 | 6 | -0.052 [-0.116, 0.003] | 0.199 (not significant) |

Random Forest **only partially corroborates** the PLS result: the cell track is again
permutation-significant (p=0.030 — real CV R² ≈0.0005, essentially flat but still
distinguishable from the null's more-negative distribution), while the supernatant track
is **not** significant with RF (p=0.199, observed R² falls within the null's 95% range) —
unlike PLS, which was significant in both tracks.

## Headline

Consistent with the color-phenotype analysis's univariate finding: no single feature, in
either the cell or supernatant fraction, survives FDR correction; permutation testing
individually confirms the top nominal p-values are not degenerate (not the same as
surviving correction); species-stratified holdout replication is at or below chance. The
power analysis shows this design (n=264/track) can only rule out associations ≥ρ≈0.34 —
the observed max |ρ|≈0.22 in both tracks means the univariate null result is
**underpowered, not disproven**, for a true modest effect, exactly mirroring the
color-phenotype analysis's own calibrated conclusion.

**Multivariate picture is more nuanced.** The **cell** track shows a weak but
permutation-significant joint signal under **both** PLS (p=0.017) and Random Forest
(p=0.030) — consistent across two different model classes, which is the kind of
corroboration this project's convention treats as meaningful (see the color-phenotype
analysis's own "don't trust a single multivariate test class" practice). The
**supernatant** track is significant under PLS (p=0.020) but not Random Forest (p=0.199)
— a single-model-class result, which per that same convention should be treated with more
caution (could be a PLS-specific artifact, e.g. sensitivity to a particular linear
combination that RF's tree splits don't capture the same way, or a false positive from
running 2 tracks × 2 model classes = 4 tests without correction across them).

**Bottom line:** there is tentative evidence for a weak, diffuse joint metabolite signal
associated with copper AUC in the **cell (intracellular)** fraction — too weak for any
single feature to be individually identifiable, and far too weak to be practically
predictive (CV R² ≈0), but statistically distinguishable from chance in two independent
model classes. The supernatant evidence is weaker/model-dependent. Neither should be
reported as a discovery without further validation (e.g., an independent holdout
specifically for the multivariate signal, or narrowing to a smaller, chemistry-motivated
feature set the way the color-phenotype analysis's successful torularhodin lead was
found).

## Key outputs

Per track (`outputs/cell/`, `outputs/supernatant/`): `corrected_all_correlations.csv.gz`,
`corrected_tier1_hits.csv.gz`, `permutation_null_results.csv`,
`holdout_calibration_all_features.csv.gz`, `holdout_hit_replication.csv`,
`multivariate_permutation_null.csv`, `random_forest_feature_importances.csv`,
`random_forest_permutation_null.csv`. Shared: `outputs/confound_check_summary.json`,
`outputs/corrected_correlation_summary.json`, `outputs/power_analysis_curve.csv` +
`power_analysis_summary.json`, `outputs/permutation_null_summary.json`,
`outputs/holdout_summary.json`, `outputs/multivariate_module_test_summary.json`,
`outputs/random_forest_module_test_summary.json`.

## Cross-reference with the genome-side copper signal (`scripts/09_cross_reference_genome_signal.py`)

The actual "connect metabolites to secreted-protein/genome predictions" goal this analysis
was built toward: do strains whose **cell-fraction metabolite profile** better predicts high
copper AUC (out-of-fold PLS/RF prediction from the F-004/F-005 models, in species-corrected
rank-residual space) also carry the elevated **genome methionine usage** found independently
in `Rhodotorula_Rodeo` (`.living/findings/metal-resistance.md` F-002)?

Joined on normalized strain name (241/264 cell-track strains matched to
`Rhodotorula_Rodeo`'s 251-strain genome-Met cohort — the cell track's own `strain_id`
column drops the `TFCN_`/`NRRL_` prefix that the genome side's `LOCUSTAG` keeps, so the
join uses `SAMPLE_NAME` instead; a sanity check confirmed both repos' independently-loaded
copies of `mean_auc_rate` agree exactly for the joined rows).

| Model | ρ (out-of-fold metabolite prediction vs. genome Met residual) | Permutation p |
|---|---|---|
| PLS | +0.077 | 0.240 |
| Random Forest | +0.050 | 0.438 |

**No correlation** (n=241 common strains, both null). The weak diffuse cell-fraction
metabolite signal (F-004/F-005) and the genome-side methionine signal (Rodeo F-002) do
**not** appear to be carried by the same strains — they are two independent, weak
signals, not two views of one underlying biology. This is a real, informative negative
result, not a failure to find something: it argues against a simple "high-Met strains
also have the distinctive metabolite profile" story, though it doesn't rule out that
each signal reflects something real on its own (both were already individually weak and
provisional).

## Next steps

- Resolve the `TFCN_17-332M-1`/`C_190`/`SUP_190` conflict with Christian Ona if a
  reconciled value becomes available, and re-run.
- The cross-reference above tested the specific hypothesis "same strains drive both
  signals" — it did not test whether specific metabolite features correlate with Met
  usage directly (bypassing the copper-AUC-prediction step), which remains open if a
  more exploratory link is wanted.
