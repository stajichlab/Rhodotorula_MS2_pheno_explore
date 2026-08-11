# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

## secreted_products

- **Path:** `analysis/secreted_products/` (`SECRETED_PRODUCTS.md`, `run.sh`, `scripts/01_secretion_analysis.py`)
- **Date:** 2026-07-02 · **Status:** complete (v1)
- **Question:** Which MS2 metabolites are uniquely & highly secreted (supernatant vs cell), and what biosynthetic gene families should we screen for in the genomes?
- **Inputs:** `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` (paired `C_*`/`SUP_*` columns), `input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz`
- **Method:** paired SUP-vs-C log2FC + Wilcoxon signed-rank + BH-FDR; cross-species Tau specificity; sign-flip permutation null; threshold sensitivity sweep.
- **Headline:** 292 paired strains (270 with genome link); 6,724 secreted, **130 uniquely secreted** features (all with MS2); top hits are lipid-like/glycolipid, strongly species-structured.
- **Key outputs:** `outputs/uniquely_secreted_features.csv`, `outputs/candidate_gene_families.md`, `outputs/figures/*.png`, `outputs/numbers.json`
- **Lineage:** uses the same raw MS2 table as the phase 0–3 color-phenotype pipeline, but exploits the C/SUP pairing those phases collapsed.

## secreted_products/sirius_annotation

- **Path:** `analysis/secreted_products/sirius_annotation/` (`run.sh`, `scripts/00_select_targets.py`, `scripts/01_export_mgf.py`, `scripts/02_run_sirius.sh`)
- **Date:** 2026-07-03 · **Status:** in-progress (v1, MGF export complete, SIRIUS formula/fingerprint/structure/canopus running)
- **Question:** What are the actual molecular formulas/structures/compound classes behind the 130 uniquely-secreted MS2 features (`secreted_products` open next-step #1)?
- **Inputs:** `analysis/secreted_products/outputs/uniquely_secreted_features.csv` (130 targets); `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` (for `adduct_rep_file`/per-sample peak areas); raw `mzML/C_*.mzML` and `ExFab_Supernatant/SUP_*.mzML` (parent dir, not in this repo).
- **Method:** positional-index join of targets to the aligned table → resolve one representative raw mzML per feature (`adduct_rep_file`, or max-peak-area fallback for 60/130 blank cases) → scan that mzML with `pyteomics` for the MS2 spectrum matching m/z (15 ppm) + RT (±0.15→0.5 min), keep highest-TIC match → MGF → SIRIUS `formula`/`fingerprint`/`structure`/`canopus`.
- **Headline:** 117/130 targets (90%) matched an MS2 spectrum and were exported (13 `no_ms2_match`, all in the ±0.15 window — widening to ±0.5 min recovered only 1 more); DDA duty cycle is the expected cause of the gap, not a script bug.
- **Key outputs:** `outputs/sirius_targets.csv`, `outputs/sirius_targets.mgf` (117 spectra), `outputs/mgf_export_summary.csv` (per-feature match audit), `outputs/sirius_project/` (SIRIUS results).
- **Caveats:** CSI:FingerID structure search and CANOPUS need `sirius login` (web service, needs internet) — not yet done as of v1; `formula`/`fingerprint` work offline, so molecular-formula/fragmentation-tree calls are unaffected by the login gap.
- **Lineage:** child of `secreted_products` — consumes its 130 uniquely-secreted feature list directly.

## phenotype_metabolite_association

- **Path:** `analysis/phenotype_metabolite_association/` (`PHENOTYPE_METABOLITE_ASSOCIATION.md`, `run.sh`, `scripts/01-07_*.py`)
- **Date:** 2026-08-11 · **Status:** complete (v1)
- **Question:** Are the strong pooled metabolite-color correlations reported in `docs/FEATURE_ANALYSIS.md` (legacy Phase 2 pipeline) real, once species/plate/C-SUP confounds are properly controlled for and validated with a permutation null, holdout replication, and a multivariate/module-level test?
- **Inputs:** `analysis/phase1_features_filtered.csv.gz`, `analysis/phase1_phenotype_data.csv.gz` (legacy Phase 0/1 outputs), `input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz` (for `ATTRIBUTE_species`/strain ID/C-SUP source, which the legacy pipeline's `Species` column lacked for ~55% of rows).
- **Method:** joint rank-space OLS partial-Spearman correlation (Species + Library Plate + sample_type regressed out at once, fixing a sequential-residualization bug and a mis-encoded-Plate bug in `scripts/phase2_correlation_analysis.py`); two-stage BH-FDR; strain-block permutation null (5,000 perms, respects C/SUP pairing); species-stratified strain-level 80/20 holdout; single-species (*R. mucilaginosa*) re-run + dedicated strain-level holdout to test for a within-species signal; linear (PLS) and non-linear (Random Forest) module-level tests with CV R² permutation-tested (hyperparameter-selection included in the null) to check for a joint many-feature signal; simulation-based power analysis (synthetic signal vs. real background nulls) for the Tier1 FDR criterion.
- **Headline:** **0/12,269** original Tier-1 hits survive correction (max corrected \|ρ\|=0.19 across all 22,023 tests); the 6 headline features named in `FEATURE_ANALYSIS.md` all fail permutation testing (17/17 pairs p>0.05); pooled holdout replication of the strongest corrected candidates is at chance (2/75, 2.7%); restricting to *R. mucilaginosa* alone (n=415, 210 strains, real a\*/b\* phenotype spread) does not recover a stronger signal (0 FDR hits, max \|ρ\|=0.196) and its own dedicated holdout also replicates at chance (3/75, 4.0%); neither the linear (PLS: CV R²=-0.09, p=0.56) nor non-linear (Random Forest: CV R²=-0.026, p=0.25) multivariate test finds a joint signal. **Power analysis caveat:** minimum detectable ρ at 80% power is ≈0.34 for both designs; power at ρ=0.20 (the largest ρ ever observed) is only 0-2%, so the null result confidently rules out the original ρ=0.7+ claims but is uninformative about a true ρ≈0.15-0.30 effect — underpowered, not disproven. The original pooled ρ up to 0.735 were almost entirely a species confound (Simpson's paradox).
- **Key outputs:** `outputs/corrected_tier1_hits.csv.gz`, `outputs/permutation_null_results.csv`, `outputs/holdout_summary.json`, `outputs/mucilaginosa_summary.json`, `outputs/mucilaginosa_holdout_summary.json`, `outputs/multivariate_module_test_summary.json`, `outputs/random_forest_module_test_summary.json`, `outputs/power_analysis_summary.json`, `outputs/numbers.json`.
- **Lineage:** child of the legacy Phase 0-3 pipeline (`scripts/phase0-2_*.py`, docs in `docs/`); supersedes `docs/FEATURE_ANALYSIS.md`'s confidence claims (that file now carries a caveats addendum pointing here).
