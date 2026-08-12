# Pathway-Targeted Metabolite-Phenotype Association

**Date:** 2026-08-11
**Status:** complete (v1)
**Parent analysis:** `analysis/phenotype_metabolite_association/` — this analysis exists
because that one's power analysis (script 09, finding F-005) showed the untargeted
7,341-feature scan cannot detect true effects below ρ≈0.34, and the largest ρ actually
observed there (0.20) was in the ~0-2%-power range. A targeted, literature-motivated
feature list avoids most of that multiple-testing burden and can find smaller real
effects that the untargeted scan structurally cannot.

## Question

Given a small set of features selected **before** looking at the phenotype data, purely
on chemical grounds (known Rhodotorula pigment-pathway masses), is there a detectable
metabolite-color correlation that survives the same confound correction, FDR, permutation
null, and holdout replication used throughout `phenotype_metabolite_association/`?

## Method

1. **Target compound list** (`reference_material/pigment_pathway_targets/build_target_list.py`
   → `pigment_pathway_targets.csv`): masses computed from molecular formulas (not
   hand-typed literature values) for the carotenoid biosynthesis pathway (phytoene →
   phytofluene → ζ-carotene → neurosporene → lycopene/γ-/β-carotene → torulene →
   torularhodin — the canonical fungal carotenogenesis pathway, high confidence for this
   genus) plus, per user's explicit scope choice, DOPA-melanin precursors and
   mycosporine-like amino acids (both flagged lower confidence for this species). 4
   adducts per compound (`[M+H]+`, `[M+Na]+`, `[M+2H]2+`, `[M-H2O+H]+`) = 68 target rows.
2. **Matching** (`scripts/01_match_targets.py`): 15 ppm tolerance against **all 16,332**
   raw features (not just the 7,341 that survived Phase 1's generic prevalence/CV
   filters — see below for why that mattered), matching on observed adduct.
3. **Corrected correlation** (`scripts/02_targeted_correlation.py`): identical joint
   rank-space partial-correlation method as `phenotype_metabolite_association/scripts/02`
   (Species + Library Plate + sample_type regressed out at once), but pulling raw peak
   areas directly from the aligned table (since 2/3 matches aren't in
   `features_cleaned.csv.gz`) and normalizing by total-ion-signal across all 16,332
   features per sample (not the QC-survivor-only sum Phase 1 used, which isn't defined
   for excluded features). Two-stage BH-FDR over only 9 tests (3 features × 3 phenotypes).
4. **Permutation null** (`scripts/02`, same call): strain-block permutation (5,000
   permutations), identical design to `phenotype_metabolite_association/scripts/03`.
5. **Holdout replication** (`scripts/03_holdout_validation.py`): species-stratified
   strain-level 80/20 split, identical design to `phenotype_metabolite_association/scripts/04`.

## A methodologically important detour: matching only survivors found almost nothing

Matching against the 7,341 Phase-1-filtered features first (the "obvious" thing to try,
since that's the feature set every other script in this project uses) found only 1 hit —
a low-confidence MAA. Matching against the **full 16,332-row raw table** found 3,
including **two torularhodin candidates** (the diagnostic red Rhodotorula pigment), neither
of which is in the 7,341-feature survivor set. This makes sense in hindsight: Phase 1's
blanket ≥10%-of-590-samples prevalence filter will discard a pigment that's only produced
by a subset of species/strains — exactly the kind of feature this analysis is looking for.
**Lesson for future targeted analyses on this dataset: always search the full raw table,
never the QC-filtered subset.**

## Results

| Compound | Adduct | ppm error | RT (min) | Detection rate | Phenotype | ρ (corrected) | q (FDR) | Permutation p |
|---|---|---|---|---|---|---|---|---|
| **torularhodin** (raw_position=12635) | [M+H]+ | 11.51 | 6.36 | 32.0% | **b\*** | **0.218** | **4×10⁻⁶** | **0.0002** |
| torularhodin (12635) | [M+H]+ | 11.51 | 6.36 | 32.0% | a\* | 0.135 | 0.008 | 0.025 |
| torularhodin (12635) | [M+H]+ | 11.51 | 6.36 | 32.0% | L\* | -0.128 | 0.009 | 0.032 |
| shinorine (raw_position=15041) | [M+H]+ | 10.33 | 0.53 | 68.2% | b\* | -0.115 | 0.018 | 0.007 |
| shinorine (15041) | [M+H]+ | 10.33 | 0.53 | 68.2% | a\* | 0.098 | 0.043 | 0.021 |
| torularhodin (raw_position=11564) | [M+H]+ | 2.46 | 5.70 | 16.9% | all 3 | \|ρ\|≤0.044 | n.s. | p≥0.38 |

5/9 tests significant at q<0.05 and permutation p<0.05 (vs. 0/22,023 in the untargeted
scan's Tier1 criterion, and 0/17 for the untargeted scan's headline features under
permutation).

### Holdout replication (strain-level 80/20 split, `scripts/03`)

Binary replication criterion (same sign + test p<0.05): **0/9**. But this masks an
important distinction from the untargeted analysis's holdout failures, where effect
sizes typically flipped sign or magnitude entirely between train/test. Here:

| Phenotype | ρ (train, n=447) | p (train) | ρ (test, n=103) | p (test) |
|---|---|---|---|---|
| torularhodin (12635) × **b\*** | 0.199 | 3.5×10⁻⁵ | **0.197** | 0.060 |
| torularhodin (12635) × a\* | 0.157 | 0.001 | -0.097 | 0.357 |
| torularhodin (12635) × L\* | -0.143 | 0.003 | -0.102 | 0.332 |
| shinorine (15041) × a\* | 0.072 | 0.14 | 0.200 | 0.056 |
| shinorine (15041) × b\* | -0.142 | 0.003 | -0.052 | 0.625 |

**torularhodin × b\* is the standout: the train and test effect sizes are essentially
identical (0.199 vs. 0.197, same sign)** — the holdout test at n=103 simply doesn't have
enough power to push p below 0.05 for an effect of that size (consistent with the parent
analysis's own power-analysis logic: smaller n → higher detection floor). This is a much
stronger form of evidence than the binary pass/fail suggests. The other two torularhodin
associations (a*, L*) and the shinorine associations show weaker or sign-inconsistent
behavior between train/test and should be treated with more caution — a* in particular
flips sign.

## Interpretation

**Torularhodin (raw_position=12635) correlating with b\* (yellow-blue axis) is the
strongest, most biologically coherent, and best-replicated candidate found across this
entire project's phenotype-metabolite work.** It:
- Was pre-registered by chemistry (a known Rhodotorula pigment, not selected from the
  correlation results), avoiding the winner's-curse problem that undermined the
  untargeted scan's "top nominal" candidates.
- Survives FDR correction and a strain-block permutation null on the full 550-sample data.
- Shows near-identical effect size in an independent held-out set of strains, with the
  binary non-replication being a holdout-sample-size power issue, not a sign/magnitude
  inconsistency.
- Points in the chemically sensible direction: torularhodin is a red-orange pigment: b*
  is the yellow-blue CIELab axis, so this reads as "more torularhodin → shifts on the
  yellow-blue axis," directionally consistent with an orange/red carotenoid's
  contribution to perceived hue.

**This should still be treated as a strong lead, not a confirmed result**, because:
- It is one candidate out of a short list, in a dataset where 8/9 targeted tests were
  weaker or didn't hold up as cleanly — this isn't "everything replicated," it's "one
  specific, chemically-motivated candidate replicated unusually well."
- The identification is mass-only (11.5 ppm, wider than ideal, though within the 15 ppm
  tolerance); MS2 spectral confirmation against a torularhodin reference/database has not
  been done (see Next Steps).
- Only 32% detection rate — consistent with a strain/species-specific pigment, but also
  means the correlation is driven by a subset of samples; worth a spot-check of which
  species/strains carry it.
- This is a short (~9 min) reversed-phase run, not a dedicated carotenoid (C30) column;
  RT 6.36 min (in the later half of the gradient) is at least directionally consistent
  with a hydrophobic carotenoid, but retention behavior alone doesn't confirm identity.

## Next Steps

1. **MS2/spectral confirmation**: this feature has an MS2 spectrum (`has_ms2=True`) —
   run it through SIRIUS (once the login issue in
   `analysis/secreted_products/sirius_annotation/` is fixed) or compare fragmentation
   pattern manually against published torularhodin MS/MS data.
2. **Spot-check which strains/species carry it** (32% detection rate) — does it track
   with known torularhodin-producing species, or is it broader/narrower than expected?
3. **Re-examine the other torularhodin candidate** (raw_position=11564, tighter 2.46 ppm
   match but null correlation) — is it a genuine isomer/different ionization state, an
   isotope peak, or a chance mass coincidence? Worth resolving before treating 12635 as
   "the" torularhodin peak.
4. If confirmed, this reopens the genetics side (`docs/CAROTENOID_GENE_BRIGHTNESS_ANALYSIS.md`,
   `docs/BRIGHTNESS_GENETICS_STRATEGY.md`) with an actual metabolite anchor rather than
   gene copy-number alone.

## Key outputs

- `reference_material/pigment_pathway_targets/pigment_pathway_targets.csv` — target
  compound × adduct list with computed masses.
- `outputs/pathway_target_matches.csv`, `outputs/pathway_target_feature_list.csv` —
  matching results.
- `outputs/targeted_correlation_results.csv`, `outputs/targeted_correlation_summary.json`,
  `outputs/targeted_permutation_results.csv`.
- `outputs/targeted_holdout_results.csv`, `outputs/targeted_holdout_summary.json`.
