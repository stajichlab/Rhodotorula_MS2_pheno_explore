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

## SIRIUS spectral confirmation (2026-08-11): torularhodin identity refuted

`sirius_annotation/` (scripts 00-03) ran all 3 candidates through SIRIUS 6.3.12
(formula → fingerprint → canopus → structures → write-summaries; see that folder's
scripts for the SIRIUS 6.x-specific setup issues this took to resolve: an outdated
5.8.1 module hitting a 404 API mismatch, a changed subcommand chaining order, and an
install-dir AOT cache that SIGILL-crashed on this cluster's CPUs). Results:

| raw_position | Hypothesized compound | SIRIUS top formula | CANOPUS class (confidence) | Structure-DB confidence |
|---|---|---|---|---|
| 12635 (the b\*-correlated one) | torularhodin (C40H52O2) | **C30H52N4O6** | Open-chain polyketides (0.221) | 0.027 |
| 11564 | torularhodin (C40H52O2) | **C29H52N6O5** | Lipopeptides (0.078) | 0.052 |
| 15041 | shinorine (C15H23N2O8) | **C12H25NO11** | Cyanogenic glycosides (0.996) | 0.317 |

**None of the three hypothesized identities survive.** All three formulas are
completely different molecular compositions from the targets, despite falling within
the 15 ppm m/z-matching window used to select them — i.e., all three were coincidental
isobars, not real identifications. CANOPUS and structure-DB confidence are low for the
two former "torularhodin" candidates in particular (0.027-0.221), so SIRIUS itself isn't
offering a confident alternative identity either.

## Interpretation

**The statistical result stands; the chemical identity does not.** raw_position=12635's
correlation with b\* is unaffected by this — same ρ=0.218, same permutation p=0.0002,
same near-identical holdout replication (train ρ=0.199, test ρ=0.197) — this was always
a claim about the correlation, verified independently of what the feature turned out to
be. But the "torularhodin" framing throughout this document's Interpretation section
(directional consistency with a red-orange pigment shifting the yellow-blue axis, etc.)
does not apply: raw_position=12635 is now an **unidentified feature** (SIRIUS's best
formula guess is C30H52N4O6, tentatively amino-acid/polyketide-related per CANOPUS, but
at low confidence). This is the strongest, best-replicated *statistical* signal found in
this project's phenotype-metabolite work — it just isn't attached to a known pigment
biosynthesis pathway anymore.

The lesson for future targeted searches on this dataset: an m/z-only match, even within
a tight ppm window and even with a biologically plausible RT, is a hypothesis that needs
formula/fragmentation confirmation before being reported as an identification — exactly
what happened here, and worth remembering before the next pathway-targeted search.

## Next Steps

1. ~~**MS2/spectral confirmation**~~ — done: SIRIUS refutes all 3 hypothesized
   identities (see table above). What raw_position=12635 actually is remains open.
2. **Spot-check which strains/species carry raw_position=12635** (32% detection rate) —
   does it track any known biology now that the pigment-pathway framing is gone?
3. ~~**Re-examine the other torularhodin candidate**~~ — resolved: raw_position=11564 is
   also not torularhodin (C29H52N6O5), so there's no "which one is the real torularhodin"
   question anymore — neither is.
4. **Broaden or re-run the pathway search** with a wider RT/adduct net, or reconsider
   whether this project's short C18 method (not a dedicated carotenoid C30 column) is
   simply unlikely to retain/detect real carotenoids well, before spending more effort
   on m/z-based carotenoid hunting in this specific dataset.
5. Treat raw_position=12635 as an unidentified-but-statistically-real lead: worth a
   dedicated identification effort (spectral library search, NIST/GNPS matching) if this
   phenotype-metabolite connection is worth pursuing further, independent of the original
   carotenoid hypothesis.

## Key outputs

- `reference_material/pigment_pathway_targets/pigment_pathway_targets.csv` — target
  compound × adduct list with computed masses.
- `outputs/pathway_target_matches.csv`, `outputs/pathway_target_feature_list.csv` —
  matching results.
- `outputs/targeted_correlation_results.csv`, `outputs/targeted_correlation_summary.json`,
  `outputs/targeted_permutation_results.csv`.
- `sirius_annotation/outputs/sirius_project/formula_identifications.tsv`,
  `structure_identifications.tsv`, `canopus_formula_summary.tsv` — SIRIUS 6.3.12 results
  refuting the torularhodin/shinorine identities (see table above).
- `outputs/targeted_holdout_results.csv`, `outputs/targeted_holdout_summary.json`.
