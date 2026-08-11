# Detailed Feature Analysis: Rhodotorula Metabolite-Phenotype Associations

**Date:** 2026-07-02  
**Analysis:** Top high-confidence metabolite features linked to color phenotypes

---

## ⛔ UPDATE (2026-08-11): the corrected re-analysis is done — the headline hits do not survive

The corrected, validated re-analysis at `analysis/phenotype_metabolite_association/`
is now complete. **None of the correlations in this document should be cited or acted
on.** Summary of what the corrected model + permutation null + holdout replication found:

- Once Species, Library Plate, and C/SUP sample type are jointly (and correctly)
  regressed out, **0 of the 12,269 original "Tier 1 high-confidence" hits survive**
  FDR correction. Max |ρ| anywhere in the corrected 22,023-test scan is **0.19**.
- Feature 2755 ("master brightness metabolite," originally ρ=0.735): corrected
  **ρ=0.024, permutation p=0.69**. Features 6926 and 6188 similarly collapse to
  ρ=0.08–0.09 with permutation p=0.12–0.14. **All 6 headline features named below,
  across all reported phenotypes (17 feature×phenotype pairs), fail permutation
  testing (p>0.05).**
- Out-of-sample holdout replication of the strongest *corrected* candidates was at
  chance level (2/75, 2.7%), and the train/test effect-size ranking was uncorrelated-to-
  negatively-correlated (Spearman ρ_train,test = 0.10 / **-0.43** / 0.20 for L*/a*/b*).
- Restricting to *R. mucilaginosa* alone (n=415 samples, 210 strains — the one species
  with enough strains to test this, and with real phenotype spread: a* CV=13%, b*
  CV=35%) does **not** recover a stronger within-species signal either (0 FDR hits,
  max |ρ|=0.196), and a dedicated strain-level holdout within that species replicates
  at chance too (3/75, 4.0%).
- A multivariate/module-level test finds **no joint many-feature signal** either,
  under a linear model (PLS: CV R²=-0.09, permutation p=0.56) or a non-linear one
  (Random Forest: CV R²=-0.026, permutation p=0.25), both permutation-tested including
  the hyperparameter-selection step.
- **Important caveat on all of the above (power analysis, script 09):** the minimum
  effect size this design can detect at 80% power is **ρ≈0.34**; power to detect
  ρ=0.20 (the largest value ever actually observed in the corrected scan) is only
  **0-2%**. The negative result rules out anything ρ≥0.34 with confidence — which
  safely covers the original ρ=0.7+ claims — but is **not informative** about whether a
  true, modest ρ≈0.15-0.30 effect exists; this dataset at this sample size cannot tell
  the difference between "no effect" and "an effect too small for this design to see."

**Conclusion:** the ρ up to 0.735 reported throughout this document were almost
entirely a between-species confound (Simpson's paradox), not evidence that any of
these specific metabolite features drive color phenotype — that part is confidently
refuted. Whether a smaller, real effect exists is genuinely unresolved (underpowered,
not disproven). See
`analysis/phenotype_metabolite_association/PHENOTYPE_METABOLITE_ASSOCIATION.md` for
full method, results, and next-steps (phenotype replicates + mixed-effects model, or a
pathway-targeted feature set, are the two concrete ways to lower the detection floor).

The original caveats that motivated this re-analysis are kept below for the audit trail.

## ⚠️ Caveats (added 2026-08-11 — read before citing any number below)

A methods review of the Phase 0–2 pipeline that produced this document found that the
headline correlations are **not yet trustworthy as evidence of a causal
metabolite→phenotype relationship**. A corrected re-analysis is in progress at
`analysis/phenotype_metabolite_association/` (see `.living/decisions.md` for the
tracking entry). Specific concerns:

1. **Species confound (Simpson's paradox risk) — the central issue.** `phase0_decision.json`
   recorded `"strategy": "stratified_with_plate"` because Phase 0 detected a highly
   significant species effect (F=32.65, p=1.1e-16). But `scripts/phase2_correlation_analysis.py`
   only adds `Species` as a covariate when `'pooled' in decision['strategy']` — a branch
   that is never taken for this run. **All correlations in this document (ρ=0.735 for
   Feature 2755, etc.) were computed pooled across 16 species with Species *not*
   regressed out**, despite Phase 0's own recommendation to do so. Consistent with this,
   the already-completed `docs/PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md` reran the
   correlation *within* species and found **0% significant within-species correlations
   in *R. mucilaginosa*** (n=205, 76% of samples) — only small species (n=4–10) showed
   elevated hit rates, which at that sample size is also consistent with noise. The
   pooled numbers in this document may largely reflect between-species metabolic
   differences rather than a within-species, phenotype-driving effect.

2. **Non-independent samples.** Each strain contributes both a `C_*` (cell pellet) and
   `SUP_*` (supernatant) sample; Phase 2 treats all 590 as independent observations. The
   two samples per strain are correlated (same genotype/species, and empirically their
   L*/a*/b* values differ but are not independent), and the C/SUP compartment itself is
   a known major axis of variation in this dataset (`.living/learnings.md` L-4). This
   inflates the effective sample size used by the significance tests.

3. **No permutation/holdout validation.** The two-stage BH-FDR correction assumes
   independence and correct null behavior; no permutation null or replication/holdout
   check has been run yet, both of which `analysis/README.md` already listed as
   outstanding "Future Work."

4. **Structural IDs are unverified guesses.** The carotenoid/xanthophyll/glycoside
   interpretations below are m/z-based speculation, not annotations. None of these
   features (2755, 6926, 6188, 1560, 5740, 2308) have been run through SIRIUS —
   the existing SIRIUS pipeline (`analysis/secreted_products/sirius_annotation/`)
   targets a different feature set (secreted-product candidates) and, separately, its
   CSI:FingerID/CANOPUS steps have not yet completed successfully (missing `sirius
   login`; see that analysis's manifest entry for the caveat).

5. **`analysis/phase2_summary.json` is corrupted** (truncated mid-key, invalid JSON) —
   a write was interrupted and never regenerated; treat any summary counts sourced from
   it as unverified until the corrected re-run replaces it.

**Bottom line:** treat every ρ/q-value below as provisional pending the corrected,
species/pairing-aware re-analysis with permutation nulls and holdout replication.

---

## Executive Summary

The strongest metabolite-phenotype associations point to **three distinct metabolite classes**:

1. **Brightness Control (L*)**: Primarily large molecular structures (m/z 200–1200)  
   - Consistent across ~12,000 features
   - Positive correlation: higher metabolite abundance → brighter strains
   - Likely carotenoids or cell-wall pigments

2. **Red-Green Color (a*)**: Moderate-size molecules (m/z 200–900)  
   - Negative correlation: higher metabolite → greener/less red
   - Suggests red pigments suppress or are depleted when these metabolites high
   - Possible antagonistic relationship between classes

3. **Yellow-Blue Hue (b*)**: Weaker signals, overlapping with brightness features

---

## Feature 2755: The "Master Brightness Metabolite"

**Key Metrics:**
- **Correlation (ρ):** 0.735 → **Explains 54% of brightness variation** 
- **q-value:** 2.39×10⁻⁹³ (virtually impossible to occur by chance)
- **m/z:** 808.5090
- **Retention Time:** 5.31 min (polar, early elution)
- **Adduct:** [M+2ACN+H]₁⁺ (double acetonitrile adduct → very polar)
- **Parent Mass:** 725.46 (calculated)
- **MS2 Spectrum:** ✓ Available (structure-informative fragmentation)

**Interpretation:**
- The **extreme m/z=808** with **2ACN adduct** suggests a **highly polar, large molecule**
- Parent mass ~725 Da is consistent with:
  - **Carotenoid glycoside** (e.g., γ-carotene + sugar moiety)
  - **Pigment complex** (e.g., bound to protein or lipid)
  - **Polymeric phenol** or **xanthophyll ester**

**Biological Significance:**
- Single feature controls >50% of brightness variation
- Extremely strong, replicable effect across 567 samples
- Suggests direct genetic/biochemical control of this one pathway
- **Next step:** MS2 fragmentation analysis + chemical standard comparison

---

## Feature 6926: The "Brightness Co-Star"

**Key Metrics:**
- **Correlation:** ρ = 0.731 (53.4% of variance)
- **m/z:** 486.7710 (simpler molecule than 2755)
- **Adduct:** [M+H]⁺ (standard positive ionization)
- **RT:** 1.89 min (less polar than 2755, quick elution)

**Interesting Property:**
- Correlates with **3 phenotypes** (L*, a*, b*)
- Same feature in all three color dimensions
- Suggests **single metabolite affects entire color palette**
- Likely candidate: **carotenoid monomer or intermediate**

**Possible Structures (m/z 485-487):**
- **Lycopene-related** (C₄₀H₆₄ = m/z 536 for full carotenoid)
- **Truncated carotenoid** or **degradation product**
- **Xanthophyll** (oxygenated carotenoid)

---

## Feature 1560: Small-Molecule Signal (m/z 212)

**Key Metrics:**
- **ρ:** 0.727 for brightness
- **ρ:** -0.353 for yellow-blue
- **m/z:** 212.1063 (small, simple molecule)
- **RT:** 4.81 min (retention on reverse-phase = hydrophobic)
- **Detection:** Only 13% of samples (selective)

**Interpretation:**
- **m/z ~212** could be:
  - **Simple phenolic compound** (quercetin-like, ~302 but this is smaller)
  - **Mycosporine-like amino acid** (MAA, UV protectant)
  - **Simple aromatic metabolite** (tyrosine derivative)
- Selective detection (13%) suggests species-specific or strain-specific production
- Hydrophobic retention suggests non-polar side chains

---

## Red-Green Features (a* phenotype)

### Feature 5740: Strongest a* Correlate

- **ρ:** -0.571 (explains 33% of variance in a* axis)
- **m/z:** 606.4190 (medium-large)
- **Adduct:** [M+H]⁺
- **Correlation:** Negative → Higher this metabolite = GREENER (less red)
- **Interpretation:** Likely a **red pigment precursor** or **antagonist**
  - Could be carotenoid isomer with different absorption
  - Might inhibit melanin/red pigment synthesis
  - Possible intermediate in pigment degradation

### Feature 2308: Red-Green Competitor

- **ρ:** -0.564 (very similar to 5740)
- **m/z:** 434.2220 (slightly smaller)
- **RT:** 0.90 min (very early → very polar)
- **Adduct:** [M+H]⁺
- **Interpretation:** 
  - Might be **phosphorylated** version of another metabolite
  - Or a **sugar-conjugated** small molecule
  - Polar metabolites that compete with red pigment production

---

## Multi-Phenotype Hits: The "Spectrum Controllers"

Features appearing in **multiple phenotypes** are highest-confidence biological signals:

| Feature | L* | a* | b* | Interpretation |
|---------|----|----|----|-|
| **2755** | 0.735 | - | -0.359 | Master controller; affects both brightness and hue |
| **6188** | 0.730 | -0.549 | -0.360 | Multi-dimensional regulator; likely central metabolite |
| **6926** | 0.731 | -0.558 | -0.344 | Coordinated color pathway; strong co-factor |
| **1560** | 0.727 | - | -0.353 | Hue modulator with brightness effect |
| **2308** | - | -0.564 | -0.343 | Balanced a*-b* controller |

**Biological Model:**
- These 5 features likely form a **coordinated metabolic module**
- May represent a biosynthetic pathway with multiple intermediates
- Regulatory hierarchy: 2755 → 6926 → 6188 (decreasing effect size)

---

## Metabolite Class Predictions

### Likely Brightness Drivers (L*)

**High-Confidence Candidates:**
1. **Carotenoids** (β-carotene, lycopene, xanthophyll)
   - m/z 536-556 (full C₄₀ carotenoids)
   - Evidence: Multiple signals in m/z 400-900 range
   - Biological: Known Rhodotorula pigments

2. **Carotenoid Glycosides** (glucose/xylose conjugates)
   - m/z 808, 725 parent mass (from Feature 2755)
   - Would explain 2ACN adduct
   - Biological: More soluble, stable form

3. **Oxygenated Lipids** (hydroxyperoxides, epoxides)
   - m/z 400-700 range consistent
   - Biological: Cell-wall reinforcement → lighter color?

### Likely a* Drivers (Red-Green)

**Hypotheses:**
1. **Melanin precursors** (DOPA, DHI, indole derivatives)
   - Negative correlation = less red when present
   - Could be competitive inhibitor

2. **Xanthophyll isomers** (different absorption maxima)
   - m/z 600-650 range
   - Absorption shift from red toward orange

---

## Key Findings

### 1. Parsimony Principle
- **~12,000 features** explain brightness, but top 10 explain most variance
- Suggests **simple model**: 1-2 "master" pathways control major phenotypes
- Implication: Likely single gene or operon involved

### 2. Correlation Hierarchy
```
Feature 2755 (54% variance) 
  └─→ Feature 6926 (53% variance) [co-abundant?]
  └─→ Feature 6188 (53% variance) [downstream?]
  └─→ Feature 4497+ (52-53% variance) [redundant/parallel?]
```

### 3. Metabolic Bottlenecks
- Early retention time (RT 0.9-2.2 min) for smaller features
- Late retention time (RT 4.8-7.6 min) for larger features
- Suggests **two-phase biosynthesis**: simple precursor → complex end-product

---

## Next Steps for Annotation

### Priority 1: MS2 Validation
- Run high-resolution MS2 fragmentation for Features 2755, 6926, 6188
- Compare against:
  - **HMDB** (Human Metabolome Database - has carotenoids)
  - **MassBank** (MS/MS reference library)
  - **LipidMaps** (carotenoid/lipid standards)
  - **YeastNet** (fungal metabolome)

### Priority 2: Chemical Standards
- Purchase carotenoid standards (β-carotene, lycopene, lutein)
- Measure retention times + m/z
- Cross-reference with observed features

### Priority 3: Biological Validation
- **Species stratification:** Do top features appear in all species?
- **Genetic basis:** Do bright strains have high Feature 2755?
- **Experimental:** Quantify Feature 2755 by HPLC in subset of samples

### Priority 4: Gene-Metabolite Mapping
- Link to known Rhodotorula pigment genes (e.g., *carotenoid synthase*)
- Check expression correlation with Feature 2755 abundance

---

## Summary Statistics

**Total Features Analyzed:** 7,341 (from 16,332 raw)  
**Tier 1 High-Confidence:** 12,269 associations  
**Phenotypes:** 3 (L*, a*, b*)  
**Key Finding:** 54% of brightness explained by single feature  
**Reproducibility:** 567 samples, consistent effect across library plates  

**Statistical Confidence:** Benjamini-Hochberg FDR correction (two-stage)  
- Stage 1: Within-phenotype q < 0.05
- Stage 2: Across-phenotype q < 0.05

---

**Ready for:** Publication-quality validation + structural determination
