# Species-Specific Feature Analysis — Summary

**Date:** 2026-07-02  
**Script:** `scripts/species_specific_features.py`  
**Input:** `results/phase0/phase0_ms2_aligned.csv.gz` (16,332 features × 590 samples)  
**Species source:** `input_data/growth_phenotype_summary_YPD2.csv.gz` (Strain ID → Species)

---

## Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Detection threshold | 0 | Peak area > 0 = detected |
| Min target prevalence | 0.80 | ≥80% of target-species samples must detect the feature |
| Max other prevalence | 0.20 | ≤20% of any other species' samples may detect it |
| Min species samples | 3 | Target species needs ≥3 samples for statistical reliability |
| FDR alpha | 0.05 | Benjamini-Hochberg correction threshold |
| Statistical test | Fisher's exact test (one-sided, greater) | Tests whether detection rate is significantly higher in target species |

---

## Results

**Total features tested:** 16,332  
**Species with sufficient samples (≥3):** 13  
**Candidate feature-species pairs (meeting prevalence criteria):** 2  
**Species-specific features (FDR q < 0.05):** 2

### Species-Specific Features

| Feature Index | Row ID | m/z | RT (min) | Target Species | N Samples (Target) | Detection Rate (Target) | Detection Rate (Others Max) | Detection Rate (Others Mean) | Specificity Score | N Detected (Target) | N Detected (Others) | Fisher p | Fisher q |
|---------------|--------|-----|----------|----------------|--------------------|-----------------------|---------------------------|------------------------------|-------------------|---------------------|---------------------|----------|----------|
| 11369 | 20534 | 1360.17 | 2.63 | R. sp. clade XI | 4 | 1.0000 | 0.1000 | 0.0018 | 0.9000 | 4 | 1 | 1.23e-9 | 2.47e-9 |
| 15641 | 47208 | 445.31 | 4.63 | R. taiwanensis | 12 | 0.8333 | 0.1794 | 0.1533 | 0.6539 | 10 | 88 | 8.76e-7 | 8.76e-7 |

### Species-Specific Features by Species

| Species | N Samples | N Species-Specific Features |
|---------|-----------|----------------------------|
| Rhodotorula mucilaginosa | 418 | 0 |
| Rhodotorula sphaerocarpa | 12 | 0 |
| Rhodotorula toruloides | 20 | 0 |
| Rhodotorula diobovata | 18 | 0 |
| Rhodotorula taiwanensis | 12 | 1 |
| Rhodotorula sp. clade I | 10 | 0 |
| Rhodotorula paludigena | 20 | 0 |
| Rhodotorula sp. clade XI | 4 | 1 |
| Rhodotorula dairenensis | 16 | 0 |
| Rhodotorula glutinis | 4 | 0 |
| Rhodotorula pacifica | 6 | 0 |
| Rhodotorula kratochvilovae | 6 | 0 |
| Rhodotorula graminis | 6 | 0 |

---

## Interpretation

Only 2 of 16,332 features (0.01%) met the strict species-specificity criteria. This indicates that very few metabolite features are truly species-specific at the binary presence/absence level. Most features are shared across Rhodotorula species, which is biologically expected for closely related taxa.

### Feature Details

**Feature 11369 (m/z 1360.17, RT 2.63 min) — R. sp. clade XI specific**
- Detected in 100% of R. sp. clade XI samples (4/4)
- Only 1 other detection: R. sp. clade I (1/10 = 10%)
- Zero detection in all other 15 species
- Mean peak area (detected): 1,271,281
- This feature was filtered out by Phase 1's zero-median filter (prevalence = 5/590 = 0.85%)

**Feature 15641 (m/z 445.31, RT 4.63 min) — R. taiwanensis enriched**
- Detected in 83.3% of R. taiwanensis samples (10/12)
- Also detected at low rates in other species: R. mucilaginosa (17.9%), R. toruloides (10%), R. diobovata (5.6%), etc.
- Mean peak area (detected, R. taiwanensis): 683,547
- This feature was also filtered out by Phase 1's zero-median filter (median = 0 across all samples)

### Note on Phase 1 Filtering

Both species-specific features were removed during Phase 1 preprocessing:
1. **Zero-median filter** removed features where the median peak area across all 590 samples is 0
2. Features present in only a small fraction of samples (like species-specific features) have a median of 0 and are filtered out

This means the Phase 2 correlation analysis (which uses Phase 1 filtered data) does not include these species-specific features. The raw Phase 0 data is required to capture all species-specific features.

---

## Output Files

| File | Description |
|------|-------------|
| `results/species_specific/sample_metadata_clean.csv` | Clean sample metadata: MS2 sample ID → strain_id, species, sample_type, library_plate |
| `results/species_specific/species_specific_features.csv` | All candidate feature-species pairs meeting prevalence criteria (2 rows) |
| `results/species_specific/species_specific_significant.csv` | FDR-significant species-specific features (2 rows, same as above) |
| `results/species_specific/species_specific_summary.json` | Summary statistics and parameters |

---

## Suggested Next Steps

1. **Relax thresholds** — Lower min target prevalence to 0.50 or raise max other prevalence to 0.30 to identify "species-enriched" features that are not strictly species-specific
2. **Separate C_ vs SUP_ analysis** — Currently both cell pellet (C_) and supernatant (SUP_) samples are pooled. Analyzing them separately may reveal species-specific signals in one sample type that are diluted by pooling
3. **Abundance-based enrichment** — Instead of binary presence/absence, compare mean peak areas across species to find features that are much more abundant in one species even if present in others
4. **Cross-reference with Phase 2 correlations** — Check whether species-enriched features (under relaxed criteria) are also correlated with color phenotypes, which could indicate confounding between species and phenotype
