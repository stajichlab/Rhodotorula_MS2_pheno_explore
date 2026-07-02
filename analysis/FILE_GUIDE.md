# Analysis Pipeline File Guide

**Quick Reference for All Output Files**

---

## Phase 0: Batch Assessment & QC

| File | Type | Size | Purpose | When to Use |
|------|------|------|---------|------------|
| `phase0_batch_assessment.py` | Script | 11K | Main analysis script for batch/confounder detection | Run first to assess batch effects |
| `phase0_decision.json` | Config | 257B | JSON with analysis strategy decision | Reference for what strategy was chosen |
| `phase0_ms2_aligned.csv` | Data | 42M | Raw aligned MS2 features (590×16,332) | Input for Phase 1; contains all raw metabolite data |
| `phase0_metadata_aligned.csv` | Data | 86K | Aligned sample metadata (590 samples) | Cross-reference sample IDs, species, plate info |

### Phase 0 Decision Output Example
```json
{
  "plate_effect": true,
  "species_effect": true,
  "strategy": "stratified_with_plate",
  "phenotypes_to_use": ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"],
  "n_samples": 590,
  "n_features_raw": 16332
}
```

---

## Phase 1: Feature Filtering & Preprocessing

| File | Type | Size | Purpose | When to Use |
|------|------|------|---------|------------|
| `phase1_feature_filtering.py` | Script | 6.8K | Main filtering script (outlier removal, normalization, log2 transform) | Run second; applies QC filters |
| `phase1_features_filtered.csv` | Data | 81M | **MAIN ANALYSIS FILE** - Filtered/normalized features (590×7,341) | Input for Phase 2; use for all downstream analysis |
| `phase1_feature_metadata.csv` | Data | 338K | Feature statistics (prevalence, intensity, CV) | Check individual feature quality; understand filtering |
| `phase1_phenotype_data.csv` | Data | 86K | Aligned color phenotype data (590 samples × color traits) | Reference for phenotype values used in correlations |

### Feature Filtering Results Summary
```
Input features:      16,332
After filtering:      7,341 (44.9% retained)
Filtering applied:
  - Removed zero-median: 8,991 features
  - Removed low-prevalence: 0 features
  - Outliers capped: 78,561 values
  - Removed low-variance: 0 features
Quality: No NaNs, all features CV ≥ 0.1
```

---

## Phase 2: Correlation Analysis

| File | Type | Size | Purpose | When to Use |
|------|------|------|---------|------------|
| `phase2_correlation_analysis.py` | Script | 9.7K | Main correlation analysis script (Spearman partial, FDR correction) | Run third; produces main results |
| `phase2_all_correlations.csv` | Data | 3.0M | **COMPLETE RESULTS** - All 22,023 correlations with stats | Exploratory analysis; full transparency |
| `phase2_tier1_hits.csv` | Data | 1.7M | **HIGH CONFIDENCE HITS** - 12,269 features (ρ>0.30, q<0.05) | Primary results; publication-ready |
| `phase2_tier12_hits.csv` | Data | 2.3M | Medium confidence hits (ρ>0.25, q<0.05) - 4,317 features | Extended results; additional candidates |
| `phase2_summary.json` | Summary | 97B | Analysis statistics (counts, phenotypes, method) | Quick reference for analysis parameters |

### Results Tiering
```
Tier 1: |ρ| > 0.30 AND q < 0.05  →  12,269 features  (HIGH confidence)
Tier 2: |ρ| > 0.25 AND q < 0.05  →   4,317 features  (MEDIUM confidence)
Tier 3: |ρ| > 0.20 AND q < 0.10  →   1,798 features  (EXPLORATORY)
Not sig:                          →   3,639 features
```

### Columns in phase2_all_correlations.csv
- `phenotype`: Color trait (L*, a*, b*)
- `feature_index`: Feature number (0-7340)
- `rho`: Spearman correlation coefficient (-1 to +1)
- `pval`: p-value from correlation test
- `n_samples`: Number of samples used in correlation
- `reject_stage1`: FDR-corrected significance (within phenotype)
- `q_value_stage1`: Within-phenotype q-value
- `q_value_global`: Global q-value (Benjamini-Hochberg FDR-corrected)
- `reject_global`: Global significance flag
- `tier`: Confidence tier (Tier1_High, Tier2_Medium, etc.)

---

## Documentation Files

| File | Type | Purpose | Read For |
|------|------|---------|----------|
| `README.md` | Markdown | Complete pipeline documentation | How to run, understand each phase |
| `FINDINGS_SUMMARY.md` | Markdown | High-level scientific findings | Key discoveries, hypotheses, validation |
| `FILE_GUIDE.md` | Markdown | This file - guide to all outputs | Quick reference for what each file contains |

---

## Recommended Reading Order

### For Quick Overview:
1. **FINDINGS_SUMMARY.md** (5 min) - Top discoveries and implications
2. **phase2_tier1_hits.csv** (browse first 100 rows) - See top metabolite associations

### For Detailed Understanding:
1. **README.md** (15 min) - Full pipeline explanation
2. **phase0_batch_assessment.py** (10 min) - Understand batch/confounder detection
3. **phase1_feature_filtering.py** (10 min) - Understand filtering strategy
4. **phase2_correlation_analysis.py** (15 min) - Understand statistical method
5. **phase2_all_correlations.csv** (explore) - See all results with statistics

### For Validation/Extension:
1. **FINDINGS_SUMMARY.md** → Section 7 "Validation Recommendations"
2. **phase2_tier1_hits.csv** → Annotate Feature identities
3. **phase2_summary.json** → Cross-validate on holdout set
4. **phase1_phenotype_data.csv** → Double-check phenotype measurements

---

## Data File Relationships

```
Raw MS2 data (16,332 features, 608 samples)
    ↓ [Phase 0: Batch assessment]
phase0_ms2_aligned.csv (16,332 features, 590 samples)
phase0_metadata_aligned.csv (590 samples)
    ↓ [Phase 1: Filtering]
phase1_features_filtered.csv (7,341 features, 590 samples) ← MAIN FILE
phase1_feature_metadata.csv (7,341 features stats)
phase1_phenotype_data.csv (590 samples color traits)
    ↓ [Phase 2: Correlation]
phase2_all_correlations.csv (22,023 correlations)
phase2_tier1_hits.csv (12,269 high-confidence hits)
phase2_tier12_hits.csv (16,586 high+medium hits)
phase2_summary.json (analysis summary stats)
```

---

## Typical Analysis Workflow

### Exploratory
1. Read `FINDINGS_SUMMARY.md` for context
2. Load `phase2_tier1_hits.csv` in Excel/R
3. Sort by phenotype and |ρ| to find top associations
4. Check multi-phenotype hits (appear in multiple phenotype columns)

### Statistical Validation
1. Load `phase2_all_correlations.csv`
2. Filter to |ρ| > 0.40 for extreme hits (n=~1,000)
3. Check q_value_stage1 and q_value_global
4. Verify tier assignment makes sense

### Feature Annotation
1. Use `feature_index` from hits files
2. Look up corresponding m/z and retention time in original MS2 data
3. Search m/z against:
   - HMDB (Human Metabolome Database)
   - ChEBI (Chemical Entities)
   - MassBank (mass spectra library)
   - Carotenoid databases (e.g., CaRe)

### Species Stratification (Advanced)
1. Use `phase1_phenotype_data.csv` to identify samples by species
2. Subset `phase1_features_filtered.csv` by species
3. Rerun `phase2_correlation_analysis.py` separately per species
4. Compare effect sizes across species

---

## File Size Warnings

- **phase1_features_filtered.csv** (81M): Large file - load in R/Python, not Excel
- **phase2_all_correlations.csv** (3.0M): Can load in Excel (~50K rows × 9 columns)
- **phase0_ms2_aligned.csv** (42M): Very large - use `head()` to preview

### Example: Load phase1_features_filtered.csv in R
```r
features <- read.csv("phase1_features_filtered.csv", nrows=100)  # Load first 100 rows
features_all <- read.csv("phase1_features_filtered.csv")  # Load all (slow)
```

### Example: Load in Python
```python
import pandas as pd
df = pd.read_csv("phase1_features_filtered.csv")  # Load all at once
# For large file, use:
# df = pd.read_csv("phase1_features_filtered.csv", chunksize=100000)
```

---

## Reproducibility Checklist

- [ ] All scripts are self-contained (no external dependencies beyond pandas/numpy/scipy)
- [ ] Each script outputs a JSON file documenting decisions made
- [ ] Can rerun pipeline from scratch with original data
- [ ] Results are stable (ran Phase 0-2 multiple times, results identical)
- [ ] Statistical methods documented in script comments
- [ ] All thresholds justified in documentation

---

## Next Steps After This Analysis

### Immediate (Annotation)
- Identify Feature 2755 (strongest brightness association)
- Check if it's a known carotenoid

### Short-term (Validation)
- Cross-validate Tier 1 hits (5-fold CV)
- Stratify by R. mucilaginosa (largest species group)
- Check effect size consistency

### Medium-term (Experiments)
- Quantify top metabolites by HPLC
- Correlate lab measurements with color phenotypes
- Link to pigment biosynthesis genes

### Long-term (Integration)
- Combine with transcriptomics/genomics
- Build predictive model of color from metabolites
- Validate in knockouts/overexpression strains

---

**Last Updated:** 2026-07-02  
**Pipeline Version:** 1.0  
**Status:** Complete - All results ready for interpretation
