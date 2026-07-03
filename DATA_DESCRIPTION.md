# Data & Script Structure: Rhodotorula MS2-Phenotype Analysis

**Quick Reference:** 590 samples × 7,341 metabolite features → 22,023 correlations with 3 color phenotypes

---

## Input Data (Raw, Read-Only)

Located: `input_data/`

### 1. **Rhodotorula_MS2_aligned_features_ms2.csv.gz** (18 MB)
**Source:** MS2 metabolomics data from mass spectrometry analysis  
**Shape:** 16,332 features × 608+ samples (includes blanks, QC standards)

**Key Columns:**
| Column | Type | Example | Notes |
|--------|------|---------|-------|
| `row ID` | String | "FT0001" | Unique feature ID from peak detection |
| `row m/z` | Float | 808.5090 | Mass-to-charge ratio (signature) |
| `row retention time` | Float | 5.31 | Liquid chromatography retention (min) |
| `row ion mobility` | Float | — | Drift time (if available) |
| `row CCS` | Float | — | Collision cross-section |
| `Gaussian_similarity` | Float | 0.95 | Peak shape quality |
| `charge` | Int | 1 | Ionization charge state |
| `adduct` | String | "[M+H]1+" | Ionization adduct type |
| `parent_mass` | Float | 725.46 | Calculated neutral mass |
| `has_ms2` | Bool | True | MS/MS fragmentation available |
| `detection_count` | Int | 567 | # samples with this feature |
| `C_1.mzML Peak area`, `C_2...` | Float | 1000–1M | Feature abundance in each sample |

**Samples:**
- Prefixes: `C_*` (strain samples), `SUP_*` (supernatant), `QC_*` (quality control), `Blank_*` (blanks)
- Total analytical samples: ~590 (608 – blanks/QC)

---

### 2. **growth_phenotype_summary_YPD2.csv.gz** (87 KB)
**Source:** Color phenotype measurements + strain metadata  
**Shape:** 608 rows × 43 columns

**Key Columns:**

| Column | Type | Example | Notes |
|--------|------|---------|-------|
| **Identifiers** | | | |
| `Strain ID` | Int | 1–335 | Unique strain identifier |
| `Strain Name` | String | "C_123" | Lab identifier |
| `Species` | String | "R. mucilaginosa" | 16 Rhodotorula species |
| `Library Plate` | Int | 1–4 | MS batch effect covariate |
| **Color Phenotypes (CIELab)** | | | |
| `Median_ColorLab_L*Mean` | Float | 40–70 | **Brightness** (0=black, 100=white) |
| `Median_ColorLab_a*Mean` | Float | −20–20 | **Red-Green axis** (−=green, +=red) |
| `Median_ColorLab_b*Mean` | Float | −10–30 | **Yellow-Blue axis** (−=blue, +=yellow) |
| **Color Phenotypes (HSV, redundant)** | | | |
| `Median_ColorHSV_HueMean` | Float | 0–360 | Hue angle |
| `Median_ColorHSV_SaturationMean` | Float | 0–100 | Color saturation |
| `Median_ColorHSV_BrightnessMean` | Float | 0–100 | Brightness (alternative) |
| **Image-Based Metrics** | | | |
| `Median_Shape_Area`, `MAD_*`, `Mean_*`, `SD_*` | Float | — | Cell morphology (not used in analysis) |
| **Metadata** | | | |
| `Incubation Temp (°C)` | Int | 25–30 | Growth condition |
| `Media` | String | "YPD" | Growth medium |

**Phenotypes Used:** `Median_ColorLab_L*Mean`, `Median_ColorLab_a*Mean`, `Median_ColorLab_b*Mean` (CIELab chosen over HSV to reduce redundancy)

---

## Phase 0: Batch Assessment & QC

**Script:** `scripts/phase0_batch_assessment.py`  
**Purpose:** Detect batch effects and confounder; decide analysis strategy  
**Output:** `analysis/phase0_*.csv.gz` + `analysis/phase0_decision.json`

### Outputs

#### phase0_decision.json
```json
{
  "plate_effect": true,              // Library Plate affects signal (F=95.16, p<1e-16)
  "species_effect": true,            // Species confounds correlations (F=32.65, p<1e-16)
  "strategy": "stratified_with_plate", // Use species-stratified approach + plate covariate
  "phenotypes_to_use": ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"],
  "n_samples": 590,                  // 608 total - blanks
  "n_features_raw": 16332
}
```

#### phase0_ms2_aligned.csv.gz
- **Shape:** 590 samples × 16,332 features  
- **Columns:** Feature annotations (m/z, RT, adduct, etc.) + peak areas per sample
- **Job:** Realign raw MS2 data to matched phenotype samples only

#### phase0_metadata_aligned.csv.gz
- **Shape:** 590 rows × 11 columns
- **Columns:** `filename`, color phenotypes (L*, a*, b*), Library Plate, Species, sample_id
- **Job:** Align phenotype metadata to same 590 samples as MS2 data

---

## Phase 1: Feature Filtering & Preprocessing

**Script:** `scripts/phase1_feature_filtering.py`  
**Purpose:** Remove noise, normalize, stabilize variance  
**Input:** `phase0_ms2_aligned.csv.gz` (590 × 16,332)  
**Output:** `results/phase1/phase1_*.csv.gz`

### Processing Steps (in order)

1. **Zero-median filter** → Remove features with no signal  
   - Removed: 8,991 features
2. **Prevalence filter** → Keep features in ≥10% of samples  
   - Removed: 0 features (all surviving step 1 passed)
3. **Outlier removal** → Cap values >3×IQR per feature  
   - Capped: 78,561 outlier values
4. **Variance filter** → Keep CV ≥ 0.1 in detected samples  
   - Removed: 0 features
5. **Total area normalization** → Scale each sample to sum = 1.0  
   - Equalizes sequencing depth / MS signal intensity
6. **Log2 transformation** → Stabilize variance (pseudocount = 1e-6)  
   - Makes abundance distributions approximately normal

### Outputs

#### phase1_features_filtered.csv.gz
- **Shape:** 590 × 7,341 (44.9% of raw features)
- **Data type:** Log2-normalized, float64
- **Rows:** Samples (same order as phase0_metadata_aligned)
- **Columns:** Feature indices (0–7340) with normalized abundances
- **Use:** Input to Phase 2 correlation analysis

#### phase1_feature_metadata.csv.gz
- **Shape:** 7,341 × 4
- **Columns:**
  - `feature_index`: Which feature in raw MS2 data (0–16331)
  - `prevalence`: Fraction of samples with detection (0.0–1.0)
  - `mean_log_intensity`: Average log2 abundance
  - `std_log_intensity`: Standard deviation across samples
- **Use:** Understand feature quality; QC check

#### phase1_phenotype_data.csv.gz
- **Shape:** 590 × 11 (same as phase0_metadata_aligned)
- **Columns:** Color phenotypes, covariates (Library Plate, Species)
- **Use:** Phenotype values for correlation computation

---

## Phase 2: Correlation Analysis (Main)

**Script:** `scripts/phase2_correlation_analysis.py`  
**Purpose:** Compute Spearman partial correlations + two-stage FDR correction  
**Input:** `phase1_features_filtered.csv.gz`, `phase1_phenotype_data.csv.gz`  
**Output:** `results/phase2/phase2_*.csv.gz`

### Statistical Method

**Correlation:** Spearman rank (robust to non-normality, outliers)  
**Partial correlation:** Controlling for Library Plate batch effect  
**N tests:** 7,341 features × 3 phenotypes = 22,023 tests  
**Multiple testing correction:** Two-stage Benjamini-Hochberg (FDR)

- **Stage 1 (Within-phenotype):** q < 0.05 within each phenotype (controls 7,341 comparisons)
- **Stage 2 (Across-phenotype):** q < 0.05 after re-correcting Stage 1 q-values (controls 3 comparisons)

### Outputs

#### phase2_all_correlations.csv.gz
- **Shape:** 22,023 rows × 10 columns
- **Rows:** One per feature-phenotype pair
- **Columns:**

| Column | Type | Range | Meaning |
|--------|------|-------|---------|
| `phenotype` | String | L*, a*, b* | Which color axis |
| `feature_index` | Int | 0–7340 | Feature ID |
| `rho` | Float | −1 to 1 | Spearman correlation coefficient |
| `pval` | Float | 1e-100 to 1.0 | Uncorrected p-value |
| `n_samples` | Int | 567 | Samples used (some missing phenotype data) |
| `reject_stage1` | Bool | T/F | Significant after Stage 1 (within-phenotype) |
| `q_value_stage1` | Float | 1e-100 to 1.0 | Within-phenotype q-value (FDR) |
| `q_value_global` | Float | 1e-100 to 1.0 | **Final q-value** (after across-phenotype correction) |
| `reject_global` | Bool | T/F | Significant at q < 0.05 (Stage 2) |
| `tier` | String | Tier1/2/3 | Confidence tier (effect size + q-value combo) |

#### phase2_tier1_hits.csv.gz
- **Shape:** 12,269 rows × 10 columns (Tier 1 only)
- **Selection:** `abs(rho) > 0.30 AND q_value_stage1 < 0.05`
- **Interpretation:** High-confidence, biologically meaningful signals
- **By phenotype:**
  - L* (brightness): 6,200 features (ρ > 0.30, explains 9%+ variance each)
  - a* (red-green): 5,846 features (|ρ| > 0.30)
  - b* (yellow-blue): 223 features (weaker signal)

#### phase2_tier12_hits.csv.gz
- **Shape:** 16,586 rows × 10 columns (Tier 1 + 2)
- **Selection:** `abs(rho) > 0.25 AND q_value_stage1 < 0.05`
- **Tier 2:** 4,317 additional features (|ρ| 0.25–0.30, exploratory)

#### phase2_summary.json
- **Contents:** Metadata about the correlation run
- **Fields:** phenotype counts, method, FDR thresholds

---

## Phase 3: Stratified Species Analysis

**Script:** `scripts/phase3_stratified_species_analysis.py`  
**Purpose:** Repeat Phase 2 separately per major species; identify species-specific signals  
**Input:** `phase1_features_filtered.csv.gz`, phenotype data  
**Output:** `results/phase3/phase3_*.csv.gz`

### Key Decision: Major Species
- Filter to species with ≥3 samples
- Result: 16 major species, 582 samples analyzed (8 excluded)

### Analysis Parts

#### Part 1: Between-Species Discriminant Features
**Job:** Find features that vary across species (regardless of phenotype)

##### Outputs
- **phase3_discriminant_features.csv.gz** (200 rows)
  - Top 200 features by variance across species
  - Columns: `feature_index`, `feature_column`, `variance_across_species`
  - Use: Screen for species effects masquerading as phenotype signals

#### Part 2: Within-Species Correlations
**Job:** Repeat Phase 2 correlation for each species separately

##### Outputs
- **phase3_within_species_correlations.csv.gz**
  - Shape: # species × # features × # phenotypes rows
  - Columns: `species`, `phenotype`, `feature`, `rho`, `pval`, `n_samples`, `q_value`, `reject`
  - Use: Check effect size consistency across species; identify species-specific hits
  
- **phase3_phenotype_by_species.csv.gz**
  - Shape: # species × # phenotypes rows
  - Columns: `species`, `phenotype`, `mean`, `std`, `n_samples`, `min`, `max`
  - Use: Understand phenotype distribution across species

- **phase3_summary.json**
  - Metadata: # species analyzed, sample counts, method

---

## Script Dependency Graph

```
input_data/
├── Rhodotorula_MS2_aligned_features_ms2.csv.gz
└── growth_phenotype_summary_YPD2.csv.gz
    ↓
[phase0_batch_assessment.py]
    ↓
results/phase0/
├── phase0_ms2_aligned.csv.gz
├── phase0_metadata_aligned.csv.gz
└── phase0_decision.json
    ↓
[phase1_feature_filtering.py]
    ↓
results/phase1/
├── phase1_features_filtered.csv.gz
├── phase1_feature_metadata.csv.gz
└── phase1_phenotype_data.csv.gz
    ├─→ [phase2_correlation_analysis.py]
    │       ↓
    │   results/phase2/
    │   ├── phase2_all_correlations.csv.gz
    │   ├── phase2_tier1_hits.csv.gz
    │   └── phase2_tier12_hits.csv.gz
    │
    └─→ [phase3_stratified_species_analysis.py]
            ↓
        results/phase3/
        ├── phase3_within_species_correlations.csv.gz
        ├── phase3_discriminant_features.csv.gz
        └── phase3_phenotype_by_species.csv.gz
```

---

## Querying the Data: Common Patterns

### Find the Strongest Features for a Phenotype
```python
import pandas as pd
import gzip

# Load Tier 1 hits
with gzip.open('results/phase2/phase2_tier1_hits.csv.gz', 'rt') as f:
    tier1 = pd.read_csv(f)

# Filter to one phenotype, sort by effect size
brightness = tier1[tier1['phenotype'] == 'Median_ColorLab_L*Mean']
top = brightness.sort_values('rho', ascending=False, key=abs).head(10)
print(top[['feature_index', 'rho', 'q_value_stage1']])
```

### Get Feature Metadata (m/z, RT, Adduct)
```python
# Load raw MS2 data
with gzip.open('input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz', 'rt') as f:
    ms2 = pd.read_csv(f)

# Look up feature 2755
feature_2755 = ms2.iloc[2755]
print(f"m/z: {feature_2755['row m/z']}")
print(f"RT: {feature_2755['row retention time']} min")
print(f"Adduct: {feature_2755['adduct']}")
print(f"MS2 available: {feature_2755['has_ms2']}")
```

### Get Feature Abundance for One Sample
```python
# Load filtered features
with gzip.open('results/phase1/phase1_features_filtered.csv.gz', 'rt') as f:
    features = pd.read_csv(f)

# Sample 0, all features
sample_0_abundances = features.iloc[0, :].values
print(f"Feature 2755 abundance (log2): {sample_0_abundances[2755]}")
```

### Check Species Distribution
```python
# Load phenotype data
with gzip.open('results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    pheno = pd.read_csv(f)

# How many samples per species?
print(pheno['Species'].value_counts())
```

### Compare Effect Sizes Across Species
```python
# Load Phase 3 results
with gzip.open('results/phase3/phase3_within_species_correlations.csv.gz', 'rt') as f:
    p3 = pd.read_csv(f)

# Filter to Feature 2755, L* phenotype
feat_2755 = p3[(p3['feature'] == 2755) & (p3['phenotype'] == 'Median_ColorLab_L*Mean')]
print(feat_2755[['species', 'rho', 'n_samples']])
```

---

## File Sizes & Compression

| File | Size (Gzipped) | Size (Uncompressed) | Notes |
|------|---------------|--------------------|-------|
| Rhodotorula_MS2_aligned_features_ms2.csv.gz | 18 MB | 42 MB | Feature matrix: 590 × 16,332 |
| growth_phenotype_summary_YPD2.csv.gz | 87 KB | 240 KB | Metadata: 608 × 43 |
| phase1_features_filtered.csv.gz | 29 MB | 81 MB | Core analysis data: 590 × 7,341 |
| phase1_feature_metadata.csv.gz | 170 KB | 338 KB | QC: 7,341 × 4 |
| phase2_all_correlations.csv.gz | 2.0 MB | 3.0 MB | Full results: 22,023 × 10 |
| phase2_tier1_hits.csv.gz | 1.7 MB | 2.3 MB | High-confidence: 12,269 × 10 |
| phase3_within_species_correlations.csv.gz | 600 KB | 1.2 MB | Species-stratified: ~16k rows × 8 |

---

## What Changed? Benjamini-Hochberg Correction

**Previous approach (Phase 3.0):** Bonferroni correction for across-phenotype tests  
**Current approach (Phase 3.1+):** Benjamini-Hochberg (FDR) correction for across-phenotype tests

**Impact:**
- Benjamini-Hochberg is less conservative (higher power) than Bonferroni
- Still controls false discovery rate at 5% level
- More appropriate for large-scale data (22k tests)
- Same Stage 1 (within-phenotype FDR), but Stage 2 uses FDR instead of FWER

---

## Reproducibility Checklist

- [ ] Run `scripts/phase0_batch_assessment.py` → `results/phase0/`
- [ ] Run `scripts/phase1_feature_filtering.py` → `results/phase1/`
- [ ] Run `scripts/phase2_correlation_analysis.py` → `results/phase2/`
- [ ] Run `scripts/phase3_stratified_species_analysis.py` → `results/phase3/`
- [ ] Verify output file counts and sizes match above
- [ ] Load `phase2_tier1_hits.csv.gz`, check Feature 2755 ρ=0.735

**Runtime:** ~20–30 minutes total (single core)  
**Dependencies:** pandas, numpy, scipy, statsmodels
