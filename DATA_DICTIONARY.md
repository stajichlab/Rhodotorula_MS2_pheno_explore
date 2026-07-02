# Data Dictionary

This document defines all columns in output CSV files across the analysis pipeline.

---

## Phase 0 Outputs

### `phase0_metadata_aligned.csv.gz`
**Purpose**: Sample metadata with species assignments and phenotype data  
**Rows**: 590 (one per sample)  
**Columns**: 11

| Column Name | Type | Description | Example |
|------------|------|-------------|---------|
| `filename` | string | Original sample filename | `Rhodo_S001.raw` |
| `Median_ColorHSV_BrightnessMean` | float | HSV brightness (0-100) | 52.3 |
| `Median_ColorHSV_HueMean` | float | HSV hue (0-360°) | 12.5 |
| `Median_ColorHSV_SaturationMean` | float | HSV saturation (0-100) | 68.2 |
| `Median_ColorLab_L*Mean` | float | **USED**: CIE Lab lightness (0-100) | 74.5 |
| `Median_ColorLab_a*Mean` | float | **USED**: CIE Lab red-green axis (-128 to +127) | 5.2 |
| `Median_ColorLab_b*Mean` | float | **USED**: CIE Lab yellow-blue axis (-128 to +127) | 8.1 |
| `Median_ColorLab_ChromaEstimatedMean` | float | Chroma (color saturation, Lab space) | 9.8 |
| `Library Plate` | string | Plate ID for technical covariate | `LP001` |
| `Species` | string | Taxonomic species assignment | `Rhodotorula mucilaginosa` |
| `sample_id` | string | Unique sample identifier | `Rhodo_001` |

**Notes**: 
- **Used phenotypes**: `Median_ColorLab_L*Mean`, `Median_ColorLab_a*Mean`, `Median_ColorLab_b*Mean`
- HSV columns not analyzed (redundant with Lab)
- Library Plate used as technical covariate

---

### `phase0_ms2_aligned.csv.gz`
**Purpose**: Aligned raw MS2 feature intensities (before filtering)  
**Rows**: 590 (samples)  
**Columns**: 7,341 (feature IDs)

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| Feature ID (e.g., `1`, `2`, `4`, ..., `16331`) | float | Log-transformed MS2 feature intensity | -3.14 |

**Notes**:
- Values are log-transformed (negative values expected, represent low abundance)
- Feature ID is position/index, not a compound name
- 590 × 7341 matrix (590 samples × 7341 original features)

---

## Phase 1 Outputs

### `phase1_phenotype_data.csv.gz`
**Purpose**: Cleaned phenotype data (Lab color space only)  
**Rows**: 590 (samples)  
**Columns**: 3

| Column Name | Type | Range | Description |
|------------|------|-------|-------------|
| `Median_ColorLab_L*Mean` | float | 61.2–78.8 | CIE Lab lightness (brightness) |
| `Median_ColorLab_a*Mean` | float | -0.48–14.4 | Red-green axis (positive=red) |
| `Median_ColorLab_b*Mean` | float | -1.10–9.08 | Yellow-blue axis (positive=yellow) |

**Notes**: Used for all phenotype-feature correlations

---

### `phase1_features_filtered.csv.gz`
**Purpose**: Filtered MS2 features (after variance/non-zero filtering)  
**Rows**: 590 (samples)  
**Columns**: 7,341 (filtered features)

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| Feature ID (e.g., `1`, `2`, `4`, ...) | float | Log-transformed intensity | -4.52 |

**Notes**:
- Same as Phase 0 but with features meeting quality thresholds
- Log-scale values (negative = low abundance)
- 7,341 features retained from 16,331 original (45%)

---

### `phase1_feature_metadata.csv.gz`
**Purpose**: Summary statistics for each filtered feature  
**Rows**: 7,341 (features)  
**Columns**: 4

| Column Name | Type | Description |
|------------|------|-------------|
| `feature_id` | int | Feature index/ID |
| `mean_abundance` | float | Mean log-intensity across samples |
| `std_abundance` | float | Standard deviation of log-intensity |
| `n_nonzero` | int | Count of samples with non-zero value |

---

## Phase 2 Outputs

### `phase2_all_correlations.csv.gz`
**Purpose**: All Spearman correlations between features and phenotypes  
**Rows**: 22,023 (feature-phenotype pairs with non-zero variance)  
**Columns**: 9

| Column Name | Type | Description |
|------------|------|-------------|
| `phenotype` | string | `Median_ColorLab_L*Mean`, `Median_ColorLab_a*Mean`, or `Median_ColorLab_b*Mean` |
| `feature_index` | int | Feature ID (column index in feature matrix) |
| `rho` | float | Spearman correlation coefficient (-1 to +1) |
| `pval` | float | Raw p-value from correlation test |
| `reject_stage1` | bool | True if significant after within-phenotype FDR |
| `q_value_stage1` | float | FDR-adjusted p-value (Benjamini-Hochberg, within-phenotype) |
| `q_value_global` | float | FDR-adjusted p-value (Benjamini-Hochberg, across phenotypes) |
| `reject_global` | bool | True if `q_value_global < 0.05` |
| `tier` | string | `Tier1_High` (ρ>0.30, q<0.05), `Tier2_Medium` (ρ>0.25), `Tier3_Exploratory` (ρ>0.20), or `Not_Significant` |

**Notes**:
- Two-stage FDR: Stage 1 within each phenotype, Stage 2 across phenotypes
- Correlation method: Spearman partial correlation (controlling for Library Plate)
- Sample size varies; typically n~590 (complete cases only)

---

### `phase2_tier1_hits.csv.gz`
**Purpose**: High-confidence correlations (effect size & significance)  
**Rows**: 6,847 (filtered from all_correlations)  
**Columns**: Same as `phase2_all_correlations.csv.gz`

**Filter criteria**: `abs(rho) > 0.30 AND q_value_stage1 < 0.05`

---

### `phase2_tier12_hits.csv.gz`
**Purpose**: High + medium confidence correlations  
**Rows**: 10,127  
**Columns**: Same as above

**Filter criteria**: `abs(rho) > 0.25 AND q_value_stage1 < 0.05`

---

### `phase2_summary.json`
**Purpose**: Summary statistics for Phase 2  
**Format**: JSON

```json
{
  "n_features_analyzed": 7341,
  "n_phenotypes": 3,
  "n_total_tests": 22023,
  "tier1_count": 6847,
  "tier2_count": 3280,
  "tier3_count": ...,
  "strategy": "pooled",
  "covariates": ["Library Plate"],
  "phenotypes_used": ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
}
```

---

## Phase 3 Outputs

### `phase3_discriminant_features.csv.gz`
**Purpose**: Features with high variance across species  
**Rows**: 200 (top discriminant features)  
**Columns**: 3

| Column Name | Type | Description |
|------------|------|-------------|
| `feature_index` | int | Feature ID |
| `feature_column` | string | Feature column name (same as feature_index) |
| `variance_across_species` | float | Variance of mean abundance across 10 species |

**Ranking**: Sorted descending by variance (top discriminant first)

---

### `phase3_within_species_correlations.csv.gz`
**Purpose**: Species-stratified metabolite-phenotype correlations  
**Rows**: 5,795 (feature-phenotype pairs tested within each species)  
**Columns**: 7

| Column Name | Type | Description |
|------------|------|-------------|
| `species` | string | Taxonomic name (e.g., `Rhodotorula mucilaginosa`) |
| `phenotype` | string | Color phenotype (`Median_ColorLab_L*Mean`, `a*Mean`, `b*Mean`) |
| `feature` | string | Feature ID |
| `rho` | float | Spearman correlation within species |
| `pval` | float | Raw p-value |
| `q_value` | float | FDR-adjusted p-value (Benjamini-Hochberg, within-species) |
| `reject` | bool | True if `q_value < 0.05` |

**Notes**:
- Only tests **discriminant features** (200 features × 3 phenotypes per species)
- Only includes species with n≥3 samples
- Library Plate included as covariate where applicable

---

### `phase3_phenotype_by_species.csv.gz`
**Purpose**: Phenotype statistics stratified by species  
**Rows**: 30 (10 species × 3 phenotypes)  
**Columns**: 7

| Column Name | Type | Description |
|------------|------|-------------|
| `species` | string | Species name |
| `phenotype` | string | Phenotype name |
| `mean` | float | Mean phenotype value across samples in species |
| `std` | float | Standard deviation |
| `n_samples` | int | Count of samples for this species |
| `min` | float | Minimum phenotype value |
| `max` | float | Maximum phenotype value |

---

### `phase3_summary.json`
**Purpose**: Summary statistics for Phase 3

```json
{
  "n_major_species": 10,
  "major_species": ["Rhodotorula mucilaginosa", "R. paludigena", ...],
  "n_samples_major": 259,
  "n_discriminant_features": 200,
  "n_within_species_correlations": 5795,
  "n_significant_within_species": 126
}
```

---

## Column Naming Conventions

- **Phenotypes**: `Median_ColorLab_*Mean` (standardized naming from input data)
- **Features**: Numeric strings (`1`, `2`, `4`, ..., `16331`)
- **Species**: Full taxonomic names (e.g., `Rhodotorula mucilaginosa`)
- **Correlation metrics**: `rho`, `pval`, `q_value` (Spearman notation)

---

## Special Values

- **Missing data**: Represented as `NaN` (handled by filtering before analysis)
- **Log scale**: Feature values are log-transformed (typically -20 to -3)
- **Correlation bounds**: `rho` always in range [-1, +1]
- **P-values**: Always in range (0, 1]
- **Q-values**: Always in range [0, 1]

---

## Reference

For full analysis methodology, see `docs/ANALYSIS_REPORT.md`

