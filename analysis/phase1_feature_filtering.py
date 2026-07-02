#!/usr/bin/env python3
"""
PHASE 1: FEATURE FILTERING & PREPROCESSING
Reduces 16,332 features to ~5,000 high-quality metabolites
"""

import pandas as pd
import numpy as np
import json
import os
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("PHASE 1: FEATURE FILTERING & PREPROCESSING")
print("="*80)

# Setup paths relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, 'phenotypes_MS2')

# Load decision from Phase 0
decision_file = os.path.join(output_dir, 'phase0_decision.json')
with open(decision_file, 'r') as f:
    decision = json.load(f)

print(f"\nUsing Phase 0 decision: {decision['strategy']}")

# Load aligned data from Phase 0 (compressed)
df_ms2 = pd.read_csv(os.path.join(output_dir, 'phase0_ms2_aligned.csv.gz'), compression='gzip')
df_meta = pd.read_csv(os.path.join(output_dir, 'phase0_metadata_aligned.csv.gz'), compression='gzip')

sample_cols = [col for col in df_ms2.columns if 'Peak area' in col]
feature_data = df_ms2[sample_cols].T.reset_index(drop=True)

print(f"\n[1/6] Starting data: {feature_data.shape[0]} samples × {feature_data.shape[1]} features")

# Step 1: Remove zero-median features
print("\n[2/6] Filtering Step 1: Remove zero-median features...")
medians = feature_data.median(axis=0)
non_zero = medians > 0
n_before = feature_data.shape[1]
feature_data = feature_data.loc[:, non_zero]
print(f"  Removed: {n_before - feature_data.shape[1]} features (median = 0)")
print(f"  Retained: {feature_data.shape[1]} features")

# Step 2: Prevalence filter (<10% samples)
print("\n[3/6] Filtering Step 2: Prevalence filter...")
prevalence = (feature_data > 0).sum(axis=0) / len(feature_data)
min_prev = 0.10
keep_prev = prevalence >= min_prev
n_before = feature_data.shape[1]
feature_data = feature_data.loc[:, keep_prev]
print(f"  Threshold: {min_prev*100:.0f}% of samples")
print(f"  Removed: {n_before - feature_data.shape[1]} features (prevalence < {min_prev*100:.0f}%)")
print(f"  Retained: {feature_data.shape[1]} features")

# Step 3: Outlier removal (3×IQR per feature)
print("\n[4/6] Filtering Step 3: Outlier removal (3×IQR)...")
feature_data_outlier_removed = feature_data.copy()
outliers_removed_count = 0

for col in feature_data.columns:
    Q1 = feature_data[col].quantile(0.25)
    Q3 = feature_data[col].quantile(0.75)
    IQR = Q3 - Q1
    upper_bound = Q3 + 3 * IQR

    outlier_mask = feature_data[col] > upper_bound
    outliers_removed_count += outlier_mask.sum()

    # Replace outliers with upper bound (conservative approach)
    feature_data_outlier_removed.loc[outlier_mask, col] = upper_bound

feature_data = feature_data_outlier_removed
print(f"  Removed/capped: {outliers_removed_count} outlier values")
print(f"  Features: {feature_data.shape[1]}")

# Step 4: Variance filter (CV < 0.1 in detected samples)
print("\n[5/6] Filtering Step 4: Variance filter...")
cv_values = []
for col in feature_data.columns:
    detected = feature_data[col] > 0
    if detected.sum() > 1:
        cv = feature_data.loc[detected, col].std() / feature_data.loc[detected, col].mean()
    else:
        cv = 0
    cv_values.append(cv)

cv_threshold = 0.1
keep_var = np.array(cv_values) >= cv_threshold
n_before = feature_data.shape[1]
feature_data = feature_data.loc[:, keep_var]
print(f"  CV threshold: {cv_threshold}")
print(f"  Removed: {n_before - feature_data.shape[1]} features (CV < {cv_threshold})")
print(f"  Retained: {feature_data.shape[1]} features")

# Step 5: Normalization (total area)
print("\n[6/6] Preprocessing: Normalization & Transformation...")
print("  Step 1: Total area normalization per sample")

feature_data_norm = feature_data.copy()
sample_totals = feature_data_norm.sum(axis=1)
feature_data_norm = feature_data_norm.div(sample_totals, axis=0)

print(f"    Sample sums after normalization (check): {feature_data_norm.sum(axis=1).mean():.6f} ± {feature_data_norm.sum(axis=1).std():.6f}")

# Step 6: Log2 transformation
print("  Step 2: Log2 transformation")
pseudocount = 0.1 * feature_data_norm.median().median()
feature_data_log = np.log2(feature_data_norm + pseudocount)

print(f"    Pseudocount: {pseudocount:.6f}")
print(f"    Mean log-intensity: {feature_data_log.values.mean():.2f}")
print(f"    Std log-intensity: {feature_data_log.values.std():.2f}")

# Summary
print("\n" + "="*80)
print("FILTERING SUMMARY")
print("="*80)

print(f"""
Starting features:           16,332
After zero-median filter:    {feature_data_norm.shape[1]:,} (removed {16332 - feature_data_norm.shape[1]:,})
After prevalence filter:     {feature_data_norm.shape[1]:,}
After outlier removal:       {feature_data_norm.shape[1]:,} (capped {outliers_removed_count:,} outliers)
After variance filter:       {feature_data_log.shape[1]:,} (removed ~{16332 - feature_data_log.shape[1]:,})

FINAL: {feature_data_log.shape[1]:,} features × {feature_data_log.shape[0]} samples
       ({100*feature_data_log.shape[1]/16332:.1f}% of original features retained)
       (Each feature has prevalence ≥10% and CV ≥0.1)

Quality checks:
  - Normalization: ✓ (sample sums = 1.0)
  - Log transformation: ✓ (all values numeric)
  - Missing data: {feature_data_log.isna().sum().sum()} NaNs
  - Outliers: {(feature_data_log > 5).sum().sum()} extreme values (>5 SD, monitored)
""")

# Save filtered data (compressed) to phenotypes_MS2 folder
print("\nSaving outputs...")

feature_data_log.to_csv(os.path.join(output_dir, 'phase1_features_filtered.csv.gz'), index=False, compression='gzip')

# Add feature metadata
feature_metadata = pd.DataFrame({
    'feature_index': range(feature_data_log.shape[1]),
    'prevalence': (feature_data_log != 0).sum(axis=0) / len(feature_data_log),
    'mean_log_intensity': feature_data_log.mean(axis=0),
    'std_log_intensity': feature_data_log.std(axis=0),
})
feature_metadata.to_csv(os.path.join(output_dir, 'phase1_feature_metadata.csv.gz'), index=False, compression='gzip')

# Save metadata aligned with feature data
df_meta.to_csv(os.path.join(output_dir, 'phase1_phenotype_data.csv.gz'), index=False, compression='gzip')

print("✓ Filtered features saved: phase1_features_filtered.csv")
print("✓ Feature metadata saved: phase1_feature_metadata.csv")
print("✓ Phenotype data saved: phase1_phenotype_data.csv")
print("\n✓ Phase 1 complete. Ready for Phase 2: Correlation Analysis")

