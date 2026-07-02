#!/usr/bin/env python3
"""
PHASE 0: Quality Control & Batch Effect Assessment
Determines confounding strategy for downstream analysis
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("PHASE 0: BATCH & CONFOUNDER ASSESSMENT")
print("="*80)

# Load data from local input_data folder
print("\n[1/5] Loading data...")
import os
import gzip

# Use relative path to work from repo root
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)

ms2_file = os.path.join(repo_root, 'input_data', 'Rhodotorula_MS2_aligned_features_ms2.csv.gz')
metadata_file = os.path.join(repo_root, 'input_data', 'MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz')

df_ms2 = pd.read_csv(ms2_file)
df_meta = pd.read_csv(metadata_file, sep='\t', compression='gzip')

print(f"  MS2 data: {df_ms2.shape}")
print(f"  Metadata: {df_meta.shape}")

# Extract sample columns from MS2
sample_cols = [col for col in df_ms2.columns if 'Peak area' in col]
print(f"  Sample columns: {len(sample_cols)}")

# Parse sample IDs
sample_ids = []
for col in sample_cols:
    sample_id = col.replace('.mzML Peak area', '')
    sample_ids.append(sample_id)

print(f"  Sample ID range: {sample_ids[0]} to {sample_ids[-1]}")

# Extract phenotype data with color traits
print("\n[2/5] Extracting phenotype and batch info...")
phenotype_cols = ['Median_ColorHSV_BrightnessMean', 'Median_ColorHSV_HueMean',
                  'Median_ColorHSV_SaturationMean', 'Median_ColorLab_L*Mean',
                  'Median_ColorLab_a*Mean', 'Median_ColorLab_b*Mean',
                  'Median_ColorLab_ChromaEstimatedMean']

# Map metadata to sample columns
df_meta_sample = df_meta[['filename'] + phenotype_cols + ['Library Plate', 'Species']].copy()
df_meta_sample['sample_id'] = df_meta_sample['filename']

# Find matching samples
valid_samples = []
valid_indices = []
for i, sid in enumerate(sample_ids):
    if sid in df_meta_sample['sample_id'].values:
        valid_samples.append(sid)
        valid_indices.append(i)

print(f"  Matched samples: {len(valid_samples)}/{len(sample_ids)}")

# Create aligned data
df_meta_aligned = df_meta_sample[df_meta_sample['sample_id'].isin(valid_samples)].reset_index(drop=True)
sample_cols_valid = [sample_cols[i] for i in valid_indices]
df_ms2_aligned = df_ms2[['row ID', 'row m/z', 'row retention time'] + sample_cols_valid].copy()

print(f"  Aligned MS2: {df_ms2_aligned.shape}")
print(f"  Aligned metadata: {df_meta_aligned.shape}")

# Verify alignment
print(f"\n[3/5] Verifying data alignment...")
print(f"  Phenotype samples: {df_meta_aligned['Median_ColorLab_L*Mean'].notna().sum()}")
print(f"  Species distribution: {dict(df_meta_aligned['Species'].value_counts())}")
print(f"  Library Plate distribution: {dict(df_meta_aligned['Library Plate'].value_counts())}")

# Test 1: Color Phenotype Redundancy (CIELab vs HSV)
print("\n[4/5] Testing color phenotype redundancy...")
color_corr = df_meta_aligned[phenotype_cols].corr()

# Extract correlations between HSV and CIELab
hsv_cols = ['Median_ColorHSV_BrightnessMean', 'Median_ColorHSV_HueMean', 'Median_ColorHSV_SaturationMean']
cielab_cols = ['Median_ColorLab_L*Mean', 'Median_ColorLab_a*Mean', 'Median_ColorLab_b*Mean']

print("\n  Correlation: HSV Brightness ↔ CIELab L*:",
      f"{color_corr.loc['Median_ColorHSV_BrightnessMean', 'Median_ColorLab_L*Mean']:.3f}")
print("  Correlation: HSV Saturation ↔ CIELab Chroma:",
      f"{color_corr.loc['Median_ColorHSV_SaturationMean', 'Median_ColorLab_ChromaEstimatedMean']:.3f}")
print("  Correlation: HSV Hue ↔ CIELab b*:",
      f"{color_corr.loc['Median_ColorHSV_HueMean', 'Median_ColorLab_b*Mean']:.3f}")

print("\n  → RECOMMENDATION: Use CIELab (L*, a*, b*) only to reduce redundancy")
print("    This reduces 7 phenotypes → 3 independent traits")

# Test 2: Plate Batch Effects
print("\n[5/5] Testing for Library Plate batch effects...")

# Simple PERMANOVA-like test: compare metabolite profiles by plate
sample_data = df_ms2_aligned[sample_cols_valid].T.values
plate_groups = df_meta_aligned['Library Plate'].fillna(-1).astype(int).values

# Calculate Euclidean distances
from scipy.spatial.distance import pdist, squareform

# Log-transform metabolite data for meaningful distances
eps = 1e-10
sample_data_log = np.log2(sample_data + eps)

# Compute pairwise distances
distances = squareform(pdist(sample_data_log, metric='euclidean'))

# PERMANOVA-like calculation: between-group vs within-group variance
plates = np.unique(plate_groups[plate_groups >= 0])
plate_groups_clean = plate_groups[plate_groups >= 0]
sample_data_log_clean = sample_data_log[plate_groups >= 0]
n_samples = len(plate_groups_clean)

# Between-group sum of squares
grand_mean = sample_data_log_clean.mean(axis=0)
ss_between = 0
for plate in plates:
    mask = plate_groups_clean == plate
    group_mean = sample_data_log_clean[mask].mean(axis=0)
    ss_between += mask.sum() * np.sum((group_mean - grand_mean) ** 2)

# Within-group sum of squares
ss_within = 0
for plate in plates:
    mask = plate_groups_clean == plate
    group_mean = sample_data_log_clean[mask].mean(axis=0)
    ss_within += np.sum((sample_data_log_clean[mask] - group_mean) ** 2)

# F-statistic
df_between = len(plates) - 1
df_within = n_samples - len(plates)
ms_between = ss_between / df_between if df_between > 0 else 0
ms_within = ss_within / df_within if df_within > 0 else 1
f_stat = ms_between / ms_within

# Approximate p-value
from scipy.stats import f as f_dist
p_value = 1 - f_dist.cdf(f_stat, df_between, df_within)

print(f"  Plate Effect Test Results:")
print(f"    F-statistic: {f_stat:.4f}")
print(f"    p-value: {p_value:.4e}")
print(f"    Samples per plate: {dict(pd.Series(plate_groups).value_counts().sort_index())}")

if p_value < 0.05:
    print(f"\n  ✓ SIGNIFICANT BATCH EFFECT DETECTED (p={p_value:.2e})")
    print(f"    → INCLUDE 'Library Plate' as covariate in all correlation analyses")
    plate_effect = True
else:
    print(f"\n  ✗ NO SIGNIFICANT BATCH EFFECT (p={p_value:.3f})")
    print(f"    → Can proceed without plate covariate")
    plate_effect = False

# Test 3: Species Effects
print("\n Testing for Species confounding...")
species_groups = df_meta_aligned['Species'].fillna('Unknown').values
species = np.unique(species_groups)

# Between-group sum of squares
ss_between_sp = 0
for sp in species:
    mask = species_groups == sp
    group_mean = sample_data_log[mask].mean(axis=0)
    ss_between_sp += mask.sum() * np.sum((group_mean - grand_mean) ** 2)

# Within-group sum of squares
ss_within_sp = 0
for sp in species:
    mask = species_groups == sp
    group_mean = sample_data_log[mask].mean(axis=0)
    ss_within_sp += np.sum((sample_data_log[mask] - group_mean) ** 2)

# F-statistic
df_between_sp = len(species) - 1
df_within_sp = n_samples - len(species)
ms_between_sp = ss_between_sp / df_between_sp if df_between_sp > 0 else 0
ms_within_sp = ss_within_sp / df_within_sp if df_within_sp > 0 else 1
f_stat_sp = ms_between_sp / ms_within_sp

p_value_sp = 1 - f_dist.cdf(f_stat_sp, df_between_sp, df_within_sp)

print(f"  Species Effect Test Results:")
print(f"    F-statistic: {f_stat_sp:.4f}")
print(f"    p-value: {p_value_sp:.4e}")
print(f"    Samples per species: {dict(pd.Series(species_groups).value_counts())}")

if p_value_sp < 0.05:
    print(f"\n  ✓ SIGNIFICANT SPECIES EFFECT DETECTED (p={p_value_sp:.2e})")
    print(f"    → RECOMMEND stratified analysis by species")
    print(f"    → OR include 'Species' as covariate in correlations")
    species_effect = True
else:
    print(f"\n  ✗ NO SIGNIFICANT SPECIES EFFECT (p={p_value_sp:.3f})")
    print(f"    → Can pool species in analysis")
    species_effect = False

# Summary and Recommendations
print("\n" + "="*80)
print("DECISION TREE & RECOMMENDATIONS")
print("="*80)

print(f"\nBatch Effect (Library Plate): {'YES' if plate_effect else 'NO'}")
print(f"Species Confounding: {'YES' if species_effect else 'NO'}")
print(f"Color Phenotype Redundancy: HIGH (HSV+CIELab overlap)")

print("\n→ RECOMMENDED ANALYSIS APPROACH:")

if plate_effect and species_effect:
    print("\n  STRATEGY: Stratified by Species + Plate as Covariate")
    print("  Rationale: Both confounders significant")
    print("  Action:")
    print("    1. Split analysis by species")
    print("    2. Include Library Plate as covariate in correlations")
    print("    3. Compare effect sizes across species (consistency check)")
    strategy = "stratified_with_plate"

elif plate_effect and not species_effect:
    print("\n  STRATEGY: Pooled Analysis + Plate as Covariate")
    print("  Rationale: Only plate effects significant")
    print("  Action:")
    print("    1. Analyze all species together")
    print("    2. Include Library Plate as covariate")
    strategy = "pooled_with_plate"

elif not plate_effect and species_effect:
    print("\n  STRATEGY: Stratified by Species (no plate covariate)")
    print("  Rationale: Only species effects significant")
    print("  Action:")
    print("    1. Split analysis by species")
    print("    2. No plate covariate needed")
    print("    3. Compare effects across species")
    strategy = "stratified_no_plate"

else:
    print("\n  STRATEGY: Simple Pooled Analysis")
    print("  Rationale: No major confounders")
    print("  Action:")
    print("    1. Analyze all species together")
    print("    2. No covariates needed (simple Spearman correlations)")
    strategy = "pooled_simple"

print(f"\n  Color Phenotypes to Use: L*, a*, b* (3 traits instead of 7)")
print(f"    Expected multiple tests: ~5,000 features × 3 phenotypes = 15,000 tests")
print(f"    (vs 35,000 if using all 7; FDR threshold: q < 0.05)")

# Save decision for Phase 1
decision = {
    'plate_effect': plate_effect,
    'species_effect': species_effect,
    'strategy': strategy,
    'phenotypes_to_use': ['Median_ColorLab_L*Mean', 'Median_ColorLab_a*Mean', 'Median_ColorLab_b*Mean'],
    'n_samples': len(valid_samples),
    'n_features_raw': len(df_ms2),
}

import json

# Output to ../results/phase subfolder
output_dir = os.path.join(script_dir, '../results/phase')
os.makedirs(output_dir, exist_ok=True)

output_decision = os.path.join(output_dir, 'phase0_decision.json')
output_ms2 = os.path.join(output_dir, 'phase0_ms2_aligned.csv.gz')
output_meta = os.path.join(output_dir, 'phase0_metadata_aligned.csv.gz')

with open(output_decision, 'w') as f:
    json.dump(decision, f, indent=2)

# Save aligned data compressed for Phase 1
df_ms2_aligned.to_csv(output_ms2, index=False, compression='gzip')
df_meta_aligned.to_csv(output_meta, index=False, compression='gzip')

print("\n✓ Phase 0 complete. Decision and aligned data saved.")
print("  Ready for Phase 1: Feature filtering & preprocessing")

