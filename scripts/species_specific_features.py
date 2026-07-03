#!/usr/bin/env python3
"""
SPECIES-SPECIFIC FEATURE ANALYSIS
Identify metabolite features present in only one species and absent/very low in others.

Approach:
  1. Build a clean sample-metadata table mapping every MS2 sample ID to its
     species (via Strain ID → growth_phenotype_summary), sample_type (C/SUP),
     strain_id, and library_plate.
  2. Using the raw Phase 0 peak-area matrix (16 332 features × 590 samples),
     compute per-species detection rates (fraction of samples with peak area > 0)
     for every feature.
  3. Score species-specificity for each feature-species pair and test
     significance with Fisher's exact test.
  4. Apply Benjamini-Hochberg FDR correction across all tests.
  5. Save results to results/species_specific/.

A feature is considered "species-specific" when:
  - detection rate in the target species ≥ 0.80
  - max detection rate in any other species ≤ 0.20
  - the target species has ≥ 3 samples
  - Fisher's exact test FDR q < 0.05
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import re
import json
import os
import warnings
warnings.filterwarnings('ignore')

# ── paths ────────────────────────────────────────────────────────────────────
script_dir  = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(script_dir)
input_dir   = os.path.join(project_dir, 'input_data')
results_dir = os.path.join(project_dir, 'results')
output_dir  = os.path.join(results_dir, 'species_specific')
os.makedirs(output_dir, exist_ok=True)

# ── parameters ───────────────────────────────────────────────────────────────
DETECTION_THRESHOLD     = 0          # peak area > 0  → detected
MIN_TARGET_PREVALENCE   = 0.80       # ≥ 80 % of target-species samples detect
MAX_OTHER_PREVALENCE    = 0.20       # ≤ 20 % of any other species detects
MIN_SPECIES_SAMPLES     = 3          # need ≥ 3 samples for reliability
FDR_ALPHA               = 0.05

# ============================================================================
# STEP 1 – Build clean sample metadata
# ============================================================================
print("=" * 80)
print("STEP 1: BUILD CLEAN SAMPLE METADATA")
print("=" * 80)

gps = pd.read_csv(os.path.join(input_dir, 'growth_phenotype_summary_YPD2.csv.gz'))
p0_ms2 = pd.read_csv(
    os.path.join(results_dir, 'phase0', 'phase0_ms2_aligned.csv.gz'),
    compression='gzip')

sample_cols = [c for c in p0_ms2.columns if 'Peak area' in c]
ms2_sids = [c.replace('.mzML Peak area', '').replace('.mzMLPeak area', '')
            for c in sample_cols]

def extract_strain_id(sid: str):
    m = re.search(r'(\d+)', sid)
    return int(m.group(1)) if m else None

def extract_sample_type(sid: str):
    return sid.split('_')[0] if '_' in sid else sid

sid_to_species   = dict(zip(gps['Strain ID'], gps['Species']))
sid_to_plate     = dict(zip(gps['Strain ID'], gps['Library Plate']))
sid_to_strain_nm = dict(zip(gps['Strain ID'], gps['Strain Name']))

rows = []
for col, sid in zip(sample_cols, ms2_sids):
    strain_id = extract_strain_id(sid)
    rows.append({
        'ms2_column':      col,
        'sample_id':       sid,
        'strain_id':       strain_id,
        'sample_type':     extract_sample_type(sid),
        'species':         sid_to_species.get(strain_id, np.nan),
        'library_plate':   sid_to_plate.get(strain_id, np.nan),
        'strain_name':     sid_to_strain_nm.get(strain_id, np.nan),
    })

clean_meta = pd.DataFrame(rows)
clean_meta.to_csv(os.path.join(output_dir, 'sample_metadata_clean.csv'),
                   index=False)

n_total     = len(clean_meta)
n_with_sp   = clean_meta['species'].notna().sum()
n_c         = (clean_meta['sample_type'] == 'C').sum()
n_sup       = (clean_meta['sample_type'] == 'SUP').sum()
species_dist = clean_meta['species'].value_counts(dropna=False)

print(f"  Total samples:           {n_total}")
print(f"  With species assignment: {n_with_sp}")
print(f"  Sample types:            C={n_c}, SUP={n_sup}")
print(f"  Species distribution:")
for sp, cnt in species_dist.items():
    print(f"    {sp}: {cnt}")
print(f"  ✓ Saved: sample_metadata_clean.csv")

# ============================================================================
# STEP 2 – Load raw peak-area matrix and align species
# ============================================================================
print("\n" + "=" * 80)
print("STEP 2: LOAD RAW PEAK-AREA MATRIX")
print("=" * 80)

peak_area = p0_ms2[sample_cols].T            # samples × features
peak_area.index = range(peak_area.shape[0])
n_features = peak_area.shape[1]
print(f"  Peak-area matrix: {peak_area.shape[0]} samples × {n_features} features")
print(f"  Zero fraction:    {(peak_area.values == 0).mean():.4f}")
print(f"  Min (nonzero):    {peak_area.values[peak_area.values > 0].min():.0f}")

# species vector aligned to peak_area rows
species_vec = clean_meta['species'].values
sample_type_vec = clean_meta['sample_type'].values

unique_species = [s for s in pd.unique(species_vec) if pd.notna(s)]
print(f"  Unique species with assignments: {len(unique_species)}")

# ============================================================================
# STEP 3 – Per-species detection rates
# ============================================================================
print("\n" + "=" * 80)
print("STEP 3: COMPUTE PER-SPECIES DETECTION RATES")
print("=" * 80)

# detection matrix: True where peak area > DETECTION_THRESHOLD
detected = (peak_area.values > DETECTION_THRESHOLD)   # samples × features (bool)

# For each species, compute detection rate per feature
detection_rates = {}    # species → array(n_features,) of detection rates
n_samples_by_sp = {}    # species → number of samples
for sp in unique_species:
    mask = species_vec == sp
    n_samples_by_sp[sp] = int(mask.sum())
    detection_rates[sp] = detected[mask].mean(axis=0)

print(f"  Computed detection rates for {len(unique_species)} species")
for sp in unique_species:
    n = n_samples_by_sp[sp]
    n_feat_high = int((detection_rates[sp] >= MIN_TARGET_PREVALENCE).sum())
    print(f"    {sp}: {n} samples, {n_feat_high} features ≥ {MIN_TARGET_PREVALENCE:.0%} prevalence")

# ============================================================================
# STEP 4 – Score species-specificity + Fisher's exact test
# ============================================================================
print("\n" + "=" * 80)
print("STEP 4: SPECIES-SPECIFICITY SCORING & STATISTICAL TESTS")
print("=" * 80)

results = []
for sp in unique_species:
    if n_samples_by_sp[sp] < MIN_SPECIES_SAMPLES:
        continue

    dr_target = detection_rates[sp]
    n_target  = n_samples_by_sp[sp]

    # other species detection rates
    other_sps = [s for s in unique_species if s != sp]
    dr_others_max  = np.zeros(n_features)
    dr_others_mean = np.zeros(n_features)
    n_others_total = 0
    for other in other_sps:
        dr_others_max  = np.maximum(dr_others_max, detection_rates[other])
        dr_others_mean += detection_rates[other] * n_samples_by_sp[other]
        n_others_total += n_samples_by_sp[other]
    if n_others_total > 0:
        dr_others_mean /= n_others_total

    specificity = dr_target - dr_others_max

    for feat_idx in range(n_features):
        dr_t = dr_target[feat_idx]
        dr_o = dr_others_max[feat_idx]
        spec = specificity[feat_idx]

        # only keep features with some signal
        if dr_t < MIN_TARGET_PREVALENCE:
            continue
        if dr_o > MAX_OTHER_PREVALENCE:
            continue

        # Fisher's exact test (2×2)
        #          | detected | not detected |
        # target   |  a       |  b           |
        # other    |  c       |  d           |
        a = int(detected[species_vec == sp, feat_idx].sum())
        b = n_target - a
        c = int(detected[species_vec != sp, feat_idx].sum())
        # only count other species samples with species assignment
        other_mask = np.array([s != sp and pd.notna(s)
                               for s in species_vec])
        d = int(other_mask.sum()) - c

        try:
            _, fisher_p = fisher_exact([[a, b], [c, d]], alternative='greater')
        except ValueError:
            fisher_p = 1.0

        results.append({
            'feature_index':           feat_idx,
            'target_species':          sp,
            'n_samples_target':        n_target,
            'detection_rate_target':   round(dr_t, 4),
            'detection_rate_others_max':  round(dr_o, 4),
            'detection_rate_others_mean': round(dr_others_mean[feat_idx], 4),
            'specificity_score':       round(spec, 4),
            'n_detected_target':       a,
            'n_detected_others':       c,
            'fisher_p':                fisher_p,
        })

print(f"  Total candidate feature-species pairs: {len(results)}")

if len(results) == 0:
    print("\n  ⚠ No features met the species-specificity criteria.")
    print("    Consider relaxing MIN_TARGET_PREVALENCE or MAX_OTHER_PREVALENCE.")
    # save empty results
    pd.DataFrame(results).to_csv(
        os.path.join(output_dir, 'species_specific_features.csv'),
        index=False)
    summary = {
        'n_features_total': n_features,
        'n_species_with_enough_samples': sum(1 for sp in unique_species
                                             if n_samples_by_sp[sp] >= MIN_SPECIES_SAMPLES),
        'n_species_specific': 0,
        'parameters': {
            'detection_threshold': DETECTION_THRESHOLD,
            'min_target_prevalence': MIN_TARGET_PREVALENCE,
            'max_other_prevalence': MAX_OTHER_PREVALENCE,
            'min_species_samples': MIN_SPECIES_SAMPLES,
            'fdr_alpha': FDR_ALPHA,
        }
    }
    with open(os.path.join(output_dir, 'species_specific_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print("\n✓ ANALYSIS COMPLETE (no species-specific features found)\n")
    import sys
    sys.exit(0)

# ============================================================================
# STEP 5 – Multiple testing correction (BH FDR)
# ============================================================================
print("\n" + "=" * 80)
print("STEP 5: MULTIPLE TESTING CORRECTION")
print("=" * 80)

results_df = pd.DataFrame(results)
reject, q_vals, _, _ = multipletests(results_df['fisher_p'],
                                      method='fdr_bh', alpha=FDR_ALPHA)
results_df['fisher_q'] = q_vals
results_df['is_species_specific'] = (
    (results_df['fisher_q'] < FDR_ALPHA) &
    (results_df['detection_rate_target'] >= MIN_TARGET_PREVALENCE) &
    (results_df['detection_rate_others_max'] <= MAX_OTHER_PREVALENCE)
)

n_specific = int(results_df['is_species_specific'].sum())
print(f"  Total tests:                {len(results_df)}")
print(f"  FDR significant (q < {FDR_ALPHA}):   {n_specific}")

# ============================================================================
# STEP 6 – Annotate with feature metadata (m/z, RT, adduct)
# ============================================================================
print("\n" + "=" * 80)
print("STEP 6: ANNOTATE WITH FEATURE METADATA")
print("=" * 80)

annotation_cols = ['row ID', 'row m/z', 'row retention time',
                   'row ion mobility', 'row CCS', 'adduct',
                   'parent_mass', 'has_ms2']
available_annot = [c for c in annotation_cols if c in p0_ms2.columns]
feature_annot = p0_ms2[available_annot].copy()
feature_annot['feature_index'] = range(len(feature_annot))

results_df = results_df.merge(feature_annot, on='feature_index', how='left')

# Sort by specificity score (descending), then Fisher q (ascending)
results_df = results_df.sort_values(
    ['specificity_score', 'fisher_q'], ascending=[False, True]
).reset_index(drop=True)

# ============================================================================
# STEP 7 – Save outputs
# ============================================================================
print("\n" + "=" * 80)
print("STEP 7: SAVE OUTPUTS")
print("=" * 80)

# All candidate results
results_df.to_csv(
    os.path.join(output_dir, 'species_specific_features.csv'),
    index=False)

# Filter to only species-specific (FDR significant)
specific_df = results_df[results_df['is_species_specific']].copy()
specific_df.to_csv(
    os.path.join(output_dir, 'species_specific_significant.csv'),
    index=False)

# Summary by species
summary_by_species = []
for sp in unique_species:
    if n_samples_by_sp[sp] < MIN_SPECIES_SAMPLES:
        continue
    sp_specific = specific_df[specific_df['target_species'] == sp]
    summary_by_species.append({
        'species': sp,
        'n_samples': n_samples_by_sp[sp],
        'n_species_specific_features': len(sp_specific),
    })
summary_by_species_df = pd.DataFrame(summary_by_species)

# Summary JSON
summary = {
    'n_features_total': int(n_features),
    'n_samples_total': int(n_total),
    'n_species_with_enough_samples': int(sum(1 for sp in unique_species
                                             if n_samples_by_sp[sp] >= MIN_SPECIES_SAMPLES)),
    'n_candidate_pairs': int(len(results_df)),
    'n_species_specific': int(n_specific),
    'species_specific_by_species': {
        row['species']: int(row['n_species_specific_features'])
        for _, row in summary_by_species_df.iterrows()
    },
    'parameters': {
        'detection_threshold': DETECTION_THRESHOLD,
        'min_target_prevalence': MIN_TARGET_PREVALENCE,
        'max_other_prevalence': MAX_OTHER_PREVALENCE,
        'min_species_samples': MIN_SPECIES_SAMPLES,
        'fdr_alpha': FDR_ALPHA,
        'multiple_testing_correction': 'Benjamini-Hochberg FDR',
        'statistical_test': "Fisher's exact test (one-sided, greater)",
    }
}

with open(os.path.join(output_dir, 'species_specific_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

# Print summary table
print("\n  Species-specific features by species:")
print(f"  {'Species':<40} {'N samples':<12} {'N features':<12}")
print(f"  {'-'*40} {'-'*12} {'-'*12}")
for _, row in summary_by_species_df.iterrows():
    print(f"  {row['species']:<40} {row['n_samples']:<12} {row['n_species_specific_features']:<12}")

print(f"\n  ✓ Saved: species_specific_features.csv ({len(results_df)} rows)")
print(f"  ✓ Saved: species_specific_significant.csv ({len(specific_df)} rows)")
print(f"  ✓ Saved: species_specific_summary.json")
print(f"  ✓ Saved: sample_metadata_clean.csv")

print("\n" + "=" * 80)
print("SPECIES-SPECIFIC FEATURE ANALYSIS COMPLETE")
print("=" * 80)
