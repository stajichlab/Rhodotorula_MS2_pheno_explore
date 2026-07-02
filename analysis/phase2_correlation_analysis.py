#!/usr/bin/env python3
"""
PHASE 2: CORRELATION ANALYSIS
Spearman partial correlations with two-stage FDR correction
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, rankdata
import json
import os
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("PHASE 2: CORRELATION ANALYSIS (SPEARMAN PARTIAL CORRELATIONS)")
print("="*80)

# Setup paths relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Load Phase 0 decision and Phase 1 data (compressed)
with open(os.path.join(script_dir, 'phase0_decision.json'), 'r') as f:
    decision = json.load(f)

features = pd.read_csv(os.path.join(script_dir, 'phase1_features_filtered.csv.gz'), compression='gzip')
metadata = pd.read_csv(os.path.join(script_dir, 'phase1_phenotype_data.csv.gz'), compression='gzip')

print(f"\nFeatures: {features.shape}")
print(f"Samples: {metadata.shape}")
print(f"Strategy: {decision['strategy']}")

# Select color phenotypes (CIELab only, not redundant HSV)
phenotype_cols = decision['phenotypes_to_use']
phenotypes = metadata[phenotype_cols].copy()

print(f"\nPhenotypes to analyze: {len(phenotype_cols)}")
for i, col in enumerate(phenotype_cols, 1):
    print(f"  {i}. {col}")

# Helper function: Spearman partial correlation
def spearman_partial_corr(X, Y, covariates=None):
    """
    Compute Spearman partial correlation between X and Y,
    controlling for covariates.
    """
    # Rank all variables (Spearman)
    X_ranked = rankdata(X)
    Y_ranked = rankdata(Y)

    if covariates is None or len(covariates) == 0:
        # Simple Spearman
        rho, pval = spearmanr(X_ranked, Y_ranked)
        return rho, pval

    # Partial correlation: regress out covariates
    from scipy.stats import linregress

    # Residuals after regressing out covariates
    X_residual = X_ranked.copy()
    Y_residual = Y_ranked.copy()

    for cov in covariates:
        # Regress X on covariate
        slope_x, intercept_x = linregress(cov, X_ranked)[:2]
        X_residual = X_ranked - (intercept_x + slope_x * cov)

        # Regress Y on covariate
        slope_y, intercept_y = linregress(cov, Y_ranked)[:2]
        Y_residual = Y_ranked - (intercept_y + slope_y * cov)

    # Correlation of residuals
    rho, pval = spearmanr(X_residual, Y_residual)
    return rho, pval

# Prepare covariates
print(f"\n[1/4] Preparing covariates...")

covariates = []
covariate_names = []

# Library Plate covariate (if needed)
if decision['plate_effect']:
    plate_numeric = pd.factorize(metadata['Library Plate'])[0].astype(float)
    covariates.append(plate_numeric)
    covariate_names.append('Library Plate')
    print(f"  ✓ Including Library Plate as covariate")

# Species covariate (if pooled analysis)
if 'pooled' in decision['strategy']:
    species_numeric = pd.factorize(metadata['Species'])[0].astype(float)
    covariates.append(species_numeric)
    covariate_names.append('Species')
    print(f"  ✓ Including Species as covariate")
else:
    print(f"  ✗ Species-stratified analysis (no species covariate)")

print(f"  Total covariates: {len(covariates)}")

# Compute correlations
print(f"\n[2/4] Computing Spearman correlations...")
print(f"  Total tests: {features.shape[1]} features × {len(phenotype_cols)} phenotypes = {features.shape[1] * len(phenotype_cols):,}")

results = []
for pheno_idx, phenotype in enumerate(phenotype_cols):
    print(f"\n  Processing phenotype {pheno_idx+1}/{len(phenotype_cols)}: {phenotype}")

    y = metadata[phenotype].values

    # Filter to samples without missing data
    valid = ~np.isnan(y)
    y_valid = y[valid]
    features_valid = features.iloc[valid].values

    X_covariates = [cov[valid] for cov in covariates]

    for feat_idx in range(features.shape[1]):
        x = features_valid[:, feat_idx]

        # Skip if all zeros/constant
        if np.std(x) == 0:
            continue

        # Compute Spearman partial correlation
        try:
            rho, pval = spearman_partial_corr(x, y_valid, X_covariates)

            results.append({
                'phenotype': phenotype,
                'feature_index': feat_idx,
                'rho': rho,
                'pval': pval,
                'n_samples': len(x)
            })
        except Exception as e:
            continue

    if (pheno_idx + 1) % 1 == 0:
        print(f"    Completed: {len(results)} correlations so far")

print(f"\n✓ Total correlations computed: {len(results):,}")

# Convert to DataFrame
df_results = pd.DataFrame(results)

# Two-stage FDR correction
print(f"\n[3/4] Multiple testing correction (Two-stage FDR)...")

# Stage 1: Within-phenotype FDR correction
stage1_results = []
for phenotype in phenotype_cols:
    df_pheno = df_results[df_results['phenotype'] == phenotype].copy()

    # Apply FDR correction within phenotype
    reject, q_values, _, _ = multipletests(df_pheno['pval'], method='fdr_bh', alpha=0.05)

    df_pheno['reject_stage1'] = reject
    df_pheno['q_value_stage1'] = q_values

    stage1_results.append(df_pheno)
    n_sig = reject.sum()
    print(f"  {phenotype}: {n_sig} features at q<0.05")

df_results = pd.concat(stage1_results, ignore_index=True)

# Stage 2: Across-phenotype Bonferroni (multiply threshold by number of phenotypes)
print(f"\n  Applying across-phenotype Bonferroni correction (×{len(phenotype_cols)})")
df_results['q_value_global'] = np.minimum(df_results['q_value_stage1'] * len(phenotype_cols), 1.0)
df_results['reject_global'] = df_results['q_value_global'] < 0.05

# Tier hits by effect size
print(f"\n[4/4] Tiering hits by confidence...")

def assign_tier(row):
    """Assign confidence tier based on effect size and significance"""
    rho_abs = abs(row['rho'])

    if rho_abs > 0.30 and row['q_value_stage1'] < 0.05:
        return 'Tier1_High'
    elif rho_abs > 0.25 and row['q_value_stage1'] < 0.05:
        return 'Tier2_Medium'
    elif rho_abs > 0.20 and row['q_value_stage1'] < 0.10:
        return 'Tier3_Exploratory'
    else:
        return 'Not_Significant'

df_results['tier'] = df_results.apply(assign_tier, axis=1)

# Summary statistics
print("\nTiered Results Summary:")
tier_counts = df_results['tier'].value_counts()
for tier in ['Tier1_High', 'Tier2_Medium', 'Tier3_Exploratory', 'Not_Significant']:
    count = tier_counts.get(tier, 0)
    print(f"  {tier:20} {count:6,} features")

# Top hits per phenotype
print("\nTop 10 hits per phenotype:")
for phenotype in phenotype_cols:
    df_pheno = df_results[df_results['phenotype'] == phenotype].sort_values('q_value_stage1')
    print(f"\n  {phenotype}:")
    for idx, row in df_pheno.head(10).iterrows():
        print(f"    Feature {int(row['feature_index']):5} | ρ={row['rho']:7.3f} | q={row['q_value_stage1']:.2e} | {row['tier']}")

# Save results (compressed)
print("\n" + "="*80)
print("SAVING RESULTS")
print("="*80)

# Full results table
df_results.to_csv(os.path.join(script_dir, 'phase2_all_correlations.csv.gz'), index=False, compression='gzip')
print("✓ All correlations saved: phase2_all_correlations.csv.gz")

# Tier 1 hits (high confidence)
tier1 = df_results[df_results['tier'] == 'Tier1_High'].sort_values('q_value_stage1')
tier1.to_csv(os.path.join(script_dir, 'phase2_tier1_hits.csv.gz'), index=False, compression='gzip')
print(f"✓ Tier 1 hits ({len(tier1)} features): phase2_tier1_hits.csv.gz")

# Tier 1+2 hits (high + medium confidence)
tier12 = df_results[df_results['tier'].isin(['Tier1_High', 'Tier2_Medium'])].sort_values('q_value_stage1')
tier12.to_csv(os.path.join(script_dir, 'phase2_tier12_hits.csv.gz'), index=False, compression='gzip')
print(f"✓ Tier 1+2 hits ({len(tier12)} features): phase2_tier12_hits.csv.gz")

# Summary statistics (convert numpy types to Python native types for JSON)
summary_stats = {
    'n_features_analyzed': int(features.shape[1]),
    'n_phenotypes': int(len(phenotype_cols)),
    'n_total_tests': int(features.shape[1] * len(phenotype_cols)),
    'tier1_count': int((df_results['tier'] == 'Tier1_High').sum()),
    'tier2_count': int((df_results['tier'] == 'Tier2_Medium').sum()),
    'tier3_count': int((df_results['tier'] == 'Tier3_Exploratory').sum()),
    'strategy': decision['strategy'],
    'covariates': covariate_names,
    'phenotypes_used': phenotype_cols,
}

with open(os.path.join(script_dir, 'phase2_summary.json'), 'w') as f:
    json.dump(summary_stats, f, indent=2)

print("\n✓ Phase 2 complete!")
print(f"\nKey findings:")
print(f"  • {summary_stats['tier1_count']} high-confidence metabolite-color associations (Tier 1)")
print(f"  • {summary_stats['tier2_count']} medium-confidence associations (Tier 2)")
print(f"  • All {len(phenotype_cols)} color phenotypes analyzed")
print(f"  • FDR-corrected p-values with multiple testing control")

print(f"\nRecommended next steps:")
print(f"  1. Review Tier 1 hits (highest confidence)")
print(f"  2. Check effect size consistency across species (if stratified)")
print(f"  3. Cross-validation with holdout set")
print(f"  4. Annotate top features (known pigments?)")
print(f"  5. Create publication figures (heatmaps, volcano plots)")

