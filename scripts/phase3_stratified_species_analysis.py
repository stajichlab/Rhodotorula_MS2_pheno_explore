#!/usr/bin/env python3
"""
PHASE 3: STRATIFIED SPECIES ANALYSIS (optimized)
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
print("PHASE 3: STRATIFIED SPECIES ANALYSIS (OPTIMIZED)")
print("="*80)

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, '../results/phase')

with open(os.path.join(output_dir, 'phase0_decision.json'), 'r') as f:
    decision = json.load(f)

features = pd.read_csv(os.path.join(output_dir, 'phase1_features_filtered.csv.gz'), compression='gzip')
metadata = pd.read_csv(os.path.join(output_dir, 'phase1_phenotype_data.csv.gz'), compression='gzip')
phenotype_cols = decision['phenotypes_to_use']

# Get feature column names (they are numeric strings)
feature_cols = features.columns.tolist()

print(f"Data: {metadata.shape[0]} samples, {features.shape[1]} features")

# Filter to major species
species_counts = metadata['Species'].value_counts()
major_species = species_counts[species_counts >= 3].index.tolist()
mask = metadata['Species'].isin(major_species)
features_major = features[mask].reset_index(drop=True)
metadata_major = metadata[mask].reset_index(drop=True)

print(f"Major species: {len(major_species)} | Samples: {len(metadata_major)}")

# ============================================================================
# PART 1: BETWEEN-SPECIES (discriminant features)
# ============================================================================
print("\n" + "="*80)
print("PART 1: BETWEEN-SPECIES ANALYSIS")
print("="*80)

species_abundance = np.array([features_major[metadata_major['Species']==s].mean(axis=0).values 
                              for s in major_species])

feature_var = species_abundance.var(axis=0)
top_disc_idx = np.argsort(feature_var)[::-1][:200]
top_disc_cols = [feature_cols[i] for i in top_disc_idx]

discriminant_df = pd.DataFrame({
    'feature_index': top_disc_idx,
    'feature_column': top_disc_cols,
    'variance_across_species': feature_var[top_disc_idx]
})
discriminant_df = discriminant_df.sort_values('variance_across_species', ascending=False)
discriminant_df.to_csv(os.path.join(output_dir, 'phase3_discriminant_features.csv'), index=False)

print(f"✓ Identified {len(discriminant_df)} top discriminant features")

# Phenotype statistics by species
pheno_by_species = []
for species in major_species:
    mask_sp = metadata_major['Species'] == species
    for pheno in phenotype_cols:
        values = metadata_major.loc[mask_sp, pheno]
        pheno_by_species.append({
            'species': species,
            'phenotype': pheno,
            'mean': values.mean(),
            'std': values.std(),
            'n_samples': int(mask_sp.sum()),
            'min': values.min(),
            'max': values.max()
        })

pheno_by_species_df = pd.DataFrame(pheno_by_species)
pheno_by_species_df.to_csv(os.path.join(output_dir, 'phase3_phenotype_by_species.csv'), index=False)
print(f"✓ Phenotype statistics by species saved")

# ============================================================================
# PART 2: WITHIN-SPECIES CORRELATIONS (discriminant features only)
# ============================================================================
print("\n" + "="*80)
print("PART 2: WITHIN-SPECIES CORRELATIONS")
print("="*80)

def spearman_partial_corr(X, Y, covariates=None):
    X_r = rankdata(X)
    Y_r = rankdata(Y)
    
    if not covariates:
        rho, pval = spearmanr(X_r, Y_r)
        return rho, pval
    
    from scipy.stats import linregress
    X_res = X_r.copy()
    Y_res = Y_r.copy()
    
    for cov in covariates:
        slope_x, intercept_x = linregress(cov, X_r)[:2]
        X_res = X_r - (intercept_x + slope_x * cov)
        slope_y, intercept_y = linregress(cov, Y_r)[:2]
        Y_res = Y_r - (intercept_y + slope_y * cov)
    
    rho, pval = spearmanr(X_res, Y_res)
    return rho, pval

plate_all = np.array(pd.factorize(metadata_major['Library Plate'])[0], dtype=float)
all_results = []

for sp_idx, species in enumerate(major_species):
    print(f"  {sp_idx+1}/{len(major_species)}: {species}")
    
    mask_sp = np.array(metadata_major['Species'] == species)
    feat_sp = features_major[mask_sp].reset_index(drop=True)
    meta_sp = metadata_major[mask_sp].reset_index(drop=True)
    plate_sp = plate_all[mask_sp]
    
    if mask_sp.sum() < 4:
        continue
    
    covs = [plate_sp] if len(np.unique(plate_sp)) > 1 else []
    results = []
    
    # Compute correlations for discriminant features only
    for feat_col in top_disc_cols:
        for pheno in phenotype_cols:
            y = meta_sp[pheno].values
            x = feat_sp[feat_col].values
            
            valid = ~np.isnan(y)
            x_v = x[valid]
            y_v = y[valid]
            
            if np.std(x_v) == 0 or len(y_v) < 4:
                continue
            
            cov_v = [c[valid] for c in covs] if covs else []
            
            try:
                rho, pval = spearman_partial_corr(x_v, y_v, cov_v)
                results.append({
                    'species': species,
                    'phenotype': pheno,
                    'feature': feat_col,
                    'rho': rho,
                    'pval': pval,
                    'n_samples': len(x_v)
                })
            except:
                pass
    
    if results:
        df_sp = pd.DataFrame(results)
        reject, q_vals, _, _ = multipletests(df_sp['pval'], method='fdr_bh', alpha=0.05)
        df_sp['q_value'] = q_vals
        df_sp['reject'] = reject
        all_results.append(df_sp)
        print(f"    {len(df_sp)} corr | {reject.sum()} sig")

df_within = pd.concat(all_results, ignore_index=True)
df_within.to_csv(os.path.join(output_dir, 'phase3_within_species_correlations.csv.gz'), 
                 index=False, compression='gzip')

summary = {
    'n_major_species': len(major_species),
    'n_discriminant_features': len(discriminant_df),
    'n_within_species_correlations': len(df_within),
    'n_significant': int((df_within['reject']).sum())
}

with open(os.path.join(output_dir, 'phase3_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print("\n✓ PHASE 3 COMPLETE\n")

