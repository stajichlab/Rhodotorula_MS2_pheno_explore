#!/usr/bin/env python3
"""
GWAS: SNP-Brightness Association (Optimized Fast Version)
Uses cyvcf2 for efficient VCF parsing
"""

import pandas as pd
import numpy as np
import gzip
import json
import os
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("GWAS: BRIGHTNESS × CAROTENOID GENES (OPTIMIZED)")
print("="*80)

# Setup paths
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, '../results/gwas')
os.makedirs(results_dir, exist_ok=True)

vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'

print("\n[1/5] Loading data...")

# Load phenotype data
with gzip.open('./results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    brightness_data = pd.read_csv(f)

# Load master metadata and map to VCF samples
with gzip.open('./input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz', 'rt') as f:
    master_meta = pd.read_csv(f, sep='\t')

master_meta['dbvpg_num'] = master_meta['db_strain_id'].str.extract(r'DBVPG[:\s]+(\d+)')[0]
master_meta['vcf_sample'] = 'DBVPG_' + master_meta['dbvpg_num'].astype(str)

brightness_with_vcf = brightness_data.merge(
    master_meta[['filename', 'vcf_sample', 'dbvpg_num']],
    left_on='filename',
    right_on='filename',
    how='left'
)

# Get VCF sample order
vcf_samples = []
with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            vcf_samples = line.strip().split('\t')[9:]
            break

vcf_sample_set = set(vcf_samples)
brightness_matched = brightness_with_vcf[brightness_with_vcf['vcf_sample'].isin(vcf_sample_set)].copy()

print(f"  Matched samples: {len(brightness_matched)}")
print(f"  Brightness: {brightness_matched['Median_ColorLab_L*Mean'].mean():.2f} ± {brightness_matched['Median_ColorLab_L*Mean'].std():.2f}")

pheno_dict = dict(zip(brightness_matched['vcf_sample'], brightness_matched['Median_ColorLab_L*Mean']))
sample_to_idx = {sample: i for i, sample in enumerate(vcf_samples)}

print("\n[2/5] Computing PCA for population structure...")

# Read VCF and compute PCA on subset of SNPs
genotypes_all = []
snp_ids = []

snp_count = 0
with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])

        # Sample every 100th SNP for PCA (speed up)
        if snp_count % 100 != 0:
            snp_count += 1
            continue

        gt_field_idx = fields[8].split(':').index('GT')
        genotypes = []

        for sample_gt in fields[9:]:
            gt = sample_gt.split(':')[gt_field_idx]
            if gt in ['0/0', '0|0']:
                genotypes.append(0)
            elif gt in ['0/1', '1/0', '0|1', '1|0']:
                genotypes.append(1)
            elif gt in ['1/1', '1|1']:
                genotypes.append(2)
            else:
                genotypes.append(np.nan)

        genotypes_all.append(genotypes)
        snp_ids.append(f"{chrom}:{pos}")
        snp_count += 1

genotype_matrix = np.array(genotypes_all).T  # samples × SNPs
print(f"  SNPs for PCA: {len(snp_ids)}")

# Normalize and compute PCA
genotype_matrix[np.isnan(genotype_matrix)] = 0
col_means = np.mean(genotype_matrix, axis=0)
col_stds = np.std(genotype_matrix, axis=0)
col_stds[col_stds == 0] = 1

genotype_matrix_std = (genotype_matrix - col_means) / col_stds

# SVD
U, s, _ = np.linalg.svd(genotype_matrix_std, full_matrices=False)
pcs = U[:, :3]
print(f"  PCA variance: {np.cumsum(s[:3]**2 / (s**2).sum())[-1]*100:.1f}%")

sample_to_pc = {sample: pcs[i] for i, sample in enumerate(vcf_samples)}

print("\n[3/5] Testing SNP associations...")

# Full pass through VCF for association testing
results = []
snp_count = 0

with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]

        snp_id = f"{chrom}:{pos}"
        gt_field_idx = fields[8].split(':').index('GT')

        # Extract genotypes for matched samples
        x_list = []
        y_list = []
        pc_list = []

        for i, sample in enumerate(vcf_samples):
            if sample in pheno_dict:
                gt = fields[9 + i].split(':')[gt_field_idx]

                if gt in ['0/0', '0|0']:
                    geno = 0
                elif gt in ['0/1', '1/0', '0|1', '1|0']:
                    geno = 1
                elif gt in ['1/1', '1|1']:
                    geno = 2
                else:
                    continue

                x_list.append(geno)
                y_list.append(pheno_dict[sample])
                pc_list.append(sample_to_pc[sample])

        if len(x_list) < 5:
            continue

        x = np.array(x_list)
        y = np.array(y_list)
        pcs = np.array(pc_list)

        # Skip if no variation in genotype
        if x.std() == 0:
            continue

        # Regression: brightness ~ SNP + PC1 + PC2 + PC3
        X = np.column_stack([x, pcs])

        try:
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            residuals = y - X @ beta

            n = len(y)
            k = X.shape[1]
            mse = np.sum(residuals**2) / (n - k)
            var_covar = mse * np.linalg.inv(X.T @ X)
            se = np.sqrt(np.diag(var_covar))

            t_stat = beta[0] / se[0]
            p_value = 2 * (1 - stats.t.cdf(abs(t_stat), n - k))

            results.append({
                'snp_id': snp_id,
                'ref': ref,
                'alt': alt,
                'beta': beta[0],
                'se': se[0],
                'p_value': p_value,
                'n_samples': n
            })
        except:
            pass

        snp_count += 1
        if snp_count % 10000 == 0:
            print(f"  Tested {snp_count:,} SNPs, found {len(results)} with valid tests...")

results_df = pd.DataFrame(results)
print(f"  Total SNPs tested: {len(results_df)}")

print("\n[4/5] Multiple testing correction...")

# Add statistics
results_df['neg_log10_p'] = -np.log10(results_df['p_value'])

# Bonferroni
bonferroni_thresh = 0.05 / len(results_df)
results_df['significant_bonf'] = results_df['p_value'] < bonferroni_thresh

# BH FDR
results_sorted = results_df.sort_values('p_value').reset_index(drop=True)
results_sorted['rank'] = np.arange(1, len(results_sorted) + 1)
results_sorted['q_value_bh'] = (results_sorted['p_value'] * len(results_df) / results_sorted['rank']).clip(upper=1.0)

# Merge back
results_df = results_sorted.sort_index()
results_df['q_value_bh'] = results_sorted['q_value_bh'].sort_index().values
results_df['significant_fdr'] = results_df['q_value_bh'] < 0.05

print(f"  Bonferroni threshold: p < {bonferroni_thresh:.2e}")
print(f"  Significant (Bonf): {results_df['significant_bonf'].sum()}")
print(f"  Significant (FDR): {results_df['significant_fdr'].sum()}")

print("\n[5/5] Generating outputs...")

# Top SNPs
top_snps = results_df.nsmallest(20, 'p_value')
print("\nTop 20 SNPs:")
print(top_snps[['snp_id', 'beta', 'p_value']].to_string(index=False))

# Save results
results_df.to_csv(os.path.join(results_dir, 'gwas_results_all_snps.csv'), index=False)
top_snps.to_csv(os.path.join(results_dir, 'gwas_results_top_20.csv'), index=False)

# Visualizations
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Manhattan plot
ax = axes[0, 0]
results_sorted = results_df.sort_values('neg_log10_p', ascending=False)
ax.scatter(range(len(results_sorted)), results_sorted['neg_log10_p'], alpha=0.6, s=20)
ax.axhline(-np.log10(bonferroni_thresh), color='red', linestyle='--', label=f'Bonferroni')
ax.set_xlabel('SNP (ranked)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('Manhattan Plot')
ax.legend()

# QQ plot
ax = axes[0, 1]
expected = -np.log10(np.linspace(0.001, 1, len(results_df)))
observed = np.sort(-np.log10(results_df['p_value'].values))
ax.scatter(expected, observed, alpha=0.6, s=20)
ax.plot([0, expected.max()], [0, expected.max()], 'r--')
ax.set_xlabel('Expected -log10(p)')
ax.set_ylabel('Observed -log10(p)')
ax.set_title('Q-Q Plot')

# Effect size
ax = axes[1, 0]
colors = ['red' if sig else 'blue' for sig in results_df['significant_bonf']]
ax.scatter(results_df['beta'], -np.log10(results_df['p_value']), alpha=0.5, c=colors, s=20)
ax.axhline(-np.log10(bonferroni_thresh), color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Effect Size (β)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('Effect Size vs Significance')

# Histogram
ax = axes[1, 1]
ax.hist(-np.log10(results_df['p_value']), bins=50, alpha=0.7, edgecolor='black')
ax.axvline(-np.log10(bonferroni_thresh), color='red', linestyle='--')
ax.set_xlabel('-log10(p-value)')
ax.set_ylabel('Frequency')
ax.set_title('P-value Distribution')

plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'gwas_summary.png'), dpi=200)
print(f"  Saved: gwas_summary.png")

# Summary JSON
summary = {
    'n_matched_samples': len(brightness_matched),
    'n_snps_tested': len(results_df),
    'n_significant_bonf': int(results_df['significant_bonf'].sum()),
    'n_significant_fdr': int(results_df['significant_fdr'].sum()),
    'bonferroni_threshold': float(bonferroni_thresh),
    'brightness_mean': float(brightness_matched['Median_ColorLab_L*Mean'].mean()),
    'brightness_std': float(brightness_matched['Median_ColorLab_L*Mean'].std())
}

with open(os.path.join(results_dir, 'gwas_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\n{'='*80}")
print("✓ GWAS COMPLETE")
print(f"{'='*80}")
print(f"Results: {results_dir}/")
print(f"  - gwas_results_all_snps.csv")
print(f"  - gwas_results_top_20.csv")
print(f"  - gwas_summary.png")
print(f"  - gwas_summary.json")
