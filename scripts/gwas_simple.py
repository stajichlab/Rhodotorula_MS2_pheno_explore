#!/usr/bin/env python3
"""
GWAS: SNP-Brightness Association (Simplified Robust Version)
"""

import pandas as pd
import numpy as np
import gzip
import json
import os
from scipy import stats
import matplotlib.pyplot as plt

print("="*80)
print("GWAS: BRIGHTNESS × SNP ASSOCIATION")
print("="*80)

# Setup
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, '../results/gwas')
os.makedirs(results_dir, exist_ok=True)

vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'

print("\n[1/4] Loading data...")

# Load phenotype data
with gzip.open('./results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    brightness_data = pd.read_csv(f)

# Load master metadata
with gzip.open('./input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz', 'rt') as f:
    master_meta = pd.read_csv(f, sep='\t')

master_meta['dbvpg_num'] = master_meta['db_strain_id'].str.extract(r'DBVPG[:\s]+(\d+)')[0]
master_meta['vcf_sample'] = 'DBVPG_' + master_meta['dbvpg_num'].astype(str)

brightness_with_vcf = brightness_data.merge(
    master_meta[['filename', 'vcf_sample']],
    left_on='filename',
    right_on='filename',
    how='left'
)

# Get VCF sample names
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

# Create phenotype dict
pheno_dict = dict(zip(brightness_matched['vcf_sample'], brightness_matched['Median_ColorLab_L*Mean']))

# Create sample index mapping
sample_idx_map = {sample: idx for idx, sample in enumerate(vcf_samples)}
matched_indices = [sample_idx_map[s] for s in brightness_matched['vcf_sample'].values]

print("\n[2/4] Reading VCF and testing associations...")

results = []
snp_count = 0
tested_count = 0

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

        snp_count += 1

        # Parse genotypes (assuming GT is first field)
        genotypes = []
        for i in matched_indices:
            gt_str = fields[9 + i]
            gt = gt_str.split(':')[0]  # Get GT field

            # Convert to allele count
            if gt in ['0/0', '0|0', './.', '.|.']:
                geno = np.nan if gt in ['./. ', '.|.'] else 0
            elif gt in ['0/1', '1/0', '0|1', '1|0']:
                geno = 1
            elif gt in ['1/1', '1|1']:
                geno = 2
            else:
                geno = np.nan

            genotypes.append(geno)

        genotypes = np.array(genotypes)

        # Skip if too many missing
        if np.isnan(genotypes).sum() > len(genotypes) * 0.5:
            continue

        # Skip if no variance
        valid_geno = genotypes[~np.isnan(genotypes)]
        if len(valid_geno) < 5 or valid_geno.std() == 0:
            continue

        # Get phenotypes
        phenotypes = brightness_matched['Median_ColorLab_L*Mean'].values

        # Simple association test
        valid_idx = ~np.isnan(genotypes)
        if valid_idx.sum() < 5:
            continue

        x = genotypes[valid_idx]
        y = phenotypes[valid_idx]

        # Linear regression
        try:
            # Simple correlation/regression
            slope, intercept, r, pval, stderr = stats.linregress(x, y)

            results.append({
                'snp_id': snp_id,
                'ref': ref,
                'alt': alt,
                'effect': slope,
                'p_value': pval,
                'r_squared': r**2,
                'n_samples': len(x)
            })
            tested_count += 1
        except:
            pass

        if snp_count % 50000 == 0:
            print(f"  Processed {snp_count:,} SNPs, tested {tested_count}, found {len(results)} associations...")

print(f"\n  Total SNPs processed: {snp_count:,}")
print(f"  SNPs tested: {tested_count:,}")
print(f"  SNPs with associations: {len(results)}")

if len(results) == 0:
    print("\n✗ ERROR: No associations found. Check VCF file format.")
    exit(1)

print("\n[3/4] Computing statistics...")

results_df = pd.DataFrame(results)
results_df['neg_log10_p'] = -np.log10(results_df['p_value'].clip(lower=1e-300))

# Multiple testing correction
bonferroni_thresh = 0.05 / len(results_df)
results_df['sig_bonf'] = results_df['p_value'] < bonferroni_thresh

# FDR
results_sorted = results_df.sort_values('p_value').reset_index(drop=True)
results_sorted['fdr'] = (results_sorted['p_value'] * len(results_df) / (np.arange(len(results)) + 1)).clip(upper=1.0)
results_df['fdr'] = results_sorted['fdr'].sort_index().values
results_df['sig_fdr'] = results_df['fdr'] < 0.05

print(f"  Tested: {len(results_df)} SNPs")
print(f"  Bonferroni threshold: p < {bonferroni_thresh:.2e}")
print(f"  Significant (Bonf): {results_df['sig_bonf'].sum()}")
print(f"  Significant (FDR): {results_df['sig_fdr'].sum()}")

# Top hits
top20 = results_df.nsmallest(20, 'p_value')
print("\n  Top 20 SNPs:")
for idx, row in top20.iterrows():
    print(f"    {row['snp_id']:20} β={row['effect']:7.4f}  p={row['p_value']:.2e}  r²={row['r_squared']:.3f}")

print("\n[4/4] Saving results...")

# Save
results_df.to_csv(os.path.join(results_dir, 'gwas_results_all_snps.csv'), index=False)
top20.to_csv(os.path.join(results_dir, 'gwas_results_top_20.csv'), index=False)

# Visualize
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Manhattan
ax = axes[0, 0]
results_sorted = results_df.sort_values('neg_log10_p', ascending=False)
ax.scatter(range(len(results_sorted)), results_sorted['neg_log10_p'], alpha=0.5, s=15)
ax.axhline(-np.log10(bonferroni_thresh), color='red', linestyle='--', label='Bonferroni')
ax.set_xlabel('SNP')
ax.set_ylabel('-log10(p)')
ax.set_title('Manhattan Plot')
ax.legend()

# QQ
ax = axes[0, 1]
pvals_sorted = np.sort(results_df['p_value'].values)
expected = -np.log10(np.linspace(1.0/len(pvals_sorted), 1, len(pvals_sorted)))
observed = -np.log10(pvals_sorted)
ax.scatter(expected, observed, alpha=0.5, s=15)
ax.plot([0, expected.max()], [0, expected.max()], 'r--', label='Expected')
ax.set_xlabel('Expected -log10(p)')
ax.set_ylabel('Observed -log10(p)')
ax.set_title('Q-Q Plot')
ax.legend()

# Effect size
ax = axes[1, 0]
colors = ['red' if sig else 'blue' for sig in results_df['sig_bonf']]
ax.scatter(results_df['effect'], results_df['neg_log10_p'], alpha=0.5, c=colors, s=15)
ax.axhline(-np.log10(bonferroni_thresh), color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Effect Size (β)')
ax.set_ylabel('-log10(p)')
ax.set_title('Effect vs. Significance')

# Histogram
ax = axes[1, 1]
ax.hist(results_df['neg_log10_p'], bins=50, alpha=0.7, edgecolor='black')
ax.axvline(-np.log10(bonferroni_thresh), color='red', linestyle='--', label='Bonf')
ax.set_xlabel('-log10(p)')
ax.set_ylabel('Count')
ax.set_title('P-value Distribution')
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'gwas_summary.png'), dpi=200)

# Summary
summary = {
    'samples': len(brightness_matched),
    'snps_tested': len(results_df),
    'significant_bonferroni': int(results_df['sig_bonf'].sum()),
    'significant_fdr': int(results_df['sig_fdr'].sum()),
    'brightness_mean': float(brightness_matched['Median_ColorLab_L*Mean'].mean()),
    'brightness_std': float(brightness_matched['Median_ColorLab_L*Mean'].std())
}

with open(os.path.join(results_dir, 'gwas_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\n{'='*80}")
print("✓ GWAS COMPLETE")
print(f"{'='*80}")
print(f"Saved to: {results_dir}/")
print(f"  - gwas_results_all_snps.csv ({len(results_df)} SNPs)")
print(f"  - gwas_results_top_20.csv")
print(f"  - gwas_summary.png")
print(f"  - gwas_summary.json")
