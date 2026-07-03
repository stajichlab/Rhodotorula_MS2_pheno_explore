#!/usr/bin/env python3
"""
GWAS Step 2: SNP-Brightness Association Testing
Tests all SNPs for association with brightness phenotype
"""

import pandas as pd
import numpy as np
import gzip
import os
import json
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("STEP 2: GWAS - SNP-BRIGHTNESS ASSOCIATION TESTING")
print("="*80)

# Setup
os.makedirs('results/gwas', exist_ok=True)
vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'

print("\n[1/5] LOADING DATA")
print("-" * 80)

# Load phenotypes
with gzip.open('./results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    brightness = pd.read_csv(f)

# Load metadata for sample mapping
with gzip.open('./input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz', 'rt') as f:
    master = pd.read_csv(f, sep='\t')

master['dbvpg_num'] = master['db_strain_id'].str.extract(r'DBVPG[:\s]+(\d+)')[0]
master['vcf_sample'] = 'DBVPG_' + master['dbvpg_num'].astype(str)

brightness_vcf = brightness.merge(master[['filename', 'vcf_sample']], on='filename', how='left')

# Get VCF samples and filter
vcf_samples = []
with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            vcf_samples = line.strip().split('\t')[9:]
            break

matched = brightness_vcf[brightness_vcf['vcf_sample'].isin(vcf_samples)].copy()

print(f"Matched samples: {len(matched)}")
print(f"Brightness: {matched['Median_ColorLab_L*Mean'].mean():.2f} ± {matched['Median_ColorLab_L*Mean'].std():.2f} L*")

pheno_dict = dict(zip(matched['vcf_sample'], matched['Median_ColorLab_L*Mean']))
sample_idx_map = {sample: idx for idx, sample in enumerate(vcf_samples)}
matched_indices = [sample_idx_map[s] for s in matched['vcf_sample'].values]

print("\n[2/5] TESTING SNP ASSOCIATIONS")
print("-" * 80)

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

        # Parse genotypes
        genotypes = []
        for i in matched_indices:
            gt_str = fields[9 + i]
            gt = gt_str.split(':')[0]

            if gt in ['0/0', '0|0']:
                geno = 0
            elif gt in ['0/1', '1/0', '0|1', '1|0']:
                geno = 1
            elif gt in ['1/1', '1|1']:
                geno = 2
            else:
                try:
                    g = int(gt)
                    geno = g if g in [0, 1, 2] else np.nan
                except:
                    geno = np.nan

            genotypes.append(geno)

        genotypes = np.array(genotypes, dtype=float)

        # Filter SNPs
        if np.isnan(genotypes).sum() > len(genotypes) * 0.5:
            continue

        valid_geno = genotypes[~np.isnan(genotypes)]
        if len(valid_geno) < 5 or valid_geno.std() == 0:
            continue

        # Test association
        phenotypes = matched['Median_ColorLab_L*Mean'].values
        valid_idx = ~np.isnan(genotypes)
        if valid_idx.sum() < 5:
            continue

        x = genotypes[valid_idx]
        y = phenotypes[valid_idx]

        try:
            slope, intercept, r, pval, stderr = stats.linregress(x, y)
            results.append({
                'snp_id': snp_id,
                'chr': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'beta': slope,
                'se': stderr,
                'p_value': pval,
                'r_squared': r**2,
                'n_samples': len(x)
            })
            tested_count += 1
        except:
            pass

        if snp_count % 100000 == 0:
            print(f"  Processed {snp_count:,} SNPs, tested {tested_count}, found {len(results)}")

print(f"\nTotal SNPs processed: {snp_count:,}")
print(f"SNPs with valid tests: {tested_count:,}")

if len(results) == 0:
    print("ERROR: No associations found!")
    exit(1)

print("\n[3/5] CALCULATING STATISTICS")
print("-" * 80)

results_df = pd.DataFrame(results)
results_df['neg_log10_p'] = -np.log10(results_df['p_value'].clip(lower=1e-300))

# Multiple testing corrections
bonferroni_thresh = 0.05 / len(results_df)
results_df['sig_bonf'] = results_df['p_value'] < bonferroni_thresh

# Benjamini-Hochberg FDR
results_sorted = results_df.sort_values('p_value').reset_index(drop=True)
results_sorted['rank'] = np.arange(1, len(results_sorted) + 1)
results_sorted['q_value'] = (results_sorted['p_value'] * len(results_df) / results_sorted['rank']).clip(upper=1.0)
results_df['q_value'] = results_sorted.sort_index()['q_value'].values
results_df['sig_fdr'] = results_df['q_value'] < 0.05

print(f"SNPs tested: {len(results_df):,}")
print(f"Bonferroni threshold: {bonferroni_thresh:.2e}")
print(f"Significant (Bonferroni): {results_df['sig_bonf'].sum()}")
print(f"Significant (FDR): {results_df['sig_fdr'].sum()}")

print("\n[4/5] SAVING RESULTS")
print("-" * 80)

# Save results
results_df.to_csv('results/gwas/gwas_all_results.csv', index=False)
results_df.nsmallest(20, 'p_value').to_csv('results/gwas/gwas_top_20.csv', index=False)

if results_df['sig_bonf'].sum() > 0:
    results_df[results_df['sig_bonf']].sort_values('p_value').to_csv('results/gwas/gwas_sig_bonf.csv', index=False)

if results_df['sig_fdr'].sum() > 0:
    results_df[results_df['sig_fdr']].sort_values('p_value').to_csv('results/gwas/gwas_sig_fdr.csv', index=False)

summary = {
    'n_samples': len(matched),
    'n_snps_tested': len(results_df),
    'n_sig_bonf': int(results_df['sig_bonf'].sum()),
    'n_sig_fdr': int(results_df['sig_fdr'].sum()),
    'bonf_threshold': float(bonferroni_thresh),
    'brightness_mean': float(matched['Median_ColorLab_L*Mean'].mean()),
    'brightness_std': float(matched['Median_ColorLab_L*Mean'].std()),
    'p_value_min': float(results_df['p_value'].min()),
    'p_value_median': float(results_df['p_value'].median())
}

with open('results/gwas/gwas_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(f"✓ Results saved to results/gwas/")
print(f"  - gwas_all_results.csv ({len(results_df):,} SNPs)")
print(f"  - gwas_top_20.csv")
if results_df['sig_bonf'].sum() > 0:
    print(f"  - gwas_sig_bonf.csv ({results_df['sig_bonf'].sum()} SNPs)")
if results_df['sig_fdr'].sum() > 0:
    print(f"  - gwas_sig_fdr.csv ({results_df['sig_fdr'].sum()} SNPs)")
print(f"  - gwas_summary.json")

print("\n[5/5] TOP RESULTS PREVIEW")
print("-" * 80)

top20 = results_df.nsmallest(20, 'p_value')
print("\nTop 20 most significant SNPs:")
print(top20[['snp_id', 'chr', 'beta', 'p_value', 'q_value']].to_string(index=False))

print("\n" + "="*80)
print("✓ GWAS ANALYSIS COMPLETE")
print("="*80)
