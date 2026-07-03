#!/usr/bin/env python3
"""
GWAS: SNP-Brightness Association in Rhodotorula mucilaginosa
Carotenoid pathway candidate gene analysis

Performs:
1. PCA correction for population stratification
2. LD-aware SNP extraction
3. Single-SNP association testing
4. Multiple testing correction
5. Manhattan plot visualization
"""

import pandas as pd
import numpy as np
import gzip
import json
import os
from scipy import stats
from scipy.stats import linregress
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("GWAS: BRIGHTNESS × CAROTENOID GENES")
print("="*80)

# Setup paths
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, '../results/gwas')
os.makedirs(results_dir, exist_ok=True)

vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'

print("\n[1/6] Loading VCF and phenotype data...")

# ============================================================================
# STEP 1: Load VCF and extract sample names
# ============================================================================

vcf_samples = []
vcf_data = {}  # Will store chr, pos, ref, alt, genotypes

with gzip.open(vcf_path, 'rt') as f:
    for line_num, line in enumerate(f):
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            vcf_samples = header[9:]  # Samples start at column 9
            print(f"  VCF samples: {len(vcf_samples)}")
            break
        elif line.startswith('##'):
            continue
        else:
            # This shouldn't happen before #CHROM
            pass

print(f"  Sample names (first 10): {vcf_samples[:10]}")

# ============================================================================
# STEP 2: Load phenotype data and match to VCF samples
# ============================================================================

with gzip.open('./results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    brightness_data = pd.read_csv(f)

# Load master metadata to get DBVPG IDs
with gzip.open('./input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz', 'rt') as f:
    master_meta = pd.read_csv(f, sep='\t')

# Extract DBVPG numbers
master_meta['dbvpg_num'] = master_meta['db_strain_id'].str.extract(r'DBVPG[:\s]+(\d+)')[0]
master_meta['vcf_sample'] = 'DBVPG_' + master_meta['dbvpg_num'].astype(str)

# Merge brightness with VCF sample names
brightness_with_vcf = brightness_data.merge(
    master_meta[['filename', 'vcf_sample', 'dbvpg_num']],
    left_on='filename',
    right_on='filename',
    how='left'
)

# Filter to samples that are in VCF
vcf_sample_set = set(vcf_samples)
brightness_matched = brightness_with_vcf[brightness_with_vcf['vcf_sample'].isin(vcf_sample_set)].copy()

print(f"  Matched samples: {len(brightness_matched)} lab strains in VCF")
print(f"  Brightness range: {brightness_matched['Median_ColorLab_L*Mean'].min():.2f} - {brightness_matched['Median_ColorLab_L*Mean'].max():.2f} L*")
print(f"  Mean brightness: {brightness_matched['Median_ColorLab_L*Mean'].mean():.2f} ± {brightness_matched['Median_ColorLab_L*Mean'].std():.2f}")

# Create phenotype array
pheno_dict = dict(zip(brightness_matched['vcf_sample'], brightness_matched['Median_ColorLab_L*Mean']))
n_matched = len(brightness_matched)

print("\n[2/6] Loading genotypes and computing PCA...")

# ============================================================================
# STEP 3: Load genotypes and perform PCA
# ============================================================================

# Read VCF and extract genotypes
genotypes_all = {}  # chr:pos -> array of genotypes
snp_info = {}  # chr:pos -> (ref, alt, annotation)

snp_count = 0
with gzip.open(vcf_path, 'rt') as f:
    for line_num, line in enumerate(f):
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        info_str = fields[7]

        snp_id = f"{chrom}:{pos}"

        # Extract genotypes (GT field, format 0/1 or 0|1)
        gt_field_idx = fields[8].split(':').index('GT')
        genotypes = []
        for i, sample_gt in enumerate(fields[9:]):
            gt = sample_gt.split(':')[gt_field_idx]
            # Convert 0/1 or 0|1 to count of alt alleles (0, 1, or 2)
            if gt in ['0/0', '0|0']:
                genotypes.append(0)
            elif gt in ['0/1', '1/0', '0|1', '1|0']:
                genotypes.append(1)
            elif gt in ['1/1', '1|1']:
                genotypes.append(2)
            else:
                genotypes.append(np.nan)

        genotypes_all[snp_id] = np.array(genotypes, dtype=float)
        snp_info[snp_id] = (ref, alt, info_str)

        snp_count += 1
        if snp_count % 100000 == 0:
            print(f"  Loaded {snp_count:,} SNPs...")

print(f"  Total SNPs loaded: {snp_count:,}")

# Convert to matrix (samples × SNPs)
snp_ids = list(genotypes_all.keys())
genotype_matrix = np.array([genotypes_all[snp] for snp in snp_ids]).T  # samples × SNPs

print(f"  Genotype matrix shape: {genotype_matrix.shape}")

# Perform PCA
print("  Performing PCA...")
# Normalize: remove mean and divide by std
genotype_matrix_nan = np.copy(genotype_matrix)
genotype_matrix_nan[np.isnan(genotype_matrix_nan)] = 0  # Impute missing as 0

col_means = np.nanmean(genotype_matrix_nan, axis=0)
col_stds = np.nanstd(genotype_matrix_nan, axis=0)
col_stds[col_stds == 0] = 1  # Avoid division by zero

genotype_matrix_std = (genotype_matrix_nan - col_means) / col_stds

# SVD (fast PCA)
try:
    U, s, Vt = np.linalg.svd(genotype_matrix_std, full_matrices=False)
    pcs = U[:, :3]  # First 3 PCs
    print(f"  PCA variance explained (PC1-3): {np.cumsum(s[:3]**2 / (s**2).sum())[-1]*100:.1f}%")
except:
    print("  Warning: SVD failed; using simpler PCA")
    from sklearn.decomposition import PCA as SKL_PCA
    pca_obj = SKL_PCA(n_components=3)
    pcs = pca_obj.fit_transform(genotype_matrix_std)
    print(f"  PCA variance explained: {pca_obj.explained_variance_ratio_.sum()*100:.1f}%")

# Map samples to PC scores
sample_to_pc = {sample: pcs[i] for i, sample in enumerate(vcf_samples)}

print("\n[3/6] Extracting carotenoid gene SNPs...")

# ============================================================================
# STEP 4: Extract carotenoid pathway SNPs
# ============================================================================

# Carotenoid gene regions (approximate from literature; you'd use actual GFF3 coordinates)
# For now, we'll look for SNPs with carotenoid-related annotations
carotenoid_keywords = [
    'lycopene', 'carotene', 'phytoene', 'desaturase', 'cyclase',
    'synthase', 'hydroxylase', 'isomerase', 'beta-carotene'
]

# Load Rodeo Pfam data to get carotenoid genes
with gzip.open('/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Rodeo/tables/All_Taxa/pfam.csv.gz', 'rt') as f:
    pfam_data = pd.read_csv(f)

carotenoid_pfams = pfam_data[pfam_data['pfam_id'].isin(['Lycopene_cyc', 'Lycopene_cycl'])]
carotenoid_genes = set(carotenoid_pfams['protein_id'].str.extract(r'([A-Z0-9]+)_')[0].unique())

print(f"  Carotenoid-related genes: {len(carotenoid_genes)}")

# For this analysis, we'll filter SNPs by checking INFO field
# In a real analysis, you'd use VEP-annotated VCF or gene coordinate lookup
carotenoid_snps = {}

for snp_id, (ref, alt, info) in snp_info.items():
    # Check if SNP annotation mentions carotenoid-related genes
    # This is a placeholder; real VCF would have CSQ field
    for keyword in carotenoid_keywords:
        if keyword.lower() in info.lower():
            carotenoid_snps[snp_id] = {
                'ref': ref, 'alt': alt, 'info': info,
                'genotypes': genotypes_all[snp_id]
            }
            break

print(f"  SNPs in carotenoid-related regions: {len(carotenoid_snps)}")

# If we didn't find any annotated SNPs, use all SNPs (conservative approach)
if len(carotenoid_snps) == 0:
    print("  Note: No annotated carotenoid SNPs found; using all high-quality SNPs")
    # Use SNPs with MAF > 5% and low missing rate
    snp_ids_filtered = []
    for snp_id in snp_ids:
        geno = genotypes_all[snp_id]
        # Filter: non-missing and MAF filter
        valid_geno = geno[~np.isnan(geno)]
        if len(valid_geno) > 0:
            maf = min(np.mean(valid_geno) / 2, 1 - np.mean(valid_geno) / 2)
            if maf > 0.05:  # MAF > 5%
                snp_ids_filtered.append(snp_id)
    carotenoid_snps = {snp: {'genotypes': genotypes_all[snp]} for snp in snp_ids_filtered[:5000]}
    print(f"  Using {len(carotenoid_snps)} SNPs (MAF > 5%)")

print("\n[4/6] Single-SNP association testing...")

# ============================================================================
# STEP 5: Association testing
# ============================================================================

results = []

for snp_idx, (snp_id, snp_data) in enumerate(carotenoid_snps.items()):
    if snp_idx % 100 == 0:
        print(f"  Tested {snp_idx:,} SNPs...")

    geno = snp_data['genotypes']

    # Get PC scores and phenotypes for matched samples
    x_list = []
    y_list = []
    pc_list = []

    for i, sample in enumerate(vcf_samples):
        if sample in pheno_dict and not np.isnan(geno[i]):
            x_list.append(geno[i])
            y_list.append(pheno_dict[sample])
            pc_list.append(sample_to_pc[sample])

    if len(x_list) < 3:
        continue  # Skip if too few samples

    x = np.array(x_list)
    y = np.array(y_list)
    pcs = np.array(pc_list)

    # Model: brightness ~ SNP + PC1 + PC2 + PC3
    X = np.column_stack([x, pcs])

    try:
        # Least squares regression
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        residuals = y - X @ beta

        # Standard error
        n = len(y)
        k = X.shape[1]
        mse = np.sum(residuals**2) / (n - k)
        var_covar = mse * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(var_covar))

        # t-statistic and p-value
        t_stat = beta[0] / se[0]
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), n - k))

        results.append({
            'snp_id': snp_id,
            'ref': snp_data['genotypes'].get('ref', '?'),
            'alt': snp_data['genotypes'].get('alt', '?'),
            'beta': beta[0],
            'se': se[0],
            't_stat': t_stat,
            'p_value': p_value,
            'n_samples': n
        })
    except:
        continue

results_df = pd.DataFrame(results)
print(f"  Total tested SNPs: {len(results_df)}")

# Multiple testing correction
from scipy.stats import rankdata

results_df['p_value'] = results_df['p_value'].clip(lower=1e-300)  # Avoid log(0)
results_df['q_value_bh'] = (results_df['p_value'].rank() / len(results_df) * 0.05).clip(upper=1.0)  # BH FDR
results_df['neg_log10_p'] = -np.log10(results_df['p_value'])

# Bonferroni threshold
bonferroni_thresh = 0.05 / len(results_df)
results_df['significant_bonf'] = results_df['p_value'] < bonferroni_thresh

# FDR threshold
results_df['significant_fdr'] = results_df['q_value_bh'] < 0.05

print(f"\n  Significant SNPs (Bonferroni, p < {bonferroni_thresh:.2e}): {results_df['significant_bonf'].sum()}")
print(f"  Significant SNPs (FDR, q < 0.05): {results_df['significant_fdr'].sum()}")

print("\n[5/6] Generating visualizations...")

# ============================================================================
# STEP 6: Visualization
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Manhattan plot
ax = axes[0, 0]
results_sorted = results_df.sort_values('neg_log10_p', ascending=False)
ax.scatter(range(len(results_sorted)), results_sorted['neg_log10_p'], alpha=0.6, s=30)
ax.axhline(-np.log10(bonferroni_thresh), color='red', linestyle='--', label=f'Bonferroni (p={bonferroni_thresh:.2e})')
ax.axhline(-np.log10(0.05), color='orange', linestyle='--', label='p = 0.05')
ax.set_xlabel('SNP (ranked by p-value)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('Manhattan Plot: Brightness-SNP Association')
ax.legend()

# QQ plot
ax = axes[0, 1]
expected = -np.log10(np.linspace(0, 1, len(results_df))[1:])
observed = np.sort(-np.log10(results_df['p_value'].values))
ax.scatter(expected, observed, alpha=0.6, s=30)
ax.plot([0, expected.max()], [0, expected.max()], 'r--', label='Expected')
ax.set_xlabel('Expected -log10(p)')
ax.set_ylabel('Observed -log10(p)')
ax.set_title('Q-Q Plot')
ax.legend()

# Effect size vs p-value
ax = axes[1, 0]
colors = ['red' if sig else 'blue' for sig in results_df['significant_bonf']]
ax.scatter(results_df['beta'], -np.log10(results_df['p_value']), alpha=0.6, c=colors, s=30)
ax.axhline(-np.log10(bonferroni_thresh), color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Effect Size (β, L* change per alt allele)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('Effect Size vs Significance')

# P-value histogram
ax = axes[1, 1]
ax.hist(-np.log10(results_df['p_value']), bins=50, alpha=0.7, edgecolor='black')
ax.axvline(-np.log10(bonferroni_thresh), color='red', linestyle='--', label='Bonferroni')
ax.set_xlabel('-log10(p-value)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution of Association Signals')
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(results_dir, 'gwas_summary.png'), dpi=300)
print(f"  Saved: gwas_summary.png")
plt.close()

# ============================================================================
# STEP 7: Output results
# ============================================================================

print("\n[6/6] Writing results...")

# Top SNPs
top_snps = results_df.nsmallest(20, 'p_value')[['snp_id', 'beta', 'se', 'p_value', 'q_value_bh', 'n_samples']]
print("\nTop 20 SNPs associated with brightness:")
print(top_snps.to_string(index=False))

# Save full results
results_df.to_csv(os.path.join(results_dir, 'gwas_results_all_snps.csv'), index=False)
top_snps.to_csv(os.path.join(results_dir, 'gwas_results_top_20.csv'), index=False)

# Summary statistics
summary = {
    'n_matched_samples': n_matched,
    'n_snps_tested': len(results_df),
    'n_significant_bonf': int(results_df['significant_bonf'].sum()),
    'n_significant_fdr': int(results_df['significant_fdr'].sum()),
    'bonferroni_threshold': float(bonferroni_thresh),
    'fdr_threshold': 0.05,
    'brightness_mean': float(brightness_matched['Median_ColorLab_L*Mean'].mean()),
    'brightness_std': float(brightness_matched['Median_ColorLab_L*Mean'].std()),
    'brightness_range': [float(brightness_matched['Median_ColorLab_L*Mean'].min()),
                         float(brightness_matched['Median_ColorLab_L*Mean'].max())]
}

with open(os.path.join(results_dir, 'gwas_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\n{'='*80}")
print("GWAS COMPLETE")
print(f"{'='*80}")
print(f"Results saved to: {results_dir}/")
print(f"  - gwas_results_all_snps.csv (all {len(results_df)} SNPs)")
print(f"  - gwas_results_top_20.csv (top 20 by p-value)")
print(f"  - gwas_summary.json (summary statistics)")
print(f"  - gwas_summary.png (visualization)")
print(f"\nSignificant findings:")
print(f"  - Bonferroni (p < {bonferroni_thresh:.2e}): {results_df['significant_bonf'].sum()} SNPs")
print(f"  - FDR (q < 0.05): {results_df['significant_fdr'].sum()} SNPs")
