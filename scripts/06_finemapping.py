#!/usr/bin/env python3
"""
GWAS Step 6: Fine-Mapping on Scaffold 10
Narrow down to likely causal variants using LD and Bayesian credible sets
"""

import pandas as pd
import numpy as np
import gzip
import os
from scipy import stats

print("="*90)
print("STEP 6: FINE-MAPPING ON SCAFFOLD 10")
print("="*90)

# Load all GWAS results
print("\n[1/6] LOADING GWAS RESULTS")
print("-" * 90)

gwas_df = pd.read_csv('results/gwas/gwas_all_results.csv')

# Focus on scaffold_10
scaffold10_snps = gwas_df[gwas_df['chr'] == 'scaffold_10'].copy()
print(f"Total SNPs on scaffold_10: {len(scaffold10_snps):,}")

# Define the region around top signal (±500kb window around scaffold_10:144566)
top_snp_pos = 144566
window = 500000

region_snps = scaffold10_snps[
    (scaffold10_snps['pos'] >= top_snp_pos - window) &
    (scaffold10_snps['pos'] <= top_snp_pos + window)
].copy()

print(f"SNPs in ±{window/1000:.0f}kb window: {len(region_snps):,}")
print(f"Region: scaffold_10:{top_snp_pos - window:,} - {top_snp_pos + window:,}")

# Sort by p-value
region_snps = region_snps.sort_values('p_value')

print(f"\nTop 10 SNPs in region:")
print(region_snps[['snp_id', 'pos', 'beta', 'p_value', 'q_value']].head(10).to_string(index=False))

# Load VCF and extract genotypes for scaffold_10 SNPs
print("\n[2/6] EXTRACTING GENOTYPES FROM VCF")
print("-" * 90)

vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'

# Load phenotype data for sample ordering
with gzip.open('./results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    brightness = pd.read_csv(f)

with gzip.open('./input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz', 'rt') as f:
    master = pd.read_csv(f, sep='\t')

master['dbvpg_num'] = master['db_strain_id'].str.extract(r'DBVPG[:\s]+(\d+)')[0]
master['vcf_sample'] = 'DBVPG_' + master['dbvpg_num'].astype(str)

brightness_vcf = brightness.merge(master[['filename', 'vcf_sample']], on='filename', how='left')

# Get VCF samples
vcf_samples = []
with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            vcf_samples = line.strip().split('\t')[9:]
            break

matched = brightness_vcf[brightness_vcf['vcf_sample'].isin(vcf_samples)].copy()
pheno_dict = dict(zip(matched['vcf_sample'], matched['Median_ColorLab_L*Mean']))
sample_idx_map = {sample: idx for idx, sample in enumerate(vcf_samples)}
matched_indices = [sample_idx_map[s] for s in matched['vcf_sample'].values]

print(f"Samples with phenotype: {len(matched)}")

# Extract genotypes for region SNPs
print(f"\nExtracting genotypes for {len(region_snps)} SNPs...")

scaffold10_positions = set(region_snps['pos'].values)
genotypes_dict = {}  # pos -> array of genotypes

with gzip.open(vcf_path, 'rt') as f:
    snp_count = 0
    found_count = 0

    for line in f:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        snp_count += 1

        # Only process scaffold_10 SNPs in region
        if chrom != 'scaffold_10' or pos not in scaffold10_positions:
            continue

        found_count += 1

        # Parse genotypes for matched samples
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
                geno = np.nan

            genotypes.append(geno)

        genotypes_dict[pos] = np.array(genotypes, dtype=float)

        if found_count % 1000 == 0:
            print(f"  Processed {snp_count:,} SNPs, found {found_count} in region")

print(f"Found {found_count} SNPs with genotype data in region")

# Compute LD matrix
print("\n[3/6] COMPUTING LINKAGE DISEQUILIBRIUM")
print("-" * 90)

region_snps['found_genotypes'] = region_snps['pos'].isin(genotypes_dict.keys())
region_snps_with_geno = region_snps[region_snps['found_genotypes']].copy()

print(f"SNPs with genotype data: {len(region_snps_with_geno)}")

if len(region_snps_with_geno) < 2:
    print("ERROR: Need at least 2 SNPs with genotypes for LD computation")
    exit(1)

# Focus on top SNPs for LD computation (top 50)
top_snps_for_ld = region_snps_with_geno.head(50).copy()

print(f"Computing LD for top {len(top_snps_for_ld)} SNPs...")

# Create genotype matrix
positions = top_snps_for_ld['pos'].values
genotype_matrix = np.array([genotypes_dict[pos] for pos in positions]).T  # samples × SNPs

# Compute LD (r-squared)
n_snps = genotype_matrix.shape[1]
ld_matrix = np.zeros((n_snps, n_snps))

for i in range(n_snps):
    for j in range(i, n_snps):
        x = genotype_matrix[:, i]
        y = genotype_matrix[:, j]

        # Filter to non-missing
        valid = ~(np.isnan(x) | np.isnan(y))

        if valid.sum() > 5:
            r, _ = stats.pearsonr(x[valid], y[valid])
            r2 = r ** 2
        else:
            r2 = np.nan

        ld_matrix[i, j] = r2
        ld_matrix[j, i] = r2

print(f"LD matrix computed ({n_snps} × {n_snps})")

# Find SNPs in strong LD with top hit
top_snp_idx = 0  # First SNP is most significant
top_snp_pos = positions[top_snp_idx]
top_snp_ld = ld_matrix[top_snp_idx, :]

print(f"\nLD with top SNP (scaffold_10:{top_snp_pos}):")
ld_df = pd.DataFrame({
    'pos': positions,
    'ld_r2': top_snp_ld
})
ld_df = ld_df.sort_values('ld_r2', ascending=False)
print(ld_df.head(15).to_string(index=False))

# Compute Bayesian credible set
print("\n[4/6] BAYESIAN CREDIBLE SET ANALYSIS")
print("-" * 90)

# Use simple Bayesian approach: posterior probability for each SNP
# Assume equal prior for all SNPs
p_values = region_snps_with_geno['p_value'].values

# Convert to Z-scores (for Bayesian calculation)
z_scores = np.sqrt(-2 * np.log(p_values))

# Calculate posterior probabilities (simple version)
# Posterior ~ likelihood × prior
# Likelihood proportional to exp(-0.5 * z²) / sqrt(2π)
likelihoods = np.exp(-0.5 * z_scores**2)

# Posterior probabilities
posterior_probs = likelihoods / likelihoods.sum()

# Sort by posterior probability
credible_set_df = pd.DataFrame({
    'snp_id': region_snps_with_geno['snp_id'].values,
    'pos': region_snps_with_geno['pos'].values,
    'p_value': region_snps_with_geno['p_value'].values,
    'q_value': region_snps_with_geno['q_value'].values,
    'beta': region_snps_with_geno['beta'].values,
    'posterior_prob': posterior_probs
})

credible_set_df = credible_set_df.sort_values('posterior_prob', ascending=False)

# Find 95% credible set
cumsum_prob = np.cumsum(credible_set_df['posterior_prob'].values)
cs_95_idx = np.where(cumsum_prob >= 0.95)[0][0] + 1

print(f"95% Credible Set: {cs_95_idx} SNPs")
print(f"\nTop SNPs in credible set:")
print(credible_set_df.head(cs_95_idx)[['snp_id', 'pos', 'posterior_prob', 'p_value', 'beta']].to_string(index=False))

# Save results
print("\n[5/6] SAVING FINE-MAPPING RESULTS")
print("-" * 90)

os.makedirs('results/gwas', exist_ok=True)

# Save full LD matrix
ld_full_df = pd.DataFrame(
    ld_matrix,
    columns=positions,
    index=positions
)
ld_full_df.to_csv('results/gwas/finemapping_ld_matrix.csv')
print(f"✓ Saved: finemapping_ld_matrix.csv")

# Save credible set
credible_set_df.to_csv('results/gwas/finemapping_credible_set.csv', index=False)
print(f"✓ Saved: finemapping_credible_set.csv")

# Save summary
finemapping_summary = {
    'locus': 'scaffold_10',
    'top_snp': f'scaffold_10:{top_snp_pos}',
    'window_size_kb': window / 1000,
    'total_snps_in_window': len(region_snps),
    'snps_with_genotypes': len(region_snps_with_geno),
    'snps_in_95_credible_set': int(cs_95_idx),
    'cumulative_posterior_prob_cs95': float(cumsum_prob[cs_95_idx-1])
}

import json
with open('results/gwas/finemapping_summary.json', 'w') as f:
    json.dump(finemapping_summary, f, indent=2)

print(f"✓ Saved: finemapping_summary.json")

# Print final summary
print("\n[6/6] FINE-MAPPING SUMMARY")
print("-" * 90)

print(f"""
SCAFFOLD 10 FINE-MAPPING RESULTS:

Locus: scaffold_10 (±{window/1000:.0f}kb around position {top_snp_pos:,})
Total SNPs tested: {len(gwas_df):,}
SNPs in region: {len(region_snps):,}
SNPs with genotypes: {len(region_snps_with_geno):,}

LD Analysis (with top SNP scaffold_10:{top_snp_pos}):
  SNPs in strong LD (r²>0.8): {(top_snp_ld > 0.8).sum()}
  SNPs in moderate LD (r²>0.5): {(top_snp_ld > 0.5).sum()}
  SNPs with any LD (r²>0.1): {(top_snp_ld > 0.1).sum()}

Bayesian Credible Set Analysis:
  95% credible set size: {cs_95_idx} SNPs
  Posterior prob range: {credible_set_df['posterior_prob'].iloc[0]:.4f} - {credible_set_df['posterior_prob'].iloc[cs_95_idx-1]:.4f}
  Top SNP posterior prob: {credible_set_df['posterior_prob'].iloc[0]:.4f}

Key Finding:
  The {cs_95_idx} SNPs in the 95% credible set are equally likely to be causal
  Focus validation efforts on these top candidates
""")

print("="*90)
print("✓ FINE-MAPPING COMPLETE")
print("="*90)

print(f"\nTop 10 candidate causal SNPs (sorted by posterior probability):")
for idx, row in credible_set_df.head(10).iterrows():
    print(f"  {row['snp_id']:30s}  Posterior: {row['posterior_prob']:.4f}  p={row['p_value']:.2e}  β={row['beta']:.4f}")
