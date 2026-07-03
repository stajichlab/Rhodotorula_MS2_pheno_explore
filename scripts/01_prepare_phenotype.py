#!/usr/bin/env python3
"""
GWAS Step 1: Prepare Phenotype File for GWAS Analysis
Matches brightness measurements to VCF sample names
"""

import pandas as pd
import gzip
import os

print("="*80)
print("STEP 1: PREPARE PHENOTYPE FILE FOR GWAS")
print("="*80)

# Load data
print("\nLoading brightness measurements...")
with gzip.open('./results/phase1/phase1_phenotype_data.csv.gz', 'rt') as f:
    brightness = pd.read_csv(f)

print(f"  Loaded {len(brightness)} samples")
print(f"  Brightness range: {brightness['Median_ColorLab_L*Mean'].min():.2f} - {brightness['Median_ColorLab_L*Mean'].max():.2f} L*")

# Load master metadata
print("\nLoading master metadata...")
with gzip.open('./input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz', 'rt') as f:
    master = pd.read_csv(f, sep='\t')

print(f"  Loaded {len(master)} samples")

# Map to VCF sample names
print("\nMapping to VCF sample IDs...")
master['dbvpg_num'] = master['db_strain_id'].str.extract(r'DBVPG[:\s]+(\d+)')[0]
master['vcf_sample'] = 'DBVPG_' + master['dbvpg_num'].astype(str)

brightness_vcf = brightness.merge(master[['filename', 'vcf_sample']], on='filename', how='left')

# Get VCF sample list
print("\nReading VCF sample names...")
vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'
vcf_samples = []
with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'):
            vcf_samples = line.strip().split('\t')[9:]
            break

print(f"  Found {len(vcf_samples)} VCF samples")

# Filter matched samples
matched = brightness_vcf[brightness_vcf['vcf_sample'].isin(vcf_samples)].copy()

print(f"\nMatched samples: {len(matched)}")
print(f"  Brightness: {matched['Median_ColorLab_L*Mean'].mean():.2f} ± {matched['Median_ColorLab_L*Mean'].std():.2f} L*")

# Create PLINK phenotype file (FID IID Phenotype)
print("\nCreating phenotype file for PLINK...")
pheno = matched[['vcf_sample', 'Median_ColorLab_L*Mean']].copy()
pheno.insert(0, 'FID', '0')  # PLINK uses FID=0 for const-fid
pheno.columns = ['FID', 'IID', 'Brightness']

os.makedirs('results/gwas', exist_ok=True)
pheno.to_csv('results/gwas/brightness_pheno.txt', sep='\t', index=False, header=True)

print(f"✓ Phenotype file created: results/gwas/brightness_pheno.txt")
print(f"  Samples: {len(pheno)}")
print(f"\nFirst 5 rows:")
print(pheno.head())
