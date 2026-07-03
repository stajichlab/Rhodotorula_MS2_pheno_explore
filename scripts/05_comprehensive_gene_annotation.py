#!/usr/bin/env python3
"""
GWAS Step 5: Comprehensive Gene Annotation & Analysis
Maps all significant SNPs to genes and creates detailed gene function analysis
"""

import pandas as pd
import gzip
import os
from collections import defaultdict

print("="*80)
print("STEP 5: COMPREHENSIVE GENE ANNOTATION & ANALYSIS")
print("="*80)

# Load gene information
print("\n[1/5] LOADING GENE ANNOTATION")
print("-" * 80)

gene_file = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.genes.csv'

genes = []
with open(gene_file, 'r') as f:
    for line in f:
        fields = line.strip().split(',')
        scaffold = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        info = fields[3]

        gene_id = info.split(';')[0]
        gene_name = None
        if 'Name=' in info:
            gene_name = info.split('Name=')[1]

        genes.append({
            'scaffold': scaffold,
            'start': start,
            'end': end,
            'gene_id': gene_id,
            'gene_name': gene_name,
            'length': abs(end - start)
        })

genes_df = pd.DataFrame(genes)

# Create index for fast lookup
genes_by_scaffold = defaultdict(list)
for _, row in genes_df.iterrows():
    genes_by_scaffold[row['scaffold']].append(row)

for scaffold in genes_by_scaffold:
    genes_by_scaffold[scaffold].sort(key=lambda x: x['start'])

print(f"Loaded {len(genes_df)} genes")

# Known pathway genes
pathway_genes_db = {
    'PSY': ['phytoene synthase', 'phytoene'],
    'PDS': ['phytoene desaturase', 'desaturase'],
    'ZDS': ['zeta-carotene desaturase', 'zeta-carotene'],
    'Lycopene_cyc': ['lycopene cyclase', 'cyclase'],
    'BCH': ['beta-carotene hydroxylase', 'hydroxylase'],
    'GGPS': ['geranylgeranyl pyrophosphate synthase', 'prenyltransferase'],
    'ABCA': ['ABC transporter'],
    'SQS': ['squalene synthase'],
    'ERG': ['ergosterol biosynthesis', 'sterol'],
}

print(f"Pathway genes tracked: {list(pathway_genes_db.keys())}")

# Function to map SNP to genes
def map_snp_to_genes(scaffold, pos, genes_by_scaffold):
    overlapping = []
    nearby = []

    if scaffold in genes_by_scaffold:
        for gene in genes_by_scaffold[scaffold]:
            if gene['start'] <= pos <= gene['end']:
                overlapping.append((gene, 'coding'))
            else:
                # Find nearby genes (promoter region: 2000bp upstream)
                if gene['end'] < pos < gene['end'] + 2000:
                    distance = pos - gene['end']
                    nearby.append((gene, f'downstream ({distance}bp)', distance))
                elif gene['start'] - 2000 < pos < gene['start']:
                    distance = gene['start'] - pos
                    nearby.append((gene, f'promoter ({distance}bp)', distance))

    # Return overlapping genes if found, otherwise nearby
    if overlapping:
        return overlapping
    elif nearby:
        nearby.sort(key=lambda x: x[2])
        return [nearby[0][:2]]
    else:
        return []

# Load and annotate Bonferroni SNPs
print("\n[2/5] ANNOTATING BONFERRONI-SIGNIFICANT SNPs (n=12)")
print("-" * 80)

bonf_snps = pd.read_csv('results/gwas/gwas_sig_bonf.csv')

bonf_annotated = []
for idx, snp in bonf_snps.iterrows():
    genes_hit = map_snp_to_genes(snp['chr'], int(snp['pos']), genes_by_scaffold)

    if genes_hit:
        for gene, location in genes_hit:
            bonf_annotated.append({
                'significance': 'Bonferroni',
                'snp_id': snp['snp_id'],
                'chr': snp['chr'],
                'pos': snp['pos'],
                'beta': snp['beta'],
                'p_value': snp['p_value'],
                'q_value': snp['q_value'],
                'gene_id': gene['gene_id'],
                'gene_name': gene['gene_name'],
                'gene_length': gene['length'],
                'location': location
            })
    else:
        bonf_annotated.append({
            'significance': 'Bonferroni',
            'snp_id': snp['snp_id'],
            'chr': snp['chr'],
            'pos': snp['pos'],
            'beta': snp['beta'],
            'p_value': snp['p_value'],
            'q_value': snp['q_value'],
            'gene_id': None,
            'gene_name': None,
            'gene_length': None,
            'location': 'intergenic'
        })

bonf_df = pd.DataFrame(bonf_annotated)

print(f"Bonferroni SNPs annotated:")
print(f"  In coding regions: {(bonf_df['location'] == 'coding').sum()}")
print(f"  In promoter regions: {bonf_df['location'].str.contains('promoter', na=False).sum()}")
print(f"  Intergenic: {(bonf_df['location'] == 'intergenic').sum()}")
print(f"  Unique genes: {bonf_df['gene_id'].nunique()}")

# Load and annotate FDR SNPs (sample of top 100)
print("\n[3/5] ANNOTATING FDR-SIGNIFICANT SNPs (sample of top 100/20386)")
print("-" * 80)

fdr_snps = pd.read_csv('results/gwas/gwas_sig_fdr.csv', nrows=100)

fdr_annotated = []
for idx, snp in fdr_snps.iterrows():
    genes_hit = map_snp_to_genes(snp['chr'], int(snp['pos']), genes_by_scaffold)

    if genes_hit:
        for gene, location in genes_hit:
            fdr_annotated.append({
                'significance': 'FDR',
                'snp_id': snp['snp_id'],
                'chr': snp['chr'],
                'pos': snp['pos'],
                'beta': snp['beta'],
                'p_value': snp['p_value'],
                'q_value': snp['q_value'],
                'gene_id': gene['gene_id'],
                'gene_name': gene['gene_name'],
                'gene_length': gene['length'],
                'location': location
            })
    else:
        fdr_annotated.append({
            'significance': 'FDR',
            'snp_id': snp['snp_id'],
            'chr': snp['chr'],
            'pos': snp['pos'],
            'beta': snp['beta'],
            'p_value': snp['p_value'],
            'q_value': snp['q_value'],
            'gene_id': None,
            'gene_name': None,
            'gene_length': None,
            'location': 'intergenic'
        })

fdr_df = pd.DataFrame(fdr_annotated)

print(f"FDR SNPs annotated (top 100):")
print(f"  In coding regions: {(fdr_df['location'] == 'coding').sum()}")
print(f"  In promoter regions: {fdr_df['location'].str.contains('promoter', na=False).sum()}")
print(f"  Intergenic: {(fdr_df['location'] == 'intergenic').sum()}")
print(f"  Unique genes: {fdr_df['gene_id'].nunique()}")

# Combine all annotated SNPs
all_annotated = pd.concat([bonf_df, fdr_df], ignore_index=True)

# Save results
print("\n[4/5] SAVING ANNOTATED RESULTS")
print("-" * 80)

os.makedirs('results/gwas', exist_ok=True)

bonf_df.to_csv('results/gwas/gwas_bonf_annotated.csv', index=False)
fdr_df.to_csv('results/gwas/gwas_fdr_annotated_top100.csv', index=False)
all_annotated.to_csv('results/gwas/gwas_all_annotated.csv', index=False)

print(f"✓ Saved: gwas_bonf_annotated.csv (12 SNPs)")
print(f"✓ Saved: gwas_fdr_annotated_top100.csv (100 SNPs)")
print(f"✓ Saved: gwas_all_annotated.csv (combined)")

# Create gene summary
print("\n[5/5] GENE SUMMARY & ANALYSIS")
print("-" * 80)

# Get unique genes hit by Bonferroni SNPs
bonf_genes = bonf_df[bonf_df['gene_id'].notna()].groupby('gene_id').agg({
    'gene_name': 'first',
    'snp_id': 'count',
    'beta': 'mean',
    'p_value': 'min'
}).rename(columns={'snp_id': 'n_snps'}).reset_index()

bonf_genes = bonf_genes.sort_values('p_value')

print("\nGENES HIT BY BONFERRONI-SIGNIFICANT SNPs (n={}):".format(len(bonf_genes)))
print(bonf_genes[['gene_id', 'gene_name', 'n_snps', 'beta']].to_string(index=False))

# Get unique genes hit by FDR SNPs (top 100)
fdr_genes = fdr_df[fdr_df['gene_id'].notna()].groupby('gene_id').agg({
    'gene_name': 'first',
    'snp_id': 'count',
    'beta': 'mean',
    'p_value': 'min'
}).rename(columns={'snp_id': 'n_snps'}).reset_index()

fdr_genes = fdr_genes.sort_values('p_value')

print("\nGENES HIT BY FDR-SIGNIFICANT SNPs (top 100, n={}):".format(len(fdr_genes)))
print(fdr_genes[['gene_id', 'gene_name', 'n_snps', 'beta']].head(20).to_string(index=False))

# Summary statistics
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

print(f"\nBonferroni SNPs:")
print(f"  Total SNPs: {len(bonf_df)}")
print(f"  Unique genes: {bonf_df['gene_id'].nunique()}")
print(f"  Avg effect size: {bonf_df['beta'].mean():.4f}")

print(f"\nFDR SNPs (top 100):")
print(f"  Total SNPs: {len(fdr_df)}")
print(f"  Unique genes: {fdr_df['gene_id'].nunique()}")
print(f"  Avg effect size: {fdr_df['beta'].mean():.4f}")

# Location distribution
print(f"\nSNP location distribution (Bonferroni):")
print(bonf_df['location'].value_counts())

print("\n✓ COMPREHENSIVE GENE ANNOTATION COMPLETE")
