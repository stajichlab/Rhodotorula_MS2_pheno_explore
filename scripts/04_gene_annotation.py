#!/usr/bin/env python3
"""
GWAS Step 4: Gene Annotation
Maps significant SNPs to genes and annotates with functional information
"""

import pandas as pd
import gzip
import os
from collections import defaultdict

print("="*80)
print("STEP 4: GENE ANNOTATION FOR SIGNIFICANT SNPs")
print("="*80)

# Load gene information
print("\n[1/4] LOADING GENE ANNOTATION")
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

        # Parse gene ID and name
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
            'strand': '+' if start < end else '-'
        })

genes_df = pd.DataFrame(genes)
print(f"Loaded {len(genes_df)} genes from reference genome")
print(f"Scaffolds: {genes_df['scaffold'].nunique()}")

# Create index for fast lookup
genes_by_scaffold = defaultdict(list)
for _, row in genes_df.iterrows():
    genes_by_scaffold[row['scaffold']].append(row)

# Sort each scaffold's genes by start position
for scaffold in genes_by_scaffold:
    genes_by_scaffold[scaffold].sort(key=lambda x: x['start'])

print(f"Indexed genes by scaffold")

# Load Pfam annotations for carotenoid pathway identification
print("\n[2/4] LOADING PROTEIN ANNOTATIONS")
print("-" * 80)

carotenoid_genes = {
    'PSY': ['PF03659'],  # Phytoene synthase
    'PDS': ['PF06857'],  # Phytoene desaturase
    'ZDS': ['PF06856'],  # Zeta-carotene desaturase
    'Lycopene_cyc': ['PF05155'],  # Lycopene cyclase
    'BCH': ['PF04551'],  # Beta-carotene hydroxylase
    'GGPS': ['PF00348'],  # Geranylgeranyl pyrophosphate synthase
}

print(f"Carotenoid pathway genes: {list(carotenoid_genes.keys())}")

# Load Pfam annotations
pfam_file = '/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Rodeo/tables/All_Taxa/pfam.csv.gz'

pfam_by_gene = defaultdict(list)
pfam_count = 0

if os.path.exists(pfam_file):
    print(f"\nLoading Pfam annotations...")
    with gzip.open(pfam_file, 'rt') as f:
        # Skip header
        header = next(f)
        for line in f:
            fields = line.strip().split(',')
            if len(fields) < 3:
                continue

            # Format: sample, gene_id, pfam_id, pfam_name, ...
            sample = fields[0]
            gene_id = fields[1]
            pfam_id = fields[2]
            pfam_name = fields[3] if len(fields) > 3 else ""

            # Filter for Y2510 reference strain
            if 'Y2510' in sample or 'NRRL' in sample or sample == 'Rhodotorula_mucilaginosa_NRRL_Y-2510':
                pfam_by_gene[gene_id].append({
                    'pfam_id': pfam_id,
                    'pfam_name': pfam_name
                })
                pfam_count += 1

    print(f"  Loaded {pfam_count} Pfam annotations")
else:
    print(f"  WARNING: Pfam file not found at {pfam_file}")
    print(f"  Proceeding without protein domain annotation")

# Load significant SNPs
print("\n[3/4] MAPPING SNPs TO GENES")
print("-" * 80)

# Load Bonferroni-significant SNPs
snp_file = 'results/gwas/gwas_sig_bonf.csv'
snps = pd.read_csv(snp_file)

print(f"Loaded {len(snps)} Bonferroni-significant SNPs")

# Map SNPs to genes
annotated_snps = []

for idx, snp in snps.iterrows():
    scaffold = snp['chr']
    pos = int(snp['pos'])

    # Find genes overlapping this position on the same scaffold
    overlapping_genes = []

    if scaffold in genes_by_scaffold:
        for gene in genes_by_scaffold[scaffold]:
            if gene['start'] <= pos <= gene['end']:
                overlapping_genes.append(gene)
            elif gene['start'] > pos + 1000:  # Stop searching after 1000bp past SNP
                break

    # Find nearest genes (within 2000bp) if no overlap
    if not overlapping_genes and scaffold in genes_by_scaffold:
        nearby_genes = []
        for gene in genes_by_scaffold[scaffold]:
            if abs(gene['start'] - pos) <= 2000:
                nearby_genes.append((gene, abs(gene['start'] - pos)))
            elif abs(gene['end'] - pos) <= 2000:
                nearby_genes.append((gene, abs(gene['end'] - pos)))

        nearby_genes.sort(key=lambda x: x[1])
        if nearby_genes:
            nearest = nearby_genes[0][0]
            distance = nearby_genes[0][1]
            overlapping_genes = [nearest]
            distance_label = f"promoter ({distance}bp)"
        else:
            distance_label = "intergenic"
    else:
        distance_label = "coding" if overlapping_genes else "intergenic"

    # Annotate with pathway information
    pathway_genes = []
    for gene in overlapping_genes:
        gene_id = gene['gene_id']

        # Check if this is a carotenoid pathway gene
        is_carotenoid = False
        pathway_name = ""

        if gene_id in pfam_by_gene:
            for pfam in pfam_by_gene[gene_id]:
                pfam_id = pfam['pfam_id']
                pfam_name = pfam['pfam_name']

                for pathway, pfams in carotenoid_genes.items():
                    if pfam_id in pfams:
                        is_carotenoid = True
                        pathway_name = pathway
                        break

        # Also check gene name
        if gene['gene_name']:
            for pathway in carotenoid_genes:
                if pathway.lower() in gene['gene_name'].lower():
                    is_carotenoid = True
                    pathway_name = pathway

        pathway_genes.append({
            'gene_id': gene_id,
            'gene_name': gene['gene_name'],
            'is_carotenoid': is_carotenoid,
            'pathway_gene': pathway_name,
            'location': distance_label
        })

    # Create annotated SNP record
    if overlapping_genes:
        for pg in pathway_genes:
            annotated_snps.append({
                'snp_id': snp['snp_id'],
                'chr': snp['chr'],
                'pos': snp['pos'],
                'beta': snp['beta'],
                'p_value': snp['p_value'],
                'gene_id': pg['gene_id'],
                'gene_name': pg['gene_name'],
                'location': pg['location'],
                'is_carotenoid': pg['is_carotenoid'],
                'pathway_gene': pg['pathway_gene']
            })
    else:
        annotated_snps.append({
            'snp_id': snp['snp_id'],
            'chr': snp['chr'],
            'pos': snp['pos'],
            'beta': snp['beta'],
            'p_value': snp['p_value'],
            'gene_id': None,
            'gene_name': None,
            'location': 'intergenic',
            'is_carotenoid': False,
            'pathway_gene': None
        })

annotated_df = pd.DataFrame(annotated_snps)

print(f"\nAnnotation summary:")
print(f"  SNPs in coding regions: {(annotated_df['location'] == 'coding').sum()}")
print(f"  SNPs in carotenoid pathway genes: {annotated_df['is_carotenoid'].sum()}")
print(f"  Unique genes hit: {annotated_df['gene_id'].nunique()}")

# Save results
print("\n[4/4] SAVING ANNOTATED RESULTS")
print("-" * 80)

os.makedirs('results/gwas', exist_ok=True)

annotated_df.to_csv('results/gwas/gwas_sig_bonf_annotated.csv', index=False)
print(f"✓ Saved: gwas_sig_bonf_annotated.csv")

# Create summary of carotenoid pathway hits
carotenoid_hits = annotated_df[annotated_df['is_carotenoid']]
if len(carotenoid_hits) > 0:
    carotenoid_hits.to_csv('results/gwas/gwas_carotenoid_hits.csv', index=False)
    print(f"✓ Saved: gwas_carotenoid_hits.csv ({len(carotenoid_hits)} hits)")

# Create summary table
print("\n" + "="*80)
print("ANNOTATED SIGNIFICANT SNPs (n={})".format(len(annotated_df)))
print("="*80)

summary_cols = ['snp_id', 'chr', 'beta', 'p_value', 'gene_name', 'pathway_gene', 'location']
print(annotated_df[summary_cols].to_string(index=False))

if len(carotenoid_hits) > 0:
    print("\n" + "="*80)
    print("CAROTENOID PATHWAY HITS (n={})".format(len(carotenoid_hits)))
    print("="*80)
    print(carotenoid_hits[['snp_id', 'chr', 'beta', 'gene_name', 'pathway_gene']].to_string(index=False))
else:
    print("\n" + "="*80)
    print("⚠ NO DIRECT CAROTENOID PATHWAY HITS IN SIGNIFICANT SNPs")
    print("="*80)
    print("Interpretation: SNPs may be in regulatory regions or regulatory genes")
    print("Consider FDR-significant SNPs for additional candidates")

print("\n✓ GENE ANNOTATION COMPLETE")
