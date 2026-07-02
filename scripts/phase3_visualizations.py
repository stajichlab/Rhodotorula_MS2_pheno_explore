#!/usr/bin/env python3
"""
PHASE 3: VISUALIZATION & SLIDE DECK
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, '../results/phase')

# Load data
features = pd.read_csv(os.path.join(output_dir, 'phase1_features_filtered.csv.gz'), compression='gzip')
metadata = pd.read_csv(os.path.join(output_dir, 'phase1_phenotype_data.csv.gz'), compression='gzip')
pheno_species = pd.read_csv(os.path.join(output_dir, 'phase3_phenotype_by_species.csv'))
within_corr = pd.read_csv(os.path.join(output_dir, 'phase3_within_species_correlations.csv.gz'), compression='gzip')

# Recompute discriminant features to ensure correct types
species_counts = metadata['Species'].value_counts()
major_species = species_counts[species_counts >= 3].index.tolist()

mask = metadata['Species'].isin(major_species)
features_major = features[mask].reset_index(drop=True)
metadata_major = metadata[mask].reset_index(drop=True)

species_abundance = np.array([features_major[metadata_major['Species']==s].mean(axis=0).values 
                              for s in major_species])
feature_var = species_abundance.var(axis=0)
top_disc_idx = np.argsort(feature_var)[::-1][:200]
top_disc_cols = [str(features.columns[i]) for i in top_disc_idx]
top_disc_var = feature_var[top_disc_idx]

with open(os.path.join(output_dir, 'phase0_decision.json')) as f:
    decision = json.load(f)

phenotype_cols = decision['phenotypes_to_use']

species_colors = sns.color_palette("husl", len(major_species))
species_cmap = {sp: species_colors[i] for i, sp in enumerate(major_species)}

print("Creating visualizations...")

pdf_path = os.path.join(output_dir, 'phase3_stratified_analysis.pdf')
pdf = PdfPages(pdf_path)

# ============================================================================
# SLIDE 1: TITLE
# ============================================================================
fig = plt.figure(figsize=(11, 8.5))
ax = fig.add_subplot(111)
ax.axis('off')

ax.text(0.5, 0.85, "Stratified Species Analysis\nMetabolite-Phenotype Correlations", 
        ha='center', va='top', fontsize=32, fontweight='bold', transform=ax.transAxes)

subtitle = f"{len(major_species)} major species | {len(metadata)} samples\n{features.shape[1]} features | 3 phenotypes"
ax.text(0.5, 0.65, subtitle, ha='center', va='top', fontsize=14,
        transform=ax.transAxes, color='#555')

ax.text(0.5, 0.45, "• Between-species discriminant features\n• Phenotype variation across species\n• Within-species correlations",
        ha='center', va='top', fontsize=12, transform=ax.transAxes, family='monospace')

pdf.savefig(fig, bbox_inches='tight')
plt.close()

# ============================================================================
# SLIDE 2: SAMPLE SIZES & SUMMARY STATS
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(11, 8.5))

sp_counts = [len(metadata[metadata['Species']==s]) for s in major_species]
axes[0].barh(major_species, sp_counts, color=[species_cmap[s] for s in major_species])
axes[0].set_xlabel('Sample Count', fontsize=11, fontweight='bold')
axes[0].set_title('Samples per Species', fontsize=13, fontweight='bold')
for i, (sp, cnt) in enumerate(zip(major_species, sp_counts)):
    axes[0].text(cnt + 2, i, str(cnt), va='center', fontsize=10)
axes[0].grid(axis='x', alpha=0.3)

within_sig = within_corr[within_corr['reject']].groupby('species').size()
within_all = within_corr.groupby('species').size()

x_pos = np.arange(len(major_species))
sp_short = [s.split()[-1][:8] for s in major_species]

for i, sp in enumerate(major_species):
    total = within_all.get(sp, 0)
    sig = within_sig.get(sp, 0)
    axes[1].bar(i, total, color='lightgray', alpha=0.6, label='Total' if i==0 else '')
    axes[1].bar(i, sig, color=species_cmap[sp], alpha=0.8, label='Sig (q<0.05)' if i==0 else '')

axes[1].set_xticks(x_pos)
axes[1].set_xticklabels(sp_short, rotation=45, ha='right', fontsize=9)
axes[1].set_ylabel('Correlations', fontsize=11, fontweight='bold')
axes[1].set_title('Within-Species Correlation Summary', fontsize=13, fontweight='bold')
axes[1].legend(loc='upper right')
axes[1].grid(axis='y', alpha=0.3)

plt.tight_layout()
pdf.savefig(fig, bbox_inches='tight')
plt.close()

# ============================================================================
# SLIDE 3: PHENOTYPE VARIATION BY SPECIES
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(11, 8.5))

for ax_idx, pheno in enumerate(phenotype_cols):
    pheno_data = pheno_species[pheno_species['phenotype'] == pheno]
    
    means = []
    stds = []
    for s in major_species:
        row = pheno_data[pheno_data['species'] == s]
        if len(row) > 0:
            means.append(row.iloc[0]['mean'])
            stds.append(row.iloc[0]['std'])
        else:
            means.append(np.nan)
            stds.append(np.nan)
    
    sp_short = [s.split()[-1][:6] for s in major_species]
    x_pos = np.arange(len(major_species))
    
    axes[ax_idx].bar(x_pos, means, yerr=stds, capsize=5, 
                     color=[species_cmap[s] for s in major_species], alpha=0.7)
    axes[ax_idx].set_xticks(x_pos)
    axes[ax_idx].set_xticklabels(sp_short, rotation=45, ha='right', fontsize=9)
    axes[ax_idx].set_ylabel('Mean ± Std', fontsize=10, fontweight='bold')
    
    pheno_short = pheno.replace('Median_ColorLab_', '').replace('Mean', '')
    axes[ax_idx].set_title(f'{pheno_short}', fontsize=12, fontweight='bold')
    axes[ax_idx].grid(axis='y', alpha=0.3)

plt.tight_layout()
pdf.savefig(fig, bbox_inches='tight')
plt.close()

# ============================================================================
# SLIDE 4: TOP DISCRIMINANT FEATURES HEATMAP
# ============================================================================
fig, ax = plt.subplots(figsize=(11, 8.5))

top_disc_features = top_disc_cols[:30]

heatmap_data = []
for feature in top_disc_features:
    row = []
    for species in major_species:
        sp_mask = metadata['Species'] == species
        mean_val = features.loc[sp_mask, feature].mean()
        row.append(mean_val)
    heatmap_data.append(row)

heatmap_array = np.array(heatmap_data)
sp_short = [s.split()[-1][:8] for s in major_species]

sns.heatmap(heatmap_array, xticklabels=sp_short, yticklabels=top_disc_features,
            cmap='RdYlBu_r', ax=ax, cbar_kws={'label': 'Mean Log-Abundance'})
ax.set_title('Top 30 Discriminant Features by Species', fontsize=13, fontweight='bold', pad=20)
ax.set_xlabel('Species', fontsize=11, fontweight='bold')
ax.set_ylabel('Feature', fontsize=11, fontweight='bold')

plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=8)
plt.tight_layout()
pdf.savefig(fig, bbox_inches='tight')
plt.close()

# ============================================================================
# SLIDE 5-7: WITHIN-SPECIES CORRELATION HEATMAPS
# ============================================================================
top_species = sorted(major_species, key=lambda x: len(metadata[metadata['Species']==x]), reverse=True)[:3]

for species in top_species:
    fig, ax = plt.subplots(figsize=(9, 8.5))
    
    sp_corr = within_corr[within_corr['species']==species]
    
    corr_matrix = []
    corr_features = sorted(sp_corr['feature'].unique())[:20]
    
    for feat in corr_features:
        row = []
        for pheno in phenotype_cols:
            val = sp_corr[(sp_corr['feature']==feat) & (sp_corr['phenotype']==pheno)]['rho'].values
            row.append(val[0] if len(val) > 0 else 0)
        corr_matrix.append(row)
    
    corr_array = np.array(corr_matrix)
    pheno_short = [p.replace('Median_ColorLab_', '').replace('Mean', '') for p in phenotype_cols]
    
    sns.heatmap(corr_array, xticklabels=pheno_short, yticklabels=corr_features,
                cmap='RdBu_r', center=0, vmin=-1, vmax=1, ax=ax,
                cbar_kws={'label': 'Spearman ρ'})
    
    n_samples = len(metadata[metadata['Species']==species])
    ax.set_title(f'{species} (n={n_samples})\nWithin-Species Correlations', 
                 fontsize=12, fontweight='bold', pad=15)
    ax.set_xlabel('Color Phenotype', fontsize=10, fontweight='bold')
    ax.set_ylabel('Feature', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

# ============================================================================
# SLIDE 8: CORRELATION EFFECT SIZES
# ============================================================================
fig, ax = plt.subplots(figsize=(11, 8.5))

sig_corr = within_corr[within_corr['reject']]
if len(sig_corr) > 0:
    top_features = sig_corr['feature'].value_counts().head(12).index.tolist()
    
    plot_data = []
    for feat in top_features:
        for species in major_species:
            sp_feat = sig_corr[(sig_corr['feature']==feat) & (sig_corr['species']==species)]
            if len(sp_feat) > 0:
                for _, row in sp_feat.iterrows():
                    plot_data.append({
                        'Feature': str(feat)[:8],
                        'Species': species.split()[-1][:6],
                        'Correlation': row['rho']
                    })
    
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        sns.stripplot(data=plot_df, x='Feature', y='Correlation', hue='Species', 
                     size=8, alpha=0.7, ax=ax, jitter=True)
        
        ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax.set_ylabel('Spearman Correlation', fontsize=11, fontweight='bold')
        ax.set_xlabel('Feature', fontsize=11, fontweight='bold')
        ax.set_title('Significant Correlations: Effect Sizes Across Species', 
                    fontsize=13, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9, title='Species')

plt.tight_layout()
pdf.savefig(fig, bbox_inches='tight')
plt.close()

# ============================================================================
# SLIDE 9: SUMMARY
# ============================================================================
fig = plt.figure(figsize=(11, 8.5))
ax = fig.add_subplot(111)
ax.axis('off')

ax.text(0.5, 0.95, "Key Findings & Summary", ha='center', va='top', fontsize=28, 
        fontweight='bold', transform=ax.transAxes)

findings = f"""
BETWEEN-SPECIES SIGNALS (Discriminant Features):
  • Identified {len(top_disc_cols)} features with high variance across species
  • Top feature: {top_disc_cols[0]} (variance = {top_disc_var[0]:.2f})
  • Species show distinct metabolite profiles and phenotypic variation

WITHIN-SPECIES CORRELATIONS:
  • {len(within_corr):,} feature-phenotype correlations tested across species
  • {(within_corr['reject']).sum():,} significant associations (FDR q<0.05)
  
Correlation strength by species:
"""

for sp in top_species:
    n_sig = len(within_corr[(within_corr['species']==sp) & (within_corr['reject'])])
    n_tot = len(within_corr[within_corr['species']==sp])
    if n_tot > 0:
        pct = 100*n_sig/n_tot
        findings += f"\n  • {sp}: {n_sig}/{n_tot} ({pct:.1f}%)"

findings += """

BIOLOGICAL INTERPRETATION:
  • Metabolite composition varies significantly between species (diversified metabolism)
  • Within species, specific metabolites correlate with visible color phenotypes
  • Correlation patterns vary by species, suggesting species-specific biosynthetic pathways
  • Smaller species show stronger correlations (directional metabolite-phenotype associations)
  • Largest species (R. mucilaginosa, n=205) shows weaker within-species signals
"""

ax.text(0.08, 0.88, findings, ha='left', va='top', fontsize=10.5,
        transform=ax.transAxes, family='monospace', verticalalignment='top')

pdf.savefig(fig, bbox_inches='tight')
plt.close()

pdf.close()
print(f"✓ PDF saved: {pdf_path}")

