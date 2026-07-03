#!/usr/bin/env python3
"""
GWAS Step 3: Generate Publication-Quality Visualizations
Creates Manhattan plots, Q-Q plots, and summary figures
"""

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("STEP 3: GENERATING GWAS VISUALIZATIONS")
print("="*80)

# Load results
print("\nLoading GWAS results...")
results_df = pd.read_csv('results/gwas/gwas_all_results.csv')

with open('results/gwas/gwas_summary.json', 'r') as f:
    summary = json.load(f)

bonf_thresh = summary['bonf_threshold']

print(f"  Loaded {len(results_df):,} SNPs")
print(f"  Bonferroni threshold: {bonf_thresh:.2e}")

# Create 4-panel summary figure
print("\nGenerating 4-panel summary figure...")
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel 1: Manhattan Plot
ax = axes[0, 0]
sorted_df = results_df.sort_values('neg_log10_p', ascending=False)
ax.scatter(range(len(sorted_df)), sorted_df['neg_log10_p'], alpha=0.5, s=20, color='steelblue')
ax.axhline(-np.log10(bonf_thresh), color='red', linestyle='--', linewidth=2, label=f'Bonferroni')
ax.set_xlabel('SNP rank', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Manhattan Plot', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(alpha=0.3)

# Panel 2: Q-Q Plot
ax = axes[0, 1]
pvals_sorted = np.sort(results_df['p_value'].values)
expected = -np.log10(np.linspace(1.0/len(pvals_sorted), 1, len(pvals_sorted)))
observed = -np.log10(pvals_sorted)
ax.scatter(expected, observed, alpha=0.5, s=20, color='steelblue')
max_val = max(expected.max(), observed.max())
ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Expected')
ax.set_xlabel('Expected -log10(p)', fontsize=12)
ax.set_ylabel('Observed -log10(p)', fontsize=12)
ax.set_title('Q-Q Plot', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(alpha=0.3)

# Panel 3: Effect Size vs. Significance
ax = axes[1, 0]
colors = ['red' if sig else 'steelblue' for sig in results_df['sig_bonf']]
ax.scatter(results_df['beta'], results_df['neg_log10_p'], c=colors, alpha=0.5, s=20)
ax.axhline(-np.log10(bonf_thresh), color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Effect Size (β)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Effect Size vs. Significance', fontsize=14, fontweight='bold')
ax.grid(alpha=0.3)

# Panel 4: P-value Histogram
ax = axes[1, 1]
ax.hist(results_df['neg_log10_p'], bins=100, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(-np.log10(bonf_thresh), color='red', linestyle='--', linewidth=2, label='Bonferroni')
ax.set_xlabel('-log10(p-value)', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('P-value Distribution', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('results/gwas/gwas_summary.png', dpi=200, bbox_inches='tight')
print("✓ Saved: gwas_summary.png")

# Create Manhattan plot by chromosome
print("\nGenerating chromosome-wise Manhattan plot...")
fig, ax = plt.subplots(figsize=(18, 6))

# Group by chromosome and plot
chrom_groups = results_df.groupby('chr')
colors_list = plt.cm.Set3(np.linspace(0, 1, len(chrom_groups)))

x_pos = 0
chrom_boundaries = []

for color_idx, (chrom, group) in enumerate(chrom_groups):
    group = group.sort_values('pos')
    positions = np.arange(len(group)) + x_pos
    ax.scatter(positions, group['neg_log10_p'], alpha=0.6, s=15,
               color=colors_list[color_idx], label=chrom)
    chrom_boundaries.append(x_pos + len(group) / 2)
    x_pos += len(group) + 100

ax.axhline(-np.log10(bonf_thresh), color='red', linestyle='--', linewidth=2, label='Bonferroni')
ax.set_xlabel('Chromosome', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Manhattan Plot by Chromosome', fontsize=14, fontweight='bold')
ax.set_xticks(chrom_boundaries)
ax.set_xticklabels([c for c in chrom_groups.groups.keys()], rotation=45)
ax.grid(alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('results/gwas/gwas_manhattan_by_chr.png', dpi=200, bbox_inches='tight')
print("✓ Saved: gwas_manhattan_by_chr.png")

# Effect size distribution
print("\nGenerating effect size distribution plots...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Beta distribution
ax = axes[0]
ax.hist(results_df['beta'], bins=100, alpha=0.7, color='steelblue', edgecolor='black')
ax.axvline(results_df['beta'].mean(), color='red', linestyle='--', linewidth=2, label=f'Mean={results_df["beta"].mean():.3f}')
ax.axvline(results_df['beta'].median(), color='green', linestyle='--', linewidth=2, label=f'Median={results_df["beta"].median():.3f}')
ax.set_xlabel('Effect Size (β)', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('Distribution of Effect Sizes', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(alpha=0.3)

# Volcano plot
ax = axes[1]
colors = ['red' if sig else 'steelblue' for sig in results_df['sig_bonf']]
ax.scatter(results_df['beta'], -np.log10(results_df['p_value']), c=colors, alpha=0.5, s=20)
ax.axhline(-np.log10(bonf_thresh), color='gray', linestyle='--', alpha=0.5, label='Bonferroni')
ax.axvline(0, color='gray', linestyle='-', alpha=0.3)
ax.set_xlabel('Effect Size (β)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title('Volcano Plot', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('results/gwas/gwas_effect_size.png', dpi=200, bbox_inches='tight')
print("✓ Saved: gwas_effect_size.png")

print("\n" + "="*80)
print("✓ VISUALIZATION COMPLETE")
print("="*80)
print("\nGenerated figures:")
print("  - gwas_summary.png (4-panel: Manhattan, Q-Q, effects, histogram)")
print("  - gwas_manhattan_by_chr.png (chromosome-wise Manhattan plot)")
print("  - gwas_effect_size.png (effect size distribution and volcano plot)")
