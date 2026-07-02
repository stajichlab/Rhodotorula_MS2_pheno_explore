#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# Load data
file_path = '/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/feature_extractMS2/MS2_samples_combine.extended_metadata_with_strain_traits.tsv'
df = pd.read_csv(file_path, sep='\t')

# Filter to samples with phenotype data
df_pheno = df[df['Median_Shape_Area'].notna()].copy()
print(f"Samples with phenotype data: {len(df_pheno)}")
print()

# Phenotype feature columns
shape_features = ['Median_Shape_Area', 'MAD_Shape_Area', 'Mean_Shape_Area', 'SD_Shape_Area', 'N_Configurations']
colorlab_features = ['Median_ColorLab_L*Mean', 'Median_ColorLab_a*Mean', 'Median_ColorLab_b*Mean', 'Median_ColorLab_ChromaEstimatedMean']
colorhsv_features = ['Median_ColorHSV_HueMean', 'Median_ColorHSV_SaturationMean', 'Median_ColorHSV_BrightnessMean']
other_features = ['Library Plate', 'Incubation Temp (°C)']

all_features = shape_features + colorlab_features + colorhsv_features + other_features

# Remove any columns that don't exist
all_features = [f for f in all_features if f in df_pheno.columns]

# Create correlation matrix
corr_matrix = df_pheno[all_features].corr(method='pearson')

print("="*80)
print("CORRELATIONS WITH MEDIAN_SHAPE_AREA")
print("="*80)
shape_area_corr = corr_matrix['Median_Shape_Area'].sort_values(ascending=False)
for feature, corr_val in shape_area_corr.items():
    if feature != 'Median_Shape_Area':
        print(f"{feature:40} {corr_val:7.4f}")

print()
print("="*80)
print("CORRELATIONS WITH HSV HUE (Median_ColorHSV_HueMean)")
print("="*80)
hue_corr = corr_matrix['Median_ColorHSV_HueMean'].sort_values(ascending=False)
for feature, corr_val in hue_corr.items():
    if feature != 'Median_ColorHSV_HueMean':
        print(f"{feature:40} {corr_val:7.4f}")

print()
print("="*80)
print("CORRELATIONS WITH HSV SATURATION (Median_ColorHSV_SaturationMean)")
print("="*80)
sat_corr = corr_matrix['Median_ColorHSV_SaturationMean'].sort_values(ascending=False)
for feature, corr_val in sat_corr.items():
    if feature != 'Median_ColorHSV_SaturationMean':
        print(f"{feature:40} {corr_val:7.4f}")

print()
print("="*80)
print("CORRELATIONS WITH HSV BRIGHTNESS (Median_ColorHSV_BrightnessMean)")
print("="*80)
bright_corr = corr_matrix['Median_ColorHSV_BrightnessMean'].sort_values(ascending=False)
for feature, corr_val in bright_corr.items():
    if feature != 'Median_ColorHSV_BrightnessMean':
        print(f"{feature:40} {corr_val:7.4f}")

print()
print("="*80)
print("KEY INTERACTIONS (Shape Area with Color Features)")
print("="*80)
print(f"Shape Area vs Hue:       {corr_matrix.loc['Median_Shape_Area', 'Median_ColorHSV_HueMean']:7.4f}")
print(f"Shape Area vs Saturation:{corr_matrix.loc['Median_Shape_Area', 'Median_ColorHSV_SaturationMean']:7.4f}")
print(f"Shape Area vs Brightness:{corr_matrix.loc['Median_Shape_Area', 'Median_ColorHSV_BrightnessMean']:7.4f}")
print()
print(f"Shape Area vs L*:        {corr_matrix.loc['Median_Shape_Area', 'Median_ColorLab_L*Mean']:7.4f}")
print(f"Shape Area vs a*:        {corr_matrix.loc['Median_Shape_Area', 'Median_ColorLab_a*Mean']:7.4f}")
print(f"Shape Area vs b*:        {corr_matrix.loc['Median_Shape_Area', 'Median_ColorLab_b*Mean']:7.4f}")
print(f"Shape Area vs Chroma:    {corr_matrix.loc['Median_Shape_Area', 'Median_ColorLab_ChromaEstimatedMean']:7.4f}")

# Save correlation matrix
corr_matrix.to_csv('/scratch/jstajich/25989115/claude-1181/-bigdata-stajichlab-shared-projects-Rhodotorula-Rhodotorula-Metabolites-feature-extractMS2/206d2d5d-9c2d-470e-be4a-c7cd9006d021/scratchpad/phenotype_correlations.csv')
print()
print("Correlation matrix saved to: phenotype_correlations.csv")

# Look for strong interactions
print()
print("="*80)
print("STRONG POSITIVE CORRELATIONS (r > 0.5)")
print("="*80)
strong_pos = []
for i in range(len(corr_matrix.columns)):
    for j in range(i+1, len(corr_matrix.columns)):
        val = corr_matrix.iloc[i, j]
        if val > 0.5:
            strong_pos.append((corr_matrix.columns[i], corr_matrix.columns[j], val))

strong_pos.sort(key=lambda x: x[2], reverse=True)
for feat1, feat2, val in strong_pos[:15]:
    print(f"{feat1:40} <--> {feat2:40} r={val:.4f}")

print()
print("="*80)
print("STRONG NEGATIVE CORRELATIONS (r < -0.3)")
print("="*80)
strong_neg = []
for i in range(len(corr_matrix.columns)):
    for j in range(i+1, len(corr_matrix.columns)):
        val = corr_matrix.iloc[i, j]
        if val < -0.3:
            strong_neg.append((corr_matrix.columns[i], corr_matrix.columns[j], val))

strong_neg.sort(key=lambda x: x[2])
for feat1, feat2, val in strong_neg[:15]:
    print(f"{feat1:40} <--> {feat2:40} r={val:.4f}")

