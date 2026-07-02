#!/usr/bin/env python3
import pandas as pd
import sys

# Read the MS2 metadata file (TSV)
ms2_file = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/feature_extractMS2/MS2_samples_combine.extended_metadata.tsv"
strain_file = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/feature_extractMS2/strain_summary_YPD2.csv"

print("Reading MS2 metadata file...")
ms2_df = pd.read_csv(ms2_file, sep='\t')
print(f"MS2 file shape: {ms2_df.shape}")
print(f"MS2 columns: {list(ms2_df.columns[:10])}")

print("\nReading strain summary file...")
strain_df = pd.read_csv(strain_file)
print(f"Strain file shape: {strain_df.shape}")
print(f"Strain columns: {list(strain_df.columns[:10])}")

# Display sample strain values to understand matching
print("\n--- MS2 Standardized Strain values ---")
print(ms2_df['Standardized Strain'].head(10).tolist())

print("\n--- Strain summary Strain column ---")
print(strain_df['Strain'].head(10).tolist())

# Extract strain identifiers for matching, normalizing prefixes
def normalize_strain(strain_str):
    """Normalize strain ID by removing TCFN_ or TFCN_ prefixes"""
    if pd.isna(strain_str):
        return ''
    strain_str = str(strain_str).strip()
    # Remove either TCFN_ or TFCN_ prefix if present
    if strain_str.startswith('TCFN_'):
        strain_str = strain_str[5:]  # Remove 'TCFN_'
    elif strain_str.startswith('TFCN_'):
        strain_str = strain_str[5:]  # Remove 'TFCN_'
    return strain_str

ms2_df['strain_key'] = ms2_df['Standardized Strain'].apply(normalize_strain)
strain_df['strain_key'] = strain_df['Strain'].apply(normalize_strain)

print("\n--- MS2 strain_key values ---")
print(ms2_df['strain_key'].head(15).tolist())

print("\n--- Strain summary strain_key values ---")
print(strain_df['strain_key'].head(15).tolist())

# Get the columns from strain_summary that we want to add (exclude duplicates with MS2)
existing_cols = set(ms2_df.columns)
strain_cols_to_add = [col for col in strain_df.columns if col not in existing_cols and col != 'strain_key']

print(f"\nColumns to add from strain_summary: {strain_cols_to_add}")

# Merge the dataframes
# First, get unique strains from strain_summary to avoid duplicates
strain_summary_unique = strain_df[['strain_key'] + strain_cols_to_add].drop_duplicates(subset=['strain_key'])

print(f"Unique strains in strain_summary: {len(strain_summary_unique)}")

# Keep all rows from MS2, add matching columns from strain_summary
merged_df = ms2_df.merge(
    strain_summary_unique,
    on='strain_key',
    how='left'
)

print(f"\nMerged file shape: {merged_df.shape}")
print(f"Number of rows with strain data added: {merged_df[strain_cols_to_add[0]].notna().sum()}")

# Remove the temporary strain_key column
merged_df = merged_df.drop(columns=['strain_key'])

# Save the output
output_file = "/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/feature_extractMS2/MS2_samples_combine.extended_metadata_with_strain_traits.tsv"
merged_df.to_csv(output_file, sep='\t', index=False)
print(f"\nOutput saved to: {output_file}")
print(f"Final shape: {merged_df.shape}")
print(f"\nNew columns added: {strain_cols_to_add}")
