# Project Structure Guide

## Overview
This project analyzes Rhodotorula MS2 metabolite data to identify associations between metabolite abundance and color phenotypes, stratified by species.

## Directory Organization

```
Rhodotorula_MS2_pheno_explore/
├── scripts/                          # All Python analysis pipelines
│   ├── phase0_batch_assessment.py   # Data preprocessing & quality control
│   ├── phase1_feature_filtering.py  # Feature selection & cleaning
│   ├── phase2_correlation_analysis.py  # Metabolite-phenotype correlations
│   ├── phase3_stratified_species_analysis.py  # Species-stratified analysis
│   ├── phase3_visualizations.py     # Slide deck generation
│   └── [utility scripts]            # Metadata merging, analysis tools
│
├── results/                          # All data outputs organized by phase
│   ├── phase0/                      # Aligned data, metadata
│   ├── phase1/                      # Filtered features, phenotype data
│   ├── phase2/                      # Correlations, tiered hits, summary stats
│   └── phase3/                      # Discriminant features, species-specific correlations
│
├── presentations/                    # PDF slide decks & visualization outputs
│   ├── phase0/                      # Data exploration plots
│   ├── phase2/                      # Phenotype-feature correlation plots
│   └── phase3_stratified_analysis.pdf  # Main findings presentation
│
├── docs/                             # Documentation & analysis summaries
│   ├── 00_START_HERE.md             # Entry point documentation
│   ├── ANALYSIS_REPORT.md           # Complete analysis narrative
│   ├── FINDINGS_SUMMARY.md          # Key results & interpretation
│   ├── MS2_COLOR_CORRELATION_STRATEGY.md  # Analytical approach
│   └── PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md  # Phase 3 findings
│
├── input_data/                       # Raw input files (gzip compressed)
│   ├── Rhodotorula_MS2_aligned_features_ms2.csv.gz
│   └── growth_phenotype_summary_YPD2.csv.gz
│
└── README.md, LICENSE, SETUP_INSTRUCTIONS.md  # Project metadata
```

## Data Files

### Phase 0: Data Alignment & QC
- `phase0_metadata_aligned.csv.gz` — Sample metadata with species, plate IDs
- `phase0_ms2_aligned.csv.gz` — Raw MS2 feature intensity data

### Phase 1: Feature Filtering
- `phase1_phenotype_data.csv.gz` — Clean phenotype data (Lab L*, a*, b*)
- `phase1_features_filtered.csv.gz` — Filtered metabolite features
- `phase1_feature_metadata.csv.gz` — Feature annotations & statistics

### Phase 2: Correlation Analysis
- `phase2_all_correlations.csv.gz` — All Spearman correlations (tiered by effect size)
- `phase2_tier1_hits.csv.gz` — High-confidence hits (ρ>0.30, q<0.05)
- `phase2_tier12_hits.csv.gz` — High + medium confidence (ρ>0.25, q<0.05)
- `phase2_summary.json` — Summary statistics

### Phase 3: Stratified Species Analysis
- `phase3_discriminant_features.csv.gz` — 200 features varying across species
- `phase3_within_species_correlations.csv.gz` — Species-specific correlations
- `phase3_phenotype_by_species.csv.gz` — Phenotype means/stds by species
- `phase3_summary.json` — Summary statistics

## Running the Analysis

1. **Start here**: Read `docs/00_START_HERE.md`
2. **Setup**: Run `bash SETUP_INSTRUCTIONS.md`
3. **Execute phases** (in order):
   ```bash
   cd scripts
   python3 phase0_batch_assessment.py
   python3 phase1_feature_filtering.py
   python3 phase2_correlation_analysis.py
   python3 phase3_stratified_species_analysis.py
   python3 phase3_visualizations.py  # Generates PDF presentation
   ```

## Key Findings

See `docs/PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md` for latest results:

- **Between-species**: 200 metabolites discriminate species (high variance)
- **Within-species**: 126 significant metabolite-phenotype correlations
- **Sample size effect**: Small species show strong correlations (16.5%), large species show weak (0%)
- **PDF presentation**: `presentations/phase3_stratified_analysis.pdf`

## File Format Notes

- **Data files**: All CSV files are gzip-compressed (`.csv.gz`) for storage efficiency
- **Scripts**: Python 3.9+, requires: pandas, scipy, scikit-learn, matplotlib, seaborn, statsmodels
- **Results**: JSON for numerical summaries, gzip-compressed CSV for large tables

## Questions?

See `docs/ANALYSIS_REPORT.md` for detailed methodology and interpretation.

