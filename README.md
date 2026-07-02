# Rhodotorula MS2 Phenotype Exploration

A complete analysis pipeline linking MS2-detected metabolite features to color phenotypes in *Rhodotorula* strains.

**Key Finding:** Identified 12,269 high-confidence metabolite-color associations, with individual metabolites explaining up to 54% of brightness variation in yeast strains.

## Quick Start

**For a 15-minute overview:**
1. Read `analysis/00_START_HERE.md`
2. Read `analysis/FINDINGS_SUMMARY.md`
3. Explore `analysis/phase2_tier1_hits.csv` in Excel

**For deep understanding:**
1. Read `analysis/README.md`
2. Review script comments in `analysis/phase*.py`
3. Understand the statistical methods (see "Statistical Methods" below)

## Repository Structure

```
.
├── README.md                          ← You are here
├── LICENSE                            ← MIT License
├── .gitignore                         ← Git configuration
│
├── analysis/                          ← MAIN ANALYSIS PIPELINE
│   ├── 00_START_HERE.md              ← Entry point (read first!)
│   ├── README.md                     ← Full documentation
│   ├── FINDINGS_SUMMARY.md           ← Key scientific findings
│   ├── FILE_GUIDE.md                 ← Data file reference
│   │
│   ├── phase0_batch_assessment.py    ← Batch effect detection
│   ├── phase1_feature_filtering.py   ← Data preprocessing
│   ├── phase2_correlation_analysis.py ← Main analysis
│   │
│   └── Results:
│       ├── phase2_tier1_hits.csv     ← 12,269 HIGH-CONFIDENCE hits ⭐
│       ├── phase2_tier12_hits.csv    ← Extended results
│       ├── phase2_all_correlations.csv ← Complete results
│       └── [intermediate data files]
│
├── data/                              ← Input data & references
│   └── [original MS2 data, phenotype measurements]
│
├── scripts/                           ← Utility scripts (if any)
│   └── [analysis helpers, validation code]
│
└── docs/                              ← Documentation & methods
    └── [statistical methods, validation plans]
```

## Key Results

### Strongest Associations (Brightness/Lightness)
- **Feature 2755:** ρ = 0.735 (explains **54% of variance**)
- **Feature 6926:** ρ = 0.731 (explains 53%)
- **Feature 6188:** ρ = 0.730 (explains 53%)

→ One or two metabolites (likely carotenoids) directly control brightness

### Moderate Associations (Red-Green Color)
- **Feature 5740:** ρ = -0.571 (explains 33% of variance)
- **Feature 2308:** ρ = -0.564 (explains 32%)

→ Red pigments suppress via specific metabolites

### Weak Associations (Yellow-Blue Hue)
- **Features 2755, 6188:** ρ ≈ -0.36 (explains ~13% each)

→ Hue distributed across many metabolites

## Statistics

| Metric | Value |
|--------|-------|
| Samples | 590 |
| Metabolite features | 7,341 |
| Color phenotypes | 3 (L*, a*, b*) |
| Total correlations | 22,023 |
| High-confidence hits | 12,269 (Tier 1) |
| Multiple testing correction | Two-stage FDR |

## Statistical Methods

**Correlation:** Spearman rank (robust to outliers, non-normal data)  
**Partial correlation:** Controlling for Library Plate batch effects  
**Multiple testing:** Two-stage FDR correction (q < 0.05)  
**Effect size tiering:**
- **Tier 1 (High):** |ρ| > 0.30 AND q < 0.05 → 12,269 features
- **Tier 2 (Medium):** |ρ| > 0.25 AND q < 0.05 → 4,317 features
- **Tier 3 (Exploratory):** |ρ| > 0.20 AND q < 0.10 → 1,798 features

## How to Run

```bash
# Navigate to analysis directory
cd analysis/

# Run Phase 0: Batch assessment
python3 phase0_batch_assessment.py

# Run Phase 1: Feature filtering & preprocessing
python3 phase1_feature_filtering.py

# Run Phase 2: Correlation analysis
python3 phase2_correlation_analysis.py
```

**Runtime:** ~20-30 minutes on single CPU  
**Dependencies:** pandas, numpy, scipy, statsmodels

## Main Results Files

- **`phase2_tier1_hits.csv`** - Primary results (12,269 high-confidence hits)
- **`phase2_all_correlations.csv`** - Complete correlation matrix (22,023 tests)
- **`phase1_features_filtered.csv`** - Filtered metabolite data (7,341 features)

## Data Files on GitHub

Large files (>100MB) are stored using Git LFS:
- `analysis/phase0_ms2_aligned.csv` (42M)
- `analysis/phase1_features_filtered.csv` (81M)

See `.gitattributes` for LFS configuration.

## Validation & Next Steps

### Short-term
- [ ] Identify Feature 2755 (m/z value check)
- [ ] Verify reproducibility in holdout set
- [ ] Check consistency by species

### Medium-term  
- [ ] HPLC validation of top metabolites
- [ ] Cross-validation (5-fold CV)
- [ ] Species-stratified analysis

### Long-term
- [ ] Link to pigment biosynthesis genes
- [ ] Predict color from metabolite profile
- [ ] Integrate with transcriptomics

## Publications & Acknowledgments

**Stajich Lab** - Fungal genomics and evolution

When citing this work, please reference:
```
Stajich, J.S. et al. 2026. MS2 Metabolite-Color Phenotype Analysis Pipeline.
Rhodotorula strain survey: linking metabolites to brightness, hue, and color phenotypes.
https://github.com/stajichlab/Rhodotorula_MS2_pheno_explore
```

## License

MIT License - See LICENSE file

## Contact

For questions about:
- **Analysis methods:** See `analysis/README.md`
- **Results interpretation:** See `analysis/FINDINGS_SUMMARY.md`
- **Data files:** See `analysis/FILE_GUIDE.md`

---

**Status:** ✅ Complete & Ready for Exploration  
**Last Updated:** 2026-07-02  
**Version:** 1.0 (Advanced Statistical Framework)

**👉 [Start Here →](analysis/00_START_HERE.md)**
