# Project Manifest

**Project**: Rhodotorula MS2 Metabolite-Phenotype Analysis  
**Version**: 3.0  
**Last Updated**: 2026-07-02  
**Status**: Complete (Phase 0-3)

## Input Data Inventory

| File | Size | Rows | Type | Status |
|------|------|------|------|--------|
| `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` | 49 MB | 590 samples × 16331 features | Raw MS2 | ✓ Git-tracked |
| `input_data/growth_phenotype_summary_YPD2.csv.gz` | 1.2 MB | 590 samples × 13 phenotypes | Phenotype | ✓ Git-tracked |

## Phase Outputs Inventory

### Phase 0: Data Alignment & QC
| Output File | Size | Rows | Columns | Status |
|-------------|------|------|---------|--------|
| `phase0_metadata_aligned.csv.gz` | 21 KB | 590 | 11 (sample ID, species, plate, phenotypes) | ✓ Complete |
| `phase0_ms2_aligned.csv.gz` | 17 MB | 590 | 7341 (filtered features) | ✓ Complete |
| `phase0_decision.json` | 257 B | - | Config snapshot | ✓ Complete |

**Outputs**: 3 files, 17.0 MB total  
**Reproducibility**: Regenerable from input data  
**Notes**: Species assignment, plate tracking, phenotype alignment

---

### Phase 1: Feature Filtering
| Output File | Size | Rows | Columns | Status |
|-------------|------|------|---------|--------|
| `phase1_phenotype_data.csv.gz` | 21 KB | 590 | 3 (L*, a*, b* color space) | ✓ Complete |
| `phase1_features_filtered.csv.gz` | 29 MB | 590 | 7341 (log-transformed) | ✓ Complete |
| `phase1_feature_metadata.csv.gz` | 149 KB | 7341 | 4 (mean, std, variance, n_nonzero) | ✓ Complete |

**Outputs**: 3 files, 29.2 MB total  
**Reproducibility**: Regenerable via `scripts/phase1_feature_filtering.py`  
**Notes**: Log-transformed abundances, non-zero filtering, variance tracking

---

### Phase 2: Correlation Analysis
| Output File | Size | Rows | Columns | Status |
|-------------|------|------|---------|--------|
| `phase2_all_correlations.csv.gz` | 987 KB | 22,023 | 9 (feature, phenotype, ρ, p, q, tier, etc) | ✓ Complete |
| `phase2_tier1_hits.csv.gz` | 503 KB | 6,847 | 9 (ρ > 0.30, q < 0.05) | ✓ Complete |
| `phase2_tier12_hits.csv.gz` | 661 KB | 10,127 | 9 (ρ > 0.25, q < 0.05) | ✓ Complete |
| `phase2_summary.json` | 349 B | - | Summary stats | ✓ Complete |

**Outputs**: 4 files, 2.2 MB total  
**Reproducibility**: Regenerable via `scripts/phase2_correlation_analysis.py`  
**Notes**: Two-stage Benjamini-Hochberg FDR correction (within-phenotype, then global)

---

### Phase 3: Stratified Species Analysis
| Output File | Size | Rows | Columns | Status |
|-------------|------|------|---------|--------|
| `phase3_discriminant_features.csv.gz` | 5.6 KB | 200 | 3 (feature_id, variance_across_species) | ✓ Complete |
| `phase3_within_species_correlations.csv.gz` | 51 KB | 5,795 | 7 (species, feature, phenotype, ρ, p, q) | ✓ Complete |
| `phase3_phenotype_by_species.csv.gz` | 3.7 KB | 30 | 7 (species × phenotype means/stds) | ✓ Complete |
| `phase3_summary.json` | 126 B | - | Summary stats | ✓ Complete |

**Outputs**: 4 files, 61 KB total  
**Reproducibility**: Regenerable via `scripts/phase3_stratified_species_analysis.py`  
**Notes**: Stratified by 10 major species (n ≥ 3)

---

## Presentation Outputs

| File | Purpose | Status |
|------|---------|--------|
| `presentations/phase3_stratified_analysis.pdf` | Main findings (9 slides) | ✓ Complete |
| `presentations/phase0/` | Data QC plots (8 PNG) | ✓ Complete |
| `presentations/phase2/` | Correlation visualizations (8 PNG) | ✓ Complete |

**Total**: 17 visualization files, 3.8 MB

---

## Documentation Inventory

| Document | Purpose | Status |
|----------|---------|--------|
| `docs/00_START_HERE.md` | Entry point | ✓ Current |
| `docs/ANALYSIS_REPORT.md` | Full methodology | ✓ Current |
| `docs/FINDINGS_SUMMARY.md` | Key results | ✓ Current |
| `docs/MS2_COLOR_CORRELATION_STRATEGY.md` | Analytical approach | ✓ Current |
| `docs/PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md` | Phase 3 findings | ✓ Current |
| `PROJECT_STRUCTURE.md` | This project's organization | ✓ Current |
| `MANIFEST.md` | This file | ✓ Current |
| `DATA_DICTIONARY.md` | Column definitions | ✓ Current |

---

## Quality Metrics

### Data Completeness
- Input samples: 590/590 ✓
- Input features (before filtering): 16,331
- Output features (after filtering): 7,341 (45% retained)
- Species with n≥3: 10 major species
- Phenotypes analyzed: 3 (Lab L*, a*, b*)

### Analysis Completeness
- Phase 0: ✓ Complete
- Phase 1: ✓ Complete
- Phase 2: ✓ Complete (with Benjamini-Hochberg correction)
- Phase 3: ✓ Complete (stratified by species)

### File Integrity
- All CSV files: gzip-compressed (`.csv.gz`)
- All summary stats: JSON format
- All visualizations: PNG (300 DPI) + PDF

---

## Reproduction Instructions

To regenerate all outputs from scratch:

```bash
cd /bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_MS2_pheno_explore

# Setup environment
pip install -r requirements.txt
# OR: conda env create -f environment.yml

# Run pipeline (phases 0-3 in order)
cd scripts
python3 phase0_batch_assessment.py
python3 phase1_feature_filtering.py
python3 phase2_correlation_analysis.py
python3 phase3_stratified_species_analysis.py
python3 phase3_visualizations.py

# All outputs will be in ../results/
```

**Expected runtime**: ~5-10 minutes (depends on system)  
**Disk space required**: ~60 MB (input + output)

---

## Archive & Backup Status

| Tier | Location | Backed Up | Notes |
|------|----------|-----------|-------|
| Source | Git repository (GitHub) | ✓ Yes | Code + small files |
| Results | `results/` folder | ✓ Yes | Regenerable; can delete locally |
| Input | `input_data/` | ✓ Yes | Original files; do not delete |
| Presentations | `presentations/` | ✓ Yes | Publication figures |

---

## Known Limitations & Future Work

1. **Phase 1 visualizations missing**: Only phases 0, 2, 3 have plots
2. **R. mucilaginosa diversity**: Largest species (n=205) shows no within-species correlations due to high genetic/phenotypic diversity
3. **Feature annotation**: Metabolite identities not included (MS/MS matching needed)
4. **Small species correlations**: Phases with n<5 may have inflated correlation strength

---

## Contact & Questions

For questions about this analysis, see `docs/00_START_HERE.md` or `docs/ANALYSIS_REPORT.md`.

