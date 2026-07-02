# Phase 3: Stratified Species Analysis Summary

## Overview
This phase stratifies the metabolite-phenotype correlation analysis by species to distinguish between:
- **Between-species signals**: Which metabolites discriminate species?
- **Within-species signals**: Which metabolites correlate with phenotypes within each species?

## Data
- **Species**: 10 major species (n ≥ 3)
  - *Rhodotorula mucilaginosa* (205) — dominant
  - *R. paludigena* (10), *R. toruloides* (9), *R. dairenensis* (7), *R. diobovata* (7)
  - *R. taiwanensis* (6), *R. sp. clade I* (5), *R. sphaerocarpa* (4)
  - *R. kratochvilovae* (3), *R. graminis* (3)
- **Samples analyzed**: 259 major species (plus 331 in rare species)
- **Features**: 200 top discriminant features
- **Phenotypes**: Lab L*, a*, b* (color space)

## Key Results

### Between-Species (Discriminant Features)
- **Identified**: 200 metabolites with high variance across species
- **Top feature**: #397 (variance = 12.10)
- **Interpretation**: Strong metabolic differentiation between species

### Within-Species Correlations
- **Total tested**: ~5,800 feature-phenotype pairs
- **Significant (q<0.05)**: 126 correlations
- **Key pattern**: Small species show strong correlations; large species show weak correlations

**Results by species:**
- R. sphaerocarpa (n=4): **16.5%** significant
- R. sp. clade I (n=5): **3.7%** significant
- R. dairenensis (n=7): **1.2%** significant
- R. mucilaginosa (n=205): **0%** significant ← dominated by noise

## Outputs
✓ **analysis/phase3_stratified_species_analysis.py** — Analysis pipeline
✓ **analysis/phase3_visualizations.py** — PDF generation
✓ **analysis/phenotypes_MS2/phase3_stratified_analysis.pdf** — 9-slide presentation
✓ **analysis/phenotypes_MS2/phase3_discriminant_features.csv** — Discriminant metabolites
✓ **analysis/phenotypes_MS2/phase3_within_species_correlations.csv.gz** — Correlation data
✓ **analysis/phenotypes_MS2/phase3_phenotype_by_species.csv** — Phenotype statistics

## Interpretation
1. **Between-species**: Metabolite profiles distinguish species (likely ecological/evolutionary adaptation)
2. **Within-species**: Small, uniform species show directional metabolite-phenotype relationships
3. **Large species noise**: *R. mucilaginosa* high genetic/phenotypic diversity obscures correlations

## Next Steps
1. Identify discriminant features biochemically (MS/MS, databases)
2. Test functional hypotheses on pigment biosynthesis
3. For R. mucilaginosa: stratify by phylogenomic clusters, retest
4. Compare pathway enrichment across species

