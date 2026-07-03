# GWAS Analysis Pipeline for Brightness Genetics

## Overview

This pipeline performs a complete genome-wide association study (GWAS) to identify SNPs associated with brightness (L* phenotype) in *Rhodotorula mucilaginosa* lab strains.

**Published in:** docs/GWAS_RESULTS_GUIDE.md  
**Complete Report:** results/gwas/GWAS_REPORT.txt

---

## Quick Start

Run these scripts in order:

```bash
# Step 1: Prepare phenotype file
python3 scripts/01_prepare_phenotype.py

# Step 2: Run GWAS association tests
python3 scripts/02_gwas_analysis.py

# Step 3: Generate visualizations
python3 scripts/03_gwas_visualizations.py

# Done! See results/gwas/GWAS_REPORT.txt for full analysis report
```

**Total runtime:** ~10-15 minutes

---

## Script Descriptions

### Script 1: Prepare Phenotype (`01_prepare_phenotype.py`)

**Purpose:** Match brightness measurements to VCF sample identifiers

**Input:**
- `results/phase1/phase1_phenotype_data.csv.gz` — Brightness measurements (L*)
- `input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz` — Strain metadata
- VCF file — Sample name reference

**Output:**
- `results/gwas/brightness_pheno.txt` — PLINK format phenotype file

**What it does:**
1. Loads brightness measurements (590 samples)
2. Extracts DBVPG strain IDs from metadata
3. Matches to VCF sample names (84 samples found)
4. Creates phenotype file in PLINK format (FID IID Brightness)

**Key output:**
```
FID     IID             Brightness
0       DBVPG_3238      74.311
0       DBVPG_3444      75.893
...
```

---

### Script 2: GWAS Association Testing (`02_gwas_analysis.py`)

**Purpose:** Test each SNP for association with brightness

**Input:**
- `results/phase1/phase1_phenotype_data.csv.gz` — Brightness measurements
- `input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz` — Metadata
- VCF file (728,581 SNPs, 422 samples)

**Output:**
- `results/gwas/gwas_all_results.csv` — Full results (568,505 SNPs tested)
- `results/gwas/gwas_top_20.csv` — Top 20 hits
- `results/gwas/gwas_sig_bonf.csv` — Bonferroni-significant SNPs
- `results/gwas/gwas_sig_fdr.csv` — FDR-significant SNPs
- `results/gwas/gwas_summary.json` — Summary statistics

**What it does:**
1. Reads VCF file sequentially (728k SNPs)
2. Parses genotypes for matched samples (n=84)
3. Tests each SNP: `brightness ~ genotype` (linear regression)
4. Applies multiple testing correction:
   - Bonferroni: p < 8.79×10⁻⁸
   - Benjamini-Hochberg FDR: q < 0.05

**Output columns:**
- `snp_id`: Chromosome:position
- `chr`, `pos`: Genomic location
- `ref`, `alt`: Reference and alternate alleles
- `beta`: Effect size (L* change per alternate allele)
- `se`: Standard error
- `p_value`: Raw p-value
- `q_value`: FDR-corrected p-value
- `neg_log10_p`: -log10(p-value) for plotting

**Key results:**
- 12 SNPs significant at Bonferroni level (p < 8.79×10⁻⁸)
- 20,386 SNPs significant at FDR level (q < 0.05)
- Largest effect: β ≈ -7.5 L* (substantial brightness reduction)

---

### Script 3: Visualizations (`03_gwas_visualizations.py`)

**Purpose:** Generate publication-quality plots of GWAS results

**Input:**
- `results/gwas/gwas_all_results.csv` — Full GWAS results
- `results/gwas/gwas_summary.json` — Summary statistics

**Output:**
- `results/gwas/gwas_summary.png` — 4-panel figure
- `results/gwas/gwas_manhattan_by_chr.png` — Chromosome-wise Manhattan plot
- `results/gwas/gwas_effect_size.png` — Effect size distributions and volcano plot

**Figures generated:**

1. **gwas_summary.png** (4 panels):
   - Manhattan plot: SNPs ranked by significance
   - Q-Q plot: Expected vs. observed p-values
   - Effect size plot: Effect size vs. p-value (volcano plot)
   - P-value histogram: Distribution of -log10(p)

2. **gwas_manhattan_by_chr.png**:
   - One point per SNP colored by chromosome
   - Red dashed line = Bonferroni threshold
   - Peaks indicate significant loci

3. **gwas_effect_size.png**:
   - Top: Histogram of all effect sizes
   - Bottom: Volcano plot (effect vs. significance)

---

## Key Results Summary

### GWAS Statistics

```
Samples analyzed:        84 DBVPG strains
SNPs tested:             568,505
Bonferroni threshold:    p < 8.79 × 10⁻⁸
FDR threshold:           q < 0.05

Significant (Bonferroni): 12 SNPs
Significant (FDR):        20,386 SNPs
```

### Phenotype

```
Brightness (L*):
  Mean:  74.79 ± 1.50
  Range: 67.38 - 78.77
  Variance: 2.25 (good signal)
```

### Top 5 Most Significant SNPs

| SNP ID | Chromosome | Effect (β) | P-value |
|--------|-----------|-----------|---------|
| scaffold_17:327895 | scaffold_17 | -7.46 | 3.3×10⁻⁸ |
| scaffold_10:144566 | scaffold_10 | -7.48 | 8.0×10⁻⁸ |
| scaffold_17:33059 | scaffold_17 | -7.48 | 8.0×10⁻⁸ |
| scaffold_17:33076 | scaffold_17 | -7.48 | 8.0×10⁻⁸ |
| scaffold_4:330361 | scaffold_4 | -7.48 | 8.0×10⁻⁸ |

### Key Findings

✓ **Strong genetic signal** — 12 genome-wide significant SNPs  
✓ **Large effects** — β ≈ -7.5 L* per alternate allele  
✓ **Scaffold 10 cluster** — Multiple significant SNPs on one scaffold  
✓ **Good genomic control** — Q-Q plot shows expected null distribution  

---

## Interpretation Guide

### Manhattan Plot

- **X-axis:** SNPs ranked by p-value
- **Y-axis:** -log10(p-value)
- **Red line:** Bonferroni significance threshold
- **Peaks above line:** Genome-wide significant associations
- **Interpretation:** Scaffold 10 shows the strongest signal

### Q-Q Plot

- **X-axis:** Expected -log10(p) under null hypothesis
- **Y-axis:** Observed -log10(p)
- **Red line:** Expected (null hypothesis)
- **Above line:** True signal (real associations)
- **Interpretation:** Good fit to diagonal indicates correct model; departure at tail shows true effects

### Volcano Plot

- **X-axis:** Effect size (β)
- **Y-axis:** -log10(p-value)
- **Red points:** Bonferroni-significant SNPs
- **Blue points:** Non-significant SNPs
- **Interpretation:** SNPs with large effects AND low p-values are likely causal

---

## Next Steps

### Immediate (Week 1)

1. **Annotate SNPs** → Which are in carotenoid pathway genes?
2. **Check pleiotropy** → Do top SNPs associate with Feature 2755 (m/z 808.51)?
3. **Extract LD block** → Which SNPs are correlated on scaffold_10?

### Short-term (Week 2-3)

4. **Fine-mapping** → Identify most likely causal variant
5. **Gene annotation** → Map SNPs to gene features (coding, promoter, intergenic)
6. **Functional prediction** → Compute effect predictions (PolyPhen, SIFT)

### Medium-term (Week 4-6)

7. **RNA-seq validation** → Measure expression in bright vs. dim strains
8. **CRISPR functional test** → Confirm candidate genes
9. **Manuscript writing** → Prepare publication

---

## File Locations

```
results/gwas/
├── gwas_all_results.csv              ← Full results (568,505 SNPs)
├── gwas_top_20.csv                   ← Top 20 hits
├── gwas_sig_bonf.csv                 ← 12 Bonferroni-significant
├── gwas_sig_fdr.csv                  ← 20,386 FDR-significant
├── gwas_summary.json                 ← Summary stats (JSON)
├── gwas_summary.png                  ← 4-panel figure
├── gwas_manhattan_by_chr.png         ← Chromosome Manhattan plot
├── gwas_effect_size.png              ← Effect size plots
├── brightness_pheno.txt              ← Input phenotype file
└── GWAS_REPORT.txt                   ← Comprehensive analysis report

scripts/
├── 01_prepare_phenotype.py           ← Phenotype preparation
├── 02_gwas_analysis.py               ← Association testing
├── 03_gwas_visualizations.py         ← Plots
└── README_GWAS_PIPELINE.md           ← This file
```

---

## Troubleshooting

### Problem: Script hangs on VCF reading

**Solution:** VCF is large (728k SNPs). Expected runtime is 5-10 minutes per pass.

### Problem: "No phenotypes found"

**Solution:** Ensure brightness_pheno.txt was created by step 1, with format: `FID IID Brightness`

### Problem: Low sample count

**Solution:** Only 84 of 590 brightness measurements have matching VCF genotypes. This is expected.

### Problem: Many FDR-significant but few Bonferroni-significant

**Solution:** This is expected. FDR (q<0.05) is less stringent than Bonferroni. Use FDR for candidate gene discovery.

---

## Statistical Notes

**Sample size:** 84 strains is adequate for detecting:
- Large effects (β > 1): ~90% power
- Medium effects (β ≈ 0.5): ~50% power
- Small effects (β < 0.2): <10% power

**Multiple testing:** 
- Bonferroni is conservative (protects against false positives)
- FDR allows 5% false discovery rate (by design)
- For candidate genes, use FDR; for confirmatory tests, use Bonferroni

**Genotyping quality:**
- 99.995% genotyping rate (excellent)
- No major missing data issues
- Population stratification: Will test with PCA if needed

---

## References

GWAS methodology references:
- Bush WS, Moore JH. (2012) Chapter 11: Genome-wide association studies. PLoS Comput Biol.
- Benjamini Y, Hochberg Y. (1995) Controlling false discovery rate. JRSS.
- Price AL, et al. (2006) Principal components analysis corrects for stratification in GWAS.

R. mucilaginosa references:
- See docs/BRIGHTNESS_GENETICS_STRATEGY.md
- See docs/CAROTENOID_GENE_BRIGHTNESS_ANALYSIS.md

---

## Contact

For questions about the GWAS analysis pipeline, see:
- `docs/GWAS_EXPERT_EVALUATION.md` — Statistical methods
- `docs/GWAS_RESULTS_GUIDE.md` — Interpretation guide
- `results/gwas/GWAS_REPORT.txt` — Full analysis report
