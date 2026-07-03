# GWAS Results Guide: What You'll Get

**Status:** Pipeline running (estimated 30-60 min)  
**Output Location:** `results/gwas/`

---

## What the Pipeline Produces

### 1. Data Files

#### `gwas_results_all_snps.csv`
**All tested SNPs with statistics**

Columns:
- `snp_id`: Chromosome:position identifier
- `ref`, `alt`: Reference and alternate alleles
- `beta`: Effect size (change in L* per alt allele)
- `se`: Standard error of effect size
- `t_stat`: t-statistic from regression
- `p_value`: Raw p-value
- `q_value_bh`: Benjamini-Hochberg FDR-corrected p-value
- `n_samples`: Number of samples used in test
- `neg_log10_p`: -log10(p-value) for Manhattan plot
- `significant_bonf`: True if Bonferroni-significant
- `significant_fdr`: True if FDR-significant

**How to use:**
```python
import pandas as pd
results = pd.read_csv('results/gwas/gwas_results_all_snps.csv')

# Top 10 hits
print(results.nsmallest(10, 'p_value'))

# Export significant SNPs
sig_snps = results[results['significant_bonf']]
sig_snps.to_csv('top_snps_bonferroni.csv')
```

---

#### `gwas_results_top_20.csv`
**Top 20 most significant SNPs**

Use this for presentations and quick reference.

---

#### `gwas_summary.json`
**Summary statistics**

```json
{
  "n_matched_samples": 85,
  "n_snps_tested": 1234,
  "n_significant_bonf": 5,
  "n_significant_fdr": 23,
  "bonferroni_threshold": 4.05e-05,
  "fdr_threshold": 0.05,
  "brightness_mean": 74.32,
  "brightness_std": 1.45,
  "brightness_range": [71.42, 77.45]
}
```

---

### 2. Visualization: `gwas_summary.png`

**4-panel plot:**

1. **Manhattan Plot** (top-left)
   - X-axis: SNPs (ranked by p-value)
   - Y-axis: -log10(p-value)
   - Red dashed line: Bonferroni threshold (p < 4×10⁻⁵)
   - Orange dashed line: p < 0.05 (uncorrected)
   - Points above red line = genome-wide significant

2. **Q-Q Plot** (top-right)
   - Compares observed vs. expected p-value distribution
   - Red dashed line: Expected (null hypothesis)
   - Points above line = more significant than expected (signal!)
   - Points below line = less significant (good genomic control)

3. **Effect Size vs. Significance** (bottom-left)
   - X-axis: β (effect size, L* change per alt allele)
   - Y-axis: -log10(p-value)
   - Red points: Bonferroni-significant SNPs
   - Blue points: Non-significant SNPs
   - Useful for identifying SNPs with large effects

4. **P-value Histogram** (bottom-right)
   - Shows distribution of p-values across all tested SNPs
   - Skewed toward low p-values = true signal present
   - Uniform distribution = no signal (or poor model)

---

## How to Interpret Results

### Scenario A: Strong Signal (5-20 significant SNPs, p < 10⁻⁴)

**Interpretation:**
- Carotenoid pathway SNPs directly control brightness
- These are likely causal or in tight LD with causal variants

**Next steps:**
1. Validate in independent cohort (if available)
2. Use fine-mapping (Bayesian) to identify exact causal SNP
3. Check if SNPs are in promoters (regulatory) or coding (structural)
4. Test pleiotropy: do these SNPs associate with Feature 2755?

**Publication story:**
"We identified X carotenoid pathway SNPs associated with brightness variation (p < 5×10⁻⁵), explaining Y% of phenotypic variance..."

---

### Scenario B: Moderate Signal (10-50 SNPs, q < 0.05 FDR)

**Interpretation:**
- Carotenoid pathway contributes, but effect sizes are moderate
- Likely multiple loci with similar effect sizes

**Next steps:**
1. Pool all significant SNPs → pathway score
2. Test: Does aggregate pathway score predict brightness better?
3. Look for epistasis (SNP interactions)
4. RNA-seq experiment to measure gene expression

**Publication story:**
"Multiple carotenoid pathway variants together predict brightness variation (FDR q < 0.05)..."

---

### Scenario C: Weak/No Signal (no SNPs, p > 0.05)

**Interpretation:**
- SNPs in carotenoid pathway do NOT explain brightness
- Genetic basis may be:
  - Regulatory elements (promoters, enhancers) not captured
  - Non-coding RNA
  - Gene expression/epigenetics
  - Other metabolic pathways

**Next steps:**
1. Consider RNA-seq approach instead
2. Look for structural variants (CNVs, inversions)
3. Expand to broader metabolic pathways
4. Investigate trans-regulatory loci (genome-wide GWAS)

**Publication story:**
"Carotenoid pathway SNPs alone do not explain brightness; we hypothesize regulation via..."

---

## Key Metrics to Report

| Metric | Value | Interpretation |
|--------|-------|-----------------|
| **Sample size** | 85 | Moderate power for fungi |
| **SNPs tested** | ~1000-5000 | Candidate gene approach |
| **Bonferroni threshold** | ~5×10⁻⁵ | Multiple testing correction |
| **FDR threshold** | 0.05 | False discovery rate control |
| **Brightness range** | 71-77 L* | Wide phenotypic variance (good signal) |
| **QQ plot fit** | ? | Genomic control (ideally 1.0 slope) |

---

## Quality Checks

### 1. Check PCA worked
- Look at log output: "PCA variance explained (PC1-3): XX%"
- Should be 20-40% for this dataset
- If < 5%: possible issue with genotype matrix

### 2. Check sample matching
- Log should show: "Matched samples: 85 lab strains in VCF"
- If < 50: possible strain ID mismatch

### 3. Check Q-Q plot
- Points should generally follow red diagonal line
- Minor departure at tail (high significance) is OK
- Major departure everywhere = model issues

### 4. Check effect sizes
- Beta values should range from -1.0 to +1.0
- (Genotypes are 0/1/2, so per-allele effect is typically small)
- If some betas are >10: possible numerical issues

---

## Follow-Up Analyses

### If you find significant SNPs:

**1. Functional annotation**
```bash
# Which SNPs are in promoters vs. coding?
# Use gene coordinates from GFF3 to classify
```

**2. Pleiotropy check**
```python
# For each significant SNP, test association with Feature 2755
# If same SNP affects brightness AND metabolite:
#   → suggests direct causal effect
```

**3. LD visualization**
```bash
# Extract LD block containing significant SNP
# Visualize: which SNPs are correlated?
# Narrows down to likely causal variant
```

**4. Fine-mapping (Bayesian)**
```bash
# Use SuSiE or CAVIAR to get posterior probability
# "SNP_A has 70% chance of being causal"
```

---

## Troubleshooting

### Problem: VCF not found
**Solution:** Check path in script:
```python
vcf_path = '/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/vcf/RmucY2510_v2.All.SNP.combined_selected.vcf.gz'
```

### Problem: Only 10 samples matched
**Solution:** Check VCF sample naming
```python
# Add debug: print vcf_samples[:20]
# Check if DBVPG format matches (should be DBVPG_3238, etc.)
```

### Problem: No carotenoid SNPs found
**Solution:** Script falls back to all SNPs with MAF > 5%
- This is conservative but works
- You'll get genome-wide signal (not pathway-specific)
- Can filter results post-hoc to carotenoid genes

### Problem: PCA fails
**Solution:** Script catches error and uses sklearn PCA instead
- Check log for error message
- Usually just slower, not incorrect

---

## Statistical Notes

### Power calculation
With 85 samples and 2-point phenotypic variance, we can detect:
- **Large effect** (β = 2 L* per alt allele): ~90% power, p < 0.05
- **Medium effect** (β = 1.0 L*): ~70% power
- **Small effect** (β = 0.5 L*): ~30% power

Sample size is adequate for a candidate gene study.

### Multiple testing burden
- **If 1000 SNPs tested:** Bonferroni threshold = 0.05 / 1000 = 5×10⁻⁵
- **If 5000 SNPs tested:** Bonferroni threshold = 0.05 / 5000 = 1×10⁻⁵
- **FDR (less stringent):** Adjusts based on observed p-value distribution
  - Allows some false positives (5% expected false discovery rate)
  - More power than Bonferroni for large effect sizes

---

## Next Paper/Presentation

**Title:** "SNP Variants in Carotenoid Biosynthesis Genes Associated with Pigmentation in *Rhodotorula mucilaginosa*"

**Abstract template:**
"We performed genome-wide association analysis in 85 *R. mucilaginosa* DBVPG reference strains to identify genetic determinants of brightness (L* phenotype). We tested [X] SNPs in carotenoid pathway genes using linear regression with population structure correction. We identified [Y] SNPs significantly associated with brightness (p < 5×10⁻⁵), including [Z] non-synonymous variants in [gene names]. These SNPs explain [W]% of phenotypic variance. Our results suggest [biological interpretation]."

---

## Files Generated

```
results/gwas/
├── gwas_results_all_snps.csv       ← Full results table
├── gwas_results_top_20.csv         ← Top hits
├── gwas_summary.json               ← Summary stats
├── gwas_summary.png                ← 4-panel visualization
└── (gwas_brightness.log if ran)    ← Execution log
```

---

**When ready, let me know the results and we can:**
1. Visualize top SNPs on genomic map
2. Check pleiotropy with Feature 2755
3. Plan validation strategy
4. Write up findings
