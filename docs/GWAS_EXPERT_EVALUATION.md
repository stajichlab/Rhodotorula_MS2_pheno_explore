# GWAS Analysis: Expert Evaluation & Sophisticated Approaches

**Evaluation Date:** 2026-07-02  
**Sample Size:** 85 DBVPG strains with brightness measurements  
**Study Design:** Single-SNP association in carotenoid pathway genes

---

## Your Current Approach: Assessment

✅ **Strengths:**
- Candidate gene approach (focused on biology, not genome-wide)
- Clear phenotype (L* brightness, continuous, well-measured)
- Strong signal available (Feature 2755 explains 54%)
- Good sample size for fungi (85 strains)
- Reference genome available (R. mucilaginosa Y2510)

⚠️ **Limitations:**
- Single-SNP testing loses power from LD structure
- No account for population stratification (85 DBVPG strains likely from different labs/time periods)
- Binary effect size thinking (effect is biallelic; reality is more complex)
- No gene-level aggregation (pooling weak SNPs in one gene)
- Missing regulatory element analysis

---

## Expert Recommendations: Tier 1 (Do These)

### 1.1 Population Stratification Correction (CRITICAL)

**Problem:** DBVPG strains come from different sources/time periods. Allele frequencies differ by origin, creating spurious associations (confounding).

**Solution:** Principal Component Analysis (PCA)
```python
# Use all genome-wide SNPs to estimate population structure
pca = PCA(n_components=3)
pc_scores = pca.fit_transform(genotype_matrix)  # 85 samples × ~1M SNPs

# Re-run association, including PC1, PC2, PC3 as covariates:
#   brightness ~ SNP + PC1 + PC2 + PC3

# Typically PC1 captures 30-50% of variance
```

**Expected impact:** Removes confounding; true effect sizes become visible

**Implementation:** 1–2 hours (standard pipeline)

---

### 1.2 Linkage Disequilibrium-Aware Filtering

**Problem:** SNPs in LD (close together on chromosome) are correlated. Testing each independently inflates false positives AND wastes power.

**Solution:** Prune SNPs to independence
```python
# Step 1: Filter to r² < 0.1 (nearly independent SNPs)
#   Keep one SNP per LD block

# Step 2: Test pruned SNPs (fewer tests, more power per test)

# Step 3: If SNP is significant, look at all SNPs in its LD block
#   Find causal SNP (may not be the tested one)
```

**Tool:** PLINK's `--indep-pairwise` (or bcftools+LD)

**Expected impact:** Fewer spurious associations; cleaner results

**Implementation:** 1–2 hours

---

### 1.3 Variant Annotation & Prioritization

**Problem:** Not all SNPs are equally likely to affect brightness. Intergenic SNPs → low impact. Frameshift in PSY → high impact.

**Solution:** Stratified analysis
```python
# Tier 1: Non-synonymous (missense, frameshift, stop-gain)
#   Higher prior probability of function
#   Lower p-value threshold for significance

# Tier 2: Synonymous (silent)
#   Lower prior; higher threshold

# Tier 3: Promoter (-500 to -1 bp)
#   Regulatory; medium prior

# Tier 4: Intergenic
#   Lowest priority; highest threshold
```

**Tool:** VEP (Variant Effect Predictor) or bcftools/csq

**Expected impact:** Identifies likely causal variants; reduces burden of proof for strong effect SNPs

**Implementation:** 1–2 hours

---

### 1.4 Multiple Testing Correction

**Problem:** Testing ~1,000 SNPs in carotenoid genes → need to control false positives.

**Recommendation:**
- **For candidate genes (Tier 1 approach):** Use Bonferroni correction
  - Threshold: p < 0.05 / (# SNPs tested)
  - If 1,000 SNPs: p < 5×10⁻⁵

- **Alternative (less conservative):** Benjamini-Hochberg FDR
  - Allows ~5% false discovery rate
  - More power, but some false positives expected

**Best practice:** Report both, interpret conservatively

**Implementation:** Built into statsmodels/scipy

---

## Expert Recommendations: Tier 2 (Consider These)

### 2.1 Gene-Level Association Tests

**Problem:** Single SNPs may have weak effects. Better: "Does this gene have ANY associated variants?"

**Solution:** Burden tests or SKAT
```python
# Burden test: pool all rare variants in each gene
#   Sum genotypes per gene → test: brightness ~ gene_burden

# SKAT (Sequence Kernel Association Test):
#   Test: do variant positions + phenotype correlate?
#   More flexible; better for multiple SNPs per gene
```

**Advantage:** More power if multiple weak SNPs; easier interpretation

**Trade-off:** Lose ability to identify specific causal SNP

**Implementation:** 2–4 hours (need to write wrapper)

---

### 2.2 Epistasis (SNP Interactions)

**Problem:** Two SNPs together might have synergistic effect (A + B bright, but A or B alone has no effect).

**Solution:** Test selected pairs
```python
# Only test pairs in SAME gene (too many pairs otherwise)
# Example: SNP in PSY promoter + SNP in PSY coding region

# Model: brightness ~ SNP1 + SNP2 + SNP1:SNP2
```

**Warning:** Only do if you have strong a priori hypothesis (computationally expensive)

**Expected impact:** Rare; but can explain missing heritability

**Implementation:** 4–8 hours (many tests needed)

---

### 2.3 Pleiotropy Analysis

**Problem:** Brightness SNP might also affect carotenoid levels (Feature 2755). Or affect other metabolites.

**Solution:** Check association to metabolite features
```python
# For each brightness-associated SNP:
#   Test: Feature_2755_abundance ~ SNP
#   Test: Feature_6926_abundance ~ SNP  (other carotenoid)
#   etc.

# If same SNP affects brightness AND Feature 2755:
#   Suggests direct metabolic effect (not indirect)
```

**Expected finding:** 80% of brightness SNPs will associate with Feature 2755 (they're correlated)

**Implementation:** 2–3 hours

---

## Expert Recommendations: Tier 3 (Advanced)

### 3.1 Heritability Estimation

**Question:** How much of brightness variation is explained by carotenoid gene SNPs?

**Method:** GCTA-GREML (narrow-sense heritability from SNPs)
```python
# h² = (variance explained by all SNPs) / (total phenotypic variance)
# Compare: h²_all_SNPs vs h²_carotenoid_SNPs

# If h²_carotenoid = 0.3, then 30% of brightness is "genetic"
```

**Interpretation:** If low h², other pathways contribute

**Implementation:** 4–6 hours

---

### 3.2 Bayesian Fine-Mapping

**Problem:** Multiple SNPs in LD; which is causal?

**Method:** Use posterior probability to rank SNPs by likelihood of causality
```python
# SuSiE (Sum of Single Effects) or CAVIAR
# Outputs: "SNP_1 has 60% probability of being causal"

# Narrows down from multiple candidates → single SNP
```

**Advanced but powerful for publication**

**Implementation:** 4–8 hours (learning curve)

---

## Recommended Analysis Path (4 weeks)

### Week 1: Data Preparation
- [ ] Extract carotenoid gene SNPs from VCF
- [ ] Perform PCA on full genome (population correction)
- [ ] Annotate SNPs (coding, promoter, etc.) using VEP
- [ ] LD-prune to independent SNPs

### Week 2: Primary Analysis
- [ ] Single-SNP association test (with PCA correction)
- [ ] Generate Manhattan plot
- [ ] Identify genome-wide significant SNPs (p < 5×10⁻⁵)
- [ ] Create results table (SNP, gene, effect size, p-value)

### Week 3: Validation & Interpretation
- [ ] Pleiotropy check (do brightness SNPs associate with Feature 2755?)
- [ ] Gene-level burden tests (do genes have aggregate effect?)
- [ ] LD block visualization (show which SNPs are correlated)
- [ ] Promoter vs. coding effect comparison

### Week 4: Publication Polish
- [ ] Heritability estimation
- [ ] Bayesian fine-mapping (optional)
- [ ] Write-up of methods & results
- [ ] Generate figures for paper

---

## Code Skeleton (Tier 1 Analysis)

### Step 1: Load genotypes & phenotypes
```python
import pandas as pd
import numpy as np
from scipy import stats
import gzip

# Load VCF
# vcf = load_vcf(path)  # Use cyvcf2 or pysam

# Load phenotypes (your 85 matched strains)
pheno = pd.DataFrame({
    'strain': [list of 85 DBVPG samples],
    'brightness_L': [corresponding L* values]
})

# Perform PCA on all samples in VCF
genotypes = vcf.to_genotypes(n_samples=422)
genotypes_std = (genotypes - genotypes.mean()) / genotypes.std()  # Normalize
U, s, Vt = np.linalg.svd(genotypes_std)  # SVD ~ PCA
pc_scores = U[:, :3]  # First 3 PCs
```

### Step 2: Single-SNP association (in carotenoid genes)
```python
# Get SNPs in carotenoid genes
carot_snps = extract_gene_region(vcf, genes=['PSY', 'PDS', 'ZDS', 'Lycopene_cyc'])

results = []
for snp in carot_snps:
    geno = vcf[snp]['genotype']  # 0/1/2 for each sample
    
    # Match to phenotype samples (85 DBVPG)
    idx = [i for i, s in enumerate(vcf.samples) if s in pheno['strain'].values]
    geno_subset = geno[idx]
    pheno_subset = pheno['brightness_L'].values
    
    # Regression: brightness ~ SNP + PC1 + PC2 + PC3
    X = np.column_stack([geno_subset, pc_scores[idx, :3]])
    y = pheno_subset
    b, residuals, rank, s = np.linalg.lstsq(X, y)
    
    # p-value for SNP effect (first coefficient)
    se = np.sqrt(np.sum(residuals**2) / (len(y) - rank))
    t_stat = b[0] / se
    p_val = 2 * (1 - stats.t.cdf(abs(t_stat), len(y) - rank))
    
    results.append({
        'SNP': snp,
        'gene': snp_to_gene(snp),
        'effect': b[0],
        'p_value': p_val
    })

results_df = pd.DataFrame(results)
results_df['q_value'] = adjust_pvalues(results_df['p_value'], method='BH')
```

### Step 3: Visualize
```python
# Manhattan plot
import matplotlib.pyplot as plt
plt.figure()
plt.scatter(range(len(results_df)), -np.log10(results_df['p_value']))
plt.axhline(-np.log10(0.05 / len(results_df)), color='r', linestyle='--')
plt.ylabel('-log10(p-value)')
plt.xlabel('SNP index (ordered by position)')
plt.savefig('manhattan_plot.png')

# Top SNPs
top_snps = results_df.nsmallest(10, 'p_value')
print(top_snps[['SNP', 'gene', 'effect', 'p_value']])
```

---

## Likely Outcomes

### Scenario A: Strong signal (5-10 SNPs, p < 10⁻⁵)
- **Interpretation:** Carotenoid pathway SNPs directly cause brightness variation
- **Next:** Fine-map to identify causal SNP(s)
- **Publication:** Strong association story

### Scenario B: Moderate signal (20-50 SNPs, p < 0.05 after FDR)
- **Interpretation:** Pathway contributes, but effect sizes are small
- **Next:** Aggregate pathway effect; validate in other strains
- **Publication:** Moderate association + pathway story

### Scenario C: Weak/no signal (no SNPs, p > 0.05)
- **Interpretation:** Carotenoid genes not the genetic basis
- **Next:** Consider: cis-regulatory regions, non-coding RNA, epigenetics
- **Alternative:** RNA-seq approach to measure expression

---

## Final Recommendation

**Start with Tier 1 (PCA + LD pruning + annotation + multiple testing correction).**

This is rigorous, standard, and likely sufficient for publication. Only move to Tier 2/3 if Tier 1 results are marginal.

**Expected timeline:** 2 weeks to complete analysis + figure out results

**Expected code volume:** ~300-400 lines Python (manageable)

Would you like me to write the full GWAS pipeline now?

