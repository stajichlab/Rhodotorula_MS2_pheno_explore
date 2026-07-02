# Strategy: Linking MS2 Metabolite Features to Color Phenotypes
## Rhodotorula Metabolites Project

**Objective:** Identify which metabolites (MS2 features) correlate with color phenotypes (HSV, CIELab)

---

## Data Structure Summary

### Available Data

**1. MS2 Metabolite Features**
- File: `Rhodotorula_MS2_aligned_features_ms2.csv`
- **16,332 features** (m/z × retention time aligned across all samples)
- **608 samples** (columns named: C_1...C_334, SUP_1...SUP_99, etc.)
- Peak area abundance for each feature in each sample

**2. Phenotype Data (already analyzed)**
- 584 samples with color phenotypes
- Key traits: Brightness, Hue, Saturation, L*, a*, b*, Chroma
- 96.1% match with strain trait summary

**3. Metadata**
- Sample names, species, strain, plate information
- Already merged into extended metadata file

### Integration Challenge

MS2 metabolite sample names (e.g., "C_1.mzML Peak area") need to be mapped to metadata sample names (e.g., "C_1" from filename column).

---

## Proposed Analysis Pipeline

### PHASE 1: Data Preparation & QC

**1.1 Load and harmonize data**
```
Input: MS2 features (16,332 × 608)
       Phenotype data (584 samples × 7 color traits)

Steps:
  • Extract sample columns from MS2 data
  • Parse sample IDs (remove ".mzML Peak area" suffix)
  • Match with metadata file samples
  • Create aligned matrix: 16,332 features × 584 phenotype samples
```

**1.2 Feature filtering (reduce false correlations)**
- **Remove low-abundance features:**
  - Features with median peak area = 0 across all samples (noise)
  - Features present in <10% of samples (unreliable)
  - Rationale: Only analyze meaningful metabolites
  
- **Remove low-variance features:**
  - Features with very low variance (no signal)
  - Keep those with CV > 0.1 or IQR-based filtering
  - Rationale: No correlation possible without variation

- **Expected outcome:** ~3,000-8,000 features (remove ~50-80% noise)

**1.3 Data transformation**
- **Log2 transform abundances:** log2(peak_area + 1)
  - Normalizes skewed distributions
  - Makes correlations more robust
  - Rationale: Peak areas often log-normal
  
- **Standardize (optional):** Z-score normalization
  - Ensures comparable scales
  - Helpful for visualization

---

### PHASE 2: Correlation Analysis

**2.1 Compute MS2 ↔ Phenotype correlations**

For each of ~5,000 filtered metabolite features:
- Compute Pearson correlation with each color phenotype (7 traits)
- Calculate: r-value, p-value, 95% CI

**Correlations to compute (7 target phenotypes):**
1. Brightness (HSV_BrightnessMean) - strongest predictor from prior analysis
2. Hue (HSV_HueMean)
3. Saturation (HSV_SaturationMean)
4. L* (CIELab_L* - lightness)
5. a* (CIELab_a* - red-green axis)
6. b* (CIELab_b* - yellow-blue axis)
7. Chroma (CIELab_ChromaEstimatedMean)

**Output:** 5,000 features × 7 phenotypes correlation matrix

**2.2 Multiple testing correction**
- Apply FDR correction (Benjamini-Hochberg) across all ~35,000 tests
- Set q-value threshold: 0.05 (5% expected false discovery rate)
- Rationale: Controls for multiple comparisons without being overly conservative

**2.3 Identify significant associations**
- **Strict threshold:** |r| > 0.3 AND q < 0.05
- **Exploratory threshold:** |r| > 0.2 AND q < 0.05
- Rationale: r=0.3 explains ~9% of variance (meaningful biological effect)

---

### PHASE 3: Feature-Level Analysis

**3.1 Prioritize features by robustness**
- Features correlated with **multiple phenotypes** → likely genuine signals
- Features with high |r| values → strong effect size
- Features with multiple metabolite identifications → better confidence

**3.2 Characterize feature properties**
For high-priority features, extract:
- m/z (mass-to-charge ratio)
- Retention time
- Charge state
- Adduct type (±H, ±Na, etc.)
- MS/MS data (if present)
- Abundance distribution across samples

**3.3 Group features by correlation pattern**
- **Brightness-associated:** Features positively/negatively correlated with brightness
- **Hue-associated:** Features varying with red-yellow-blue spectrum
- **Saturation-associated:** Features correlating with color intensity
- **Species-specific metabolites:** Features enriched in certain Rhodotorula species

---

### PHASE 4: Visualization & Interpretation

**4.1 Create correlation plots**

```
Plot 1: Heatmap of top metabolite features × phenotypes
- Show features with |r| > 0.25
- Hierarchical clustering: identify co-correlated features
- Color: red (positive) to blue (negative)

Plot 2: Scatter plots (top hits)
- Feature abundance vs Brightness (strongest phenotype)
- Feature abundance vs Hue (second strongest)
- Include trend lines and R² values

Plot 3: Volcano plot (for each phenotype)
- X-axis: correlation coefficient (r)
- Y-axis: -log10(q-value)
- Color by significance threshold
- Identify steep slope features (strong + reliable)

Plot 4: Feature abundance distribution by species
- Box plots of top features grouped by species
- Show if metabolites drive species color differences
```

**4.2 Metabolite identification**
- Match features to known Rhodotorula metabolites
- If MS/MS available: compare to spectral libraries
- Link to: carotenoids, melanins, other pigments

**4.3 Biological interpretation**
- What metabolites control brightness? (Likely: carotenoid abundance)
- What drives hue variation? (Yellow vs red pigments)
- Are color phenotypes determined by single or multiple metabolites?

---

## Expected Outcomes

### Positive Scenario (Good hits found)
- **~50-200 metabolites** with significant phenotype correlations
- **Top 5-10 features** strongly associated with brightness
- Clear mechanistic insights:
  - "High-brightness strains have elevated m/z 550 metabolite"
  - "Red hue strains enriched in m/z 234 feature"
  - "Species-specific metabolite explains color differences"

### Realistic Scenario (Complex relationships)
- **~500-1000 features** pass initial correlation threshold
- Multiple weak correlations (no single dominant metabolite)
- Suggests:
  - Color is polygenic/polygenic-equivalent (many small effects)
  - Metabolite interactions (not individual effects)
  - Non-linear relationships (need different approach)

### Null Scenario (No strong hits)
- Few features survive multiple testing correction
- Suggests:
  - Color phenotypes driven by non-metabolite factors (cell wall structure, pigment organization)
  - Metabolite effects masked by sample variation
  - Need higher resolution data or alternative approaches

---

## Statistical Considerations

### Sample Size
- N = 584 samples with color phenotypes
- **Sufficient for detecting r > 0.12 at α = 0.05**
- Minimum meaningful effect: r = 0.2-0.3 (explains 4-9% variance)

### Multiple Testing Burden
- ~5,000 features × 7 phenotypes = 35,000 tests
- FDR correction necessary
- Expected false positives at q < 0.05: ~1,750
- Expected true positives: need |r| > 0.25-0.3 for confidence

### Missing Data Handling
- MS2: May have zeros (not detected) vs NAs (not measured)
- Phenotype: 24 missing samples
- Strategy: Use pairwise complete observations (include sample if both MS2 AND phenotype present)

---

## Implementation Roadmap

### Step 1: Data Integration (1-2 hours)
- Load MS2 + phenotype data
- Map sample IDs
- Create aligned matrix
- QC report on missingness

### Step 2: Feature Filtering (30 mins)
- Remove zero/low-abundance features
- Apply CV threshold
- Document filtering decisions

### Step 3: Correlation Analysis (1-2 hours)
- Compute correlations
- FDR correction
- Identify significant features

### Step 4: Visualization (2-3 hours)
- Create summary plots
- Generate top-feature scatter plots
- Produce heatmaps

### Step 5: Interpretation (Variable)
- Metabolite matching
- Species-level analysis
- Biological write-up

**Total estimated effort:** 6-10 hours for complete analysis

---

## Advanced Options (if initial analysis successful)

### If strong hits identified:
1. **Partial correlation:** Control for species/plate effects
2. **Pathway analysis:** Link metabolites to biosynthetic pathways
3. **Interaction terms:** Do metabolites work synergistically?
4. **Mediation analysis:** Do metabolites mediate species → color relationships?

### If weak/no hits:
1. **Alternative features:** 
   - Use feature ratios (metabolite ratioing often more informative)
   - Focus on MS/MS fragmentation patterns (more selective)
   
2. **Different scaling:**
   - Arcsin transformation (for compositional data)
   - Robust standardization (resistant to outliers)
   
3. **Non-linear models:**
   - Spearman rank correlation (robust to outliers)
   - Local regression (LOWESS) for non-linear relationships

---

## Quality Assurance Checkpoints

- [ ] Sample ID mapping: 100% of samples matched
- [ ] Feature filtering: Document why features removed
- [ ] Correlation matrix: Check for unexpected values (all -1 to +1?)
- [ ] Multiple testing: p-values appear uniform under null?
- [ ] Top hits: Do they make biological sense?
- [ ] Reproducibility: Can rerun analysis with same results?

---

## Deliverables

1. **Filtered MS2 dataset** (features × samples, cleaned)
2. **Correlation matrix** (features × phenotypes)
3. **Top features table** (m/z, RT, r, p-value, q-value for each phenotype)
4. **Publication-quality plots** (heatmaps, volcanoes, scatter plots)
5. **Interpretation document** (biological findings + mechanisms)
6. **R/Python code** (reproducible, documented)

---

## Next Decision Point

**Would you like me to:**

1. **Implement the full pipeline** (phases 1-4)?
2. **Start with phase 1-2 only** (correlation analysis)?
3. **Focus on specific phenotypes** (e.g., brightness only)?
4. **Use alternative statistical approach** (e.g., machine learning)?

**Or refine the strategy first?** Key questions:
- Should we prioritize features with MS/MS data?
- Are you interested in species-specific metabolites?
- Do you want to control for plate/batch effects?
