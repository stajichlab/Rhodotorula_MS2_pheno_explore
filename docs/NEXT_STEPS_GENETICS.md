# Your Next Steps: From Genomes to Brightness Genetics

**Status:** You have 102 matched DBVPG strains ready for analysis  
**Biggest finding:** Gene copy number ≠ Brightness (regulation is likely key)  
**Timeline:** 2-4 weeks to definitive answer

---

## What You Have RIGHT NOW

✅ **Lab brightness phenotypes:** 590 strains, range 61–78 L* (Median_ColorLab_L*Mean)  
✅ **Rodeo genomes:** 272 genomes with Pfam domain annotations  
✅ **DBVPG mapping:** 102 strains with BOTH phenotypes AND reference genomes  
✅ **Top metabolite candidate:** Feature 2755 (m/z 808.51, explains 54% brightness)

---

## Decision Tree: Which Analysis to Pursue?

```
Do you have SNP/VCF data for your lab strains C_1–C_590?
│
├─ YES → Path A: GWAS (SNP-phenotype association)
│         Timeline: 1-2 weeks | Effort: Medium | Power: HIGHEST
│
├─ NO, but can resequence → Path B: Resequencing DBVPG subset
│         Timeline: 4-6 weeks | Effort: High $ | Power: High
│
└─ NO → Path C: RNA-seq expression analysis (recommended)
         Timeline: 2-3 weeks (design) + 6 weeks (experiment) | Effort: High | Power: High
```

**RECOMMENDED:** Start with **Path C (RNA-seq)** if you don't have SNP data. It's:
- More directly tied to your finding (gene regulation, not copy number)
- Measurable at lab scale (culture bright vs. dim strains)
- Directly validates Feature 2755 hypothesis

---

## Path A: GWAS Analysis (If You Have SNP Data)

### Quick Check
Do you have:
- [ ] BAM files (mapped reads) for your C_1–C_590 lab strains?
- [ ] VCF file with SNPs already called?
- [ ] Reference genome assembly for alignment?

### If YES, here's the workflow:

**Week 1: Data Preparation**
```python
# Extract carotenoid gene regions from VCF
# Filter to coding sequences (from GFF3)
# Annotate: synonymous vs. non-synonymous, upstream regions

carot_genes = ['PSY', 'PDS', 'ZDS', 'Lycopene_cyc', 'BCH', 'GGPS']
# Get coordinates from Rodeo annotation

# Extract SNPs in these regions
subset_vcf = vcf[vcf['CHROM:POS in carot_regions']]
```

**Week 2: Association Testing**
```python
# For each SNP in carotenoid genes:
#   brightness ~ genotype (additive model)
#   Correct for population structure (PCA)
#   FDR adjust p-values

# Output: Manhattan plot showing brightness-associated SNPs
```

**Week 3: Interpretation**
- Which SNPs are in promoters (-500 to -1)?
- Which are non-synonymous (loss-of-function)?
- Do they segregate by bright vs. dim strains?

---

## Path C: RNA-seq Expression Analysis (Recommended)

### Why This Approach

Your key finding: **Gene copy number doesn't explain brightness**

**Hypothesis:** High-brightness strains have HIGH expression of carotenoid biosynthesis genes

### Experiment Design (3 replicates each)

**Sample 1: Bright R. mucilaginosa**
- Pick 3 strains with L* > 76 from your lab data
- Examples: C_3 (75.7), C_4 (75.7), C_18 (75.5)
- Culture under standard conditions

**Sample 2: Dim R. mucilaginosa**
- Pick 3 strains with L* < 72 from your lab data
- Examples: C_189 (71.4), C_178 (73.5), C_202 (73.5)
- Same culture conditions

### What to Measure

1. **RNA-seq:** Map reads to R. mucilaginosa reference genome
2. **Quantify expression (TPM):**
   - Phytoene synthase (PSY)
   - Phytoene desaturase (PDS)
   - Zeta-carotene desaturase (ZDS)
   - Lycopene cyclase (LCY)
   - Beta-carotene hydroxylase (BCH)
   - Geranylgeranyl pyrophosphate synthase (GGPS)

3. **qPCR validation:** Top 3 hits

4. **LC-MS metabolomics:** Same strains, measure Feature 2755 abundance

### Expected Output

**Correlation matrix:**
```
           PSY      PDS      ZDS     LCY     Feature_2755   Brightness
PSY        1.0      0.8      0.6     0.5     0.7            0.6
PDS        0.8      1.0      0.7     0.6     0.6            0.5
...
Feature    0.7      0.6      ...     ...     1.0            0.75 ***
Brightness 0.6      0.5      ...     ...     0.75 ***       1.0
```

### Timeline & Cost

| Phase | Time | Cost | Notes |
|-------|------|------|-------|
| 1. Culture & RNA prep | 1 week | $100 | Standard microbiology |
| 2. RNA-seq | 2 weeks | $2–3K | Sequencing facility |
| 3. qPCR validation | 1 week | $500 | Your lab |
| 4. LC-MS metabolites | 1 week | $500 | Your MS facility |
| 5. Analysis & writing | 2 weeks | $0 | Bioinformatics |
| **Total** | **~6 weeks** | **$3–4K** | **Publishable result** |

---

## Path B: Genome Resequencing (If Budget Available)

If you want strain-level SNP data:

1. **Select 20–30 strains** (stratified: bright, medium, dim)
2. **Whole-genome resequencing** (20–50x coverage)
3. **Call SNPs** in carotenoid genes
4. **GWAS** as in Path A

**Cost:** $3–5K (sequencing only) + $1–2K (analysis)  
**Timeline:** 8–12 weeks  
**Power:** Very high (can detect rare variants)

---

## What to Do This Week

1. **Check for existing SNP data:**
   - Do you have VCF files or BAM files for strains?
   - Check `/bigdata/stajichlab/shared/projects/Rhodotorula/` for any vcf/ or bam/ directories

2. **If NO SNP data:**
   - Decide: Resequence or do RNA-seq?
   - If RNA-seq: Pick your 6 strains (3 bright, 3 dim)
   - Contact your sequencing facility about timeline

3. **If YES SNP data:**
   - Get VCF file + R. mucilaginosa reference genome
   - I can write the GWAS pipeline

4. **In parallel:**
   - Expand Pfam analysis to OTHER carotenoid genes (PSY, PDS, ZDS)
   - Script already exists; just need to filter Rodeo Pfam table

---

## Files Created

- `BRIGHTNESS_GENETICS_STRATEGY.md` — Full brainstorm of 5 research angles
- `CAROTENOID_GENE_BRIGHTNESS_ANALYSIS.md` — Initial gene copy number analysis
- `DATA_DESCRIPTION.md` — Complete data structure reference
- `FEATURE_ANALYSIS.md` — Top 5 metabolite features (Feature 2755 deep-dive)
- Visualizations: `metabolite_features_visualization.html`

---

## Key Question for You

**Which path are you taking?**

- [ ] **Path A:** I have SNP data (VCF, BAM)
- [ ] **Path B:** I want to resequence selected strains
- [ ] **Path C:** I want to do RNA-seq + metabolomics (recommended if no SNP data)
- [ ] **Path D:** I want more time to decide

Once you answer, I can:
1. Write the exact analysis scripts
2. Design the RNA-seq experiment (if Path C)
3. Set up the GWAS pipeline (if Path A)

---

## Summary: Why You're Well-Positioned

1. **Large sample size:** 590 strains with brightness measurements
2. **Reference genomes:** 272 annotated genomes to compare against
3. **Matched DBVPG subset:** 102 strains with phenotypes + genomes
4. **Strong metabolite signal:** Feature 2755 explains 54% of variance
5. **Clear hypothesis:** Gene regulation (not copy number) drives brightness

**This is publishable work.** Let's just pick the execution path and go.

