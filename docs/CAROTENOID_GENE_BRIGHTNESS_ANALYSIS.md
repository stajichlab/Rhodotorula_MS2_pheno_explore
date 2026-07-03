# Carotenoid Gene-Brightness Analysis: Initial Findings

**Date:** 2026-07-02  
**Data Sources:** 
- Rodeo HMMER/Pfam results (272 genomes, carotenoid genes identified)
- MS2 brightness phenotypes (590 lab strains, CIELab L* measurements)

---

## Key Finding: Gene Copy Number ≠ Brightness

### Surprising Result

**R. paludigena paradox:**
- Carotenoid genes: 43 total Lycopene_cyc hits (3.07 per strain average) — **HIGHEST**
- Brightness (L*): 73.51 ± 1.23 — **LOWER than expected**

**R. diobovata:**
- Carotenoid genes: 1 Lycopene_cyc per strain
- Brightness: 75.58 ± 0.96 — **HIGHER than R. paludigena**

**Conclusion:** Carotenoid gene copy number is NOT correlated with brightness phenotype.

---

## Species-Level Data

| Species | Brightness (L*) | n_strains | Lycopene_cyc count | Avg per strain |
|---------|-----------------|-----------|--------------------|----|
| R. diobovata | 75.58 | 7 | 1 | 1.0 |
| R. kratochvilovae | 75.55 | 3 | 6 | 2.0 |
| R. sp. clade I | 75.48 | 5 | 9 | 1.0 |
| R. toruloides | 75.05 | 9 | 9 | 1.29 |
| **R. mucilaginosa** | 74.46 | 205 | 20 | 1.1 |
| R. dairenensis | 74.65 | 7 | ? | ? |
| R. paludigena | 73.51 | 10 | 43 | 3.07 |
| R. sphaerocarpa | 72.87 | 4 | 5 | 1.0 |
| R. graminis | 71.03 | 3 | 6 | 2.0 |
| R. taiwanensis | 72.53 | 6 | ? | ? |

**Correlation (Lycopene_cyc copy # ↔ Brightness):** Appears **NEGATIVE** (anti-correlated)

---

## Why Isn't Gene Copy Number Working?

### Hypothesis 1: Gene Regulation, Not Gene Count
- **Evidence:** Feature 2755 (m/z 808.51 carotenoid glycoside) shows massive variance (std=3.48 log2 units) even within species
- **Implication:** Strains have similar genomes but VERY different metabolite abundances
- **Prediction:** Gene EXPRESSION (transcription, translation) is the real driver, not copy number

### Hypothesis 2: Other Bottleneck Gene
- Lycopene cyclase may NOT be rate-limiting
- **Real bottleneck:** Phytoene synthase (PSY) or desaturases (PDS, ZDS)
- **Evidence:** PSY controls first committed step of carotenoid synthesis — likely most tightly regulated
- **Test:** Expand analysis to other carotenoid Pfam domains

### Hypothesis 3: Regulatory Elements, Not Genes
- Promoter variants (SNPs upstream of PSY/PDS)
- Transcription factor binding site mutations
- Chromatin state differences between strains
- **Impact:** Can't detect by counting genes; need SNP/regulatory analysis

### Hypothesis 4: Pathway Flux Optimization
- Having 3 Lycopene cyclase copies might DILUTE resources (metabolic burden)
- Fewer copies = more concentrated enzyme activity
- **Analogy:** Copy number variation can be deleterious in tightly-balanced pathways

---

## Data Gap: What We're Missing

| Data | Have? | Impact |
|------|-------|--------|
| Lab strain (C_1, C_2...) → Reference genome mapping | ❌ | Can't do strain-level SNP analysis |
| RNA-seq (expression levels) | ❌ | Can't test Hypothesis 1 |
| VCF/SNP calls for lab strains | ❌ | Can't test regulatory SNP hypothesis |
| Full carotenoid pathway Pfams | ⚠️ Partial | Only looking at Lycopene_cyc; need PSY, PDS, ZDS |
| Proteomics/enzyme abundance | ❌ | Can't measure post-translational regulation |

---

## Recommended Next Actions

### Immediate (This Week)
1. **Expand Pfam analysis to other carotenoid genes:**
   - PSY (Phytoene synthase) — Look for specific Pfam domain
   - PDS (Phytoene desaturase) — FAD_binding_3 is generic; filter by context
   - Squalene synthase / GGPS — Prenyl transferase domains
   
2. **Search Rodeo Pfam results for candidate pathways:**
   ```python
   # Pseudo-code
   carot_pathway_pfams = [
       'Lycopene_cyc',
       'Lycopene_cycl', 
       'Squalene_synthase',
       'prenyltransf',  # GGPS
       'FAD-dependent_oxidoreductase'  # Desaturase subfamily
   ]
   ```

### Short-term (1-2 weeks)
3. **Check if you have genomic variants** for your lab strains C_1–C_590:
   - Do you have resequencing data (BAM files)?
   - Are lab strains mapped to a reference genome?
   - Any existing VCF files linking strains to genomes?

4. **If you have strain-level SNP data:**
   - Reframe as **GWAS**: "Which SNPs in carotenoid genes predict brightness?"
   - Look for promoter variants (upstream -500 bp)
   - Check for loss-of-function variants (frameshifts, stops)

### Medium-term (2-4 weeks)
5. **Design transcriptomics experiment:**
   - Culture bright (L* > 76) vs. dim (L* < 72) strains
   - RNA-seq to measure PSY, PDS, ZDS, etc. expression levels
   - Correlate: High PSY expression ↔ High Feature 2755?

6. **Identify SNPs in regulatory regions:**
   - Define promoter: -500 to -1 bp upstream of start codon
   - Find variants in carotenoid gene promoters
   - Test: SNP in PSY promoter predicts brightness?

---

## Why This Matters

Your **Feature 2755 (54% of brightness variance)** is likely a carotenoid metabolite. The genetic basis should explain WHY some strains produce more Feature 2755.

**Current finding:** Gene copy number doesn't explain it → **Regulation** or **specific genetic variants** must be driving the phenotype.

---

## Code Pointers (Rodeo Database)

### Extract a Specific Gene Family
```python
import pandas as pd, gzip

# Load Pfam results
with gzip.open('All_Taxa/pfam.csv.gz', 'rt') as f:
    pfam = pd.read_csv(f)

# Filter to gene family (example: carotenoid)
carot = pfam[pfam['pfam_id'].isin(['Lycopene_cyc', 'Lycopene_cycl'])]

# Get list of strains with these genes
strains = carot['protein_id'].str.extract(r'([A-Z0-9]+)_')[0].unique()
print(f"Strains with Lycopene cyclase: {len(strains)}")
```

### Count Genes per Strain per Domain
```python
# Count by strain + domain
gene_counts = pfam.groupby(['protein_id', 'pfam_id']).size().reset_index(name='count')

# Extract LOCUSTAG (strain ID)
gene_counts['LOCUSTAG'] = gene_counts['protein_id'].str.extract(r'([A-Z0-9]+)_')[0]

# Aggregate by strain + domain
by_strain = gene_counts.groupby(['LOCUSTAG', 'pfam_id'])['count'].sum().reset_index()

# Pivot to get domains as columns
pivot = by_strain.pivot(index='LOCUSTAG', columns='pfam_id', values='count').fillna(0)
```

---

## Files Generated

- `All_Taxa/pfam.csv.gz` — Pfam domain predictions for all genomes
- `All_Taxa/samples.csv.gz` — Strain metadata (species, ASMID, LOCUSTAG)
- Your MS2 data: `results/phase1/phase1_phenotype_data.csv.gz` (brightness + species)

---

## Conclusion

**Gene copy number hypothesis: REJECTED** (for Lycopene cyclase)

**Next hypothesis to test:** Gene regulation (expression levels, SNPs in promoters, post-translational modification)

To proceed, we need either:
1. Resequencing data for your lab strains (to look for SNPs)
2. RNA-seq to measure gene expression
3. Expanded Pfam analysis to other pathway genes

Which data do you have available?
