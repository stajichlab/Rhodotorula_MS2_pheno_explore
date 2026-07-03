# Research Angles: Connecting Genetics to Brightness Production

**Goal:** Map genetic variation → metabolite features → brightness phenotype  
**Starting Point:** Feature 2755 (m/z 808.51) explains 54% of brightness; likely carotenoid glycoside  
**Available Data:** 20+ annotated genomes, 590 strain brightness measurements, 7,341 filtered metabolites

---

## Tier 1: Direct Pathway Mapping (Highest ROI, 1-2 weeks)

### 1.1 **Carotenoid Biosynthesis Gene Catalog**

**Goal:** Find all carotenoid pathway genes in your genomes; link to brightness.

**Method:**
1. Search *R. mucilaginosa* reference genome (most complete annotation) for known carotenoid pathway genes:
   - Phytoene synthase (PSY)
   - Phytoene desaturase (PDS)
   - ζ-carotene desaturase (ZDS)
   - Lycopene cyclase (LCY)
   - β-carotene hydroxylase (BCH)
   - Geranylgeranyl pyrophosphate synthase (GGPS)

2. Extract protein sequences from `Rhodotorula_mucilaginosa_DBVPG_3045_hap1.proteins.fa`

3. Query HMMER/BLAST against all 20+ species genomes → get orthologs

4. Compare brightness in strains with:
   - Intact pathway (all genes present, single copy)
   - Gene duplications (copy number variants)
   - Gene loss/pseudogenization (frameshifts, stops)

**Output:** Table of brightness vs. carotenoid gene completeness → correlation?

**Tool Stack:** HMMER, BLAST, bedtools, bcftools

**Feasibility:** ⭐⭐⭐⭐⭐ (High) — Standard genomics, fast turnaround

---

### 1.2 **Orthologue Conservation vs. Brightness**

**Goal:** Do bright strains have conserved carotenoid genes; dim strains have divergent/lost?

**Method:**
1. Build OrthoFinder pangenome (20+ species) → clusters of orthologous genes
2. Extract brightness-associated orthogroups (carotenoid pathway)
3. Measure sequence divergence (Ka/Ks) within orthogroups across species
4. Correlate: High conservation in bright species ↔ Low conservation in dim?

**Hypothesis:** Bright species = strong purifying selection on carotenoid genes (essential for color); dim species = relaxed selection (lost function)?

**Output:** Phylogenetic tree with brightness colored at tips; branch-specific Ka/Ks for pathway genes

**Feasibility:** ⭐⭐⭐⭐ (High) — Established methodology, depends on pangenome availability

---

## Tier 2: Genomic Variation Analysis (2-4 weeks)

### 2.1 **SNP-Phenotype Association (GWAS Analog)**

**Goal:** Find carotenoid pathway SNPs associated with brightness differences *within* species.

**Method:**
1. Subset to *R. mucilaginosa* strains (largest group, 200+ samples)
2. Call SNPs from whole-genome sequences (if available) or map reads to reference
3. Extract variants in carotenoid pathway regions
4. Test: Does each SNP associate with brightness?
   - Linear regression: brightness ~ genotype (additive model)
   - Corrected for population structure (PCA covariates)
5. Manhattan plot: Show significant SNPs in pathway genes

**Key SNP Classes to Track:**
- Synonymous (neutral) vs. non-synonymous (functional)
- Promoter variants (upstream -500 to -1)
- Stop-gain/frameshift (loss-of-function)

**Output:** List of brightness-associated SNPs in carotenoid genes; effect sizes

**Feasibility:** ⭐⭐⭐ (Medium) — Requires resequencing data; standard pipeline if available

---

### 2.2 **Copy Number Variation (CNV) & Brightness**

**Goal:** Do bright strains have duplicated carotenoid genes?

**Method:**
1. Map short reads (if available) to reference genome
2. Call CNVs using depth-of-coverage (samtools/CNVnator)
3. Identify CNVs overlapping carotenoid genes
4. Test: Do strains with carotenoid gene duplications have higher brightness?
5. Compare:
   - 1 copy (baseline)
   - 2+ copies (duplication)
   - 0 copies (deletion)
6. Effect size: brightness difference per additional copy?

**Hypothesis:** Extra carotenoid synthase copies → higher carotenoid production → brighter

**Output:** CNV map; bar plot of brightness by copy number for key genes

**Feasibility:** ⭐⭐⭐ (Medium) — Easier if short-read data exists; harder if only scaffolds

---

## Tier 3: Gene Expression Integration (4-8 weeks)

### 3.1 **Transcriptomic Association: Carotenoid Gene Expression ↔ Metabolites**

**Goal:** Do strains with high carotenoid gene expression produce more Feature 2755?

**Method:**
1. If you have RNA-seq: Map reads to carotenoid pathway genes
2. Quantify expression (TPM) for each gene in each strain
3. Correlate: Expression level ↔ Feature 2755 abundance
4. Test: High carotenoid expression strains have high Feature 2755?

**Design (if no RNA-seq yet):**
- Culture bright vs. dim *R. mucilaginosa* strains under standard conditions
- RNA-seq (3 replicates each)
- qPCR validation for top hits
- Measure metabolites by LC-MS (your Feature 2755)

**Output:** Expression vs. metabolite scatter plot; Spearman ρ per gene

**Feasibility:** ⭐⭐ (Low) — Requires new RNA-seq experiment; 4-6 week timeline

---

### 3.2 **Proteomics: Carotenoid Enzyme Abundance ↔ Brightness**

**Goal:** Do bright strains accumulate more carotenoid synthase protein?

**Method:**
1. If you have proteomics data: Quantify PSY, PDS, other pathway enzymes
2. Correlate: Protein level ↔ brightness + Feature 2755
3. Check for post-translational modification (phosphorylation, ubiquitination) patterns

**Hypothesis:** Bright strains = high basal enzyme levels OR post-translationally active enzymes

**Feasibility:** ⭐⭐ (Low) — Requires proteomics; expensive/specialized

---

## Tier 4: Functional Validation (8-16 weeks)

### 4.1 **CRISPR Knockout: Carotenoid Genes → Brightness**

**Goal:** Proof-of-concept: Delete PSY or PDS in bright strain → becomes dim?

**Method:**
1. Design CRISPR guides targeting carotenoid synthase in *R. mucilaginosa*
2. Generate knockouts
3. Measure:
   - Brightness (L* phenotype)
   - Feature 2755 abundance (LC-MS)
   - Carotenoid levels (HPLC)
4. Reverse test: Overexpress carotenoid gene → brighter?

**Output:** Knockout phenotype; validation that pathway controls brightness

**Feasibility:** ⭐ (Very Low) — Requires CRISPR development; 12+ week timeline; needs biosafety approval

---

## Tier 5: Comparative Genomics (2-3 weeks)

### 5.1 **Bright vs. Dim Lineage Comparison**

**Goal:** Evolutionary perspective: Did bright and dim *Rhodotorula* lineages diverge in carotenoid pathway?

**Method:**
1. Root phylogenetic tree on all 20+ genomes
2. Identify nodes separating "bright" clades from "dim" clades
3. Extract genes that:
   - Duplicated in bright clade → lost in dim clade
   - Show positive selection (dN/dS > 1) in bright clade
   - Are branch-length specific to bright lineage

**Output:** Phylogeny with branch-colored by carotenoid gene copy number; identify "brightness acquisition" events

**Feasibility:** ⭐⭐⭐⭐ (High) — Comparative genomics standard approach

---

### 5.2 **Horizontal Gene Transfer (HGT) Screen**

**Goal:** Did carotenoid genes come from non-*Rhodotorula* source?

**Method:**
1. For each carotenoid pathway gene in *R. mucilaginosa*:
   - BLAST against NCBI NR database
   - Check top hits: all from *Rhodotorula* or from other fungi/bacteria?
2. If HGT detected:
   - Date acquisition (phylogenetic reconciliation)
   - Correlate: Did species acquiring HGT become brighter?

**Output:** Evidence of HGT events; timeline of acquisition

**Feasibility:** ⭐⭐⭐ (Medium) — Bioinformatic, but phylogenetic inference is complex

---

## Recommended Prioritization (What to do first)

### **Week 1-2: Quick Wins**
1. ✅ **1.1 Carotenoid Gene Catalog** — Identify all pathway genes; map to your genomes
2. ✅ **5.1 Bright vs. Dim Phylogeny** — Tree showing brightness coloration + carotenoid gene status

**Outcome:** Know which genes to focus on; understand evolutionary context

### **Week 3-6: Core Analysis**
3. ✅ **1.2 Orthologue Conservation** — Ka/Ks analysis; is bright = conserved?
4. ✅ **2.1 SNP-Phenotype Association** (if sequence data exists) — Find causal variants

**Outcome:** Genetic mechanisms (conservation vs. variation); specific SNPs to follow up

### **Week 7+: Deep Dives**
5. ⚠️ **2.2 CNV Analysis** — If copy number data available
6. ⚠️ **3.1 Transcriptomics** — Design/run if not already collected

---

## Data You Need (Checklist)

- [ ] Annotated genomes for all strains (GFF3 with gene names) ✅ **Have**
- [ ] Protein sequences (FASTA) ✅ **Have**
- [ ] Genome assembly quality metrics (N50, completeness)
- [ ] Short reads or BAM files (for SNP/CNV calling) — **Unknown, check**
- [ ] RNA-seq (if considering expression analysis) — **Unknown**
- [ ] Resequencing data for multiple strains — **Unknown**
- [ ] Strain assignment to species (to subset for within-species analysis) ✅ **Have**

**Next Step:** Audit what sequence data exists; prioritize Tier 1 approaches while you prepare Tier 2+3

---

## Bridge to Your Metabolites

**Key Link:** Feature 2755 (m/z 808.51 carotenoid glycoside) is your **proxy for pathway activity**.

Instead of just brightness, you can now ask:
- **Does PSY gene copy number correlate with Feature 2755 abundance?** (proximate metric)
- **Which SNPs in carotenoid genes predict high Feature 2755?** (functional prediction)
- **Did *R. mucilaginosa* acquire a carotenoid duplication event that dim species missed?** (evolutionary hypothesis)

---

## Tools & Scripts to Build

| Analysis | Tool | Input | Output |
|----------|------|-------|--------|
| Gene catalog | HMMER/BLAST | Genomes + seed proteins | `carotenoid_genes_all_species.tsv` |
| Orthologue tree | OrthoFinder | All proteomes | `Orthogroups.tsv` + phylogeny |
| SNP calling | BCFtools/GATK | BAM files | `carotenoid_pathway.vcf` |
| SNP-phenotype test | R/statsmodels | VCF + brightness data | `gwas_results.txt` |
| CNV detection | CNVnator/LUMPY | BAM files | `cnv_calls.bedpe` |
| Ka/Ks | PAML/codeml | Alignments | `branch_dn_ds.txt` |

---

## Estimated Timeline & Resource Needs

| Approach | Timeline | Compute | Cost | High Risk? |
|----------|----------|---------|------|-----------|
| 1.1 Gene catalog | 1 week | Desktop | $0 | No |
| 1.2 Orthologue conservation | 2 weeks | 4-core server | $0 | No |
| 2.1 SNP-phenotype | 2 weeks | 8-core server | $100–500 (if sequences to download) | No |
| 2.2 CNV | 1 week | 4-core server | $0–200 | No |
| 3.1 Transcriptomics | 4–6 weeks | sequencing facility | $3–5K | Yes (experimental) |
| 4.1 CRISPR validation | 12+ weeks | lab time | $10–30K | Yes (technical) |

---

## Next Actions

1. **Audit existing data:** Do you have resequencing data (BAM/VCF) for your strains?
2. **Select entry point:** Start with 1.1 (gene catalog) this week
3. **Build script:** Script to extract carotenoid pathway genes from GFF3 → map to all species
4. **Validate connection:** Confirm Feature 2755 is indeed carotenoid (MS2 fragmentation, standards)

Would you like me to:
- [ ] Script the gene catalog extraction (HMMER search)?
- [ ] Build the orthologue phylogeny with OrthoFinder?
- [ ] Check what sequence data you have available?
- [ ] Design the SNP-phenotype association test?
