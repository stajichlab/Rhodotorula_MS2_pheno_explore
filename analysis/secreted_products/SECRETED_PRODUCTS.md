# Secreted Products — uniquely & highly secreted metabolites → candidate genes

**Status:** complete (v1) · **Date:** 2026-07-02 · **Author:** Jason Stajich (+ Claude)
**Script:** `scripts/01_secretion_analysis.py` · **Reproduce:** `bash run.sh`

## Question

Identify MS2 metabolite products that are **uniquely and highly secreted** across
*Rhodotorula* strains, and suggest **candidate biosynthetic gene families** to look for in
the genome datasets to explain the differences.

## Key idea (what earlier phases missed)

The raw MS2 table contains **paired cell (`C_*`) and supernatant (`SUP_*`) samples** for
each strain. The phase 0–3 pipeline collapsed this, but it is exactly the axis that
defines *secretion*. Here secretion is a **within-strain, paired SUP-vs-C contrast**, which
also cancels much of the strain-level and plate-level nuisance variation.

## Method

- **Samples:** 590 MS2 columns matched to metadata; restricted to strains with **both** a
  cell and a supernatant sample → **292 paired strains** (**270 also have a genome link**).
- **Normalization:** within-sample total-ion-current (relative abundance, columns sum to 1);
  pseudocount `EPS = 1e-6`.
- **Secretion score:** per feature, per strain `log2((SUP+EPS)/(C+EPS))`; feature-level
  **median log2FC** across strains + **paired Wilcoxon signed-rank** (SUP vs C), **BH-FDR**.
- **Secreted** = median log2FC ≥ **1.0** (≥2×) **and** q < **0.05** **and** detected in ≥ **10%**
  of supernatant samples.
- **Uniqueness:** cross-species **Tau specificity** of mean supernatant abundance over the
  17 labeled species; **uniquely secreted** = secreted **and** Tau ≥ **0.85**.
- **Robustness:** sign-flip **permutation null** (shuffle C/SUP label within strain, 200×)
  and a **sensitivity sweep** over log2FC ∈ {0.5,1,1.5,2} × Tau ∈ {0.75,0.85,0.95}.
- Seed 20260702; every step asserts shapes/counts and logs row counts (robust-analysis).

## Results (headline numbers → `outputs/numbers.json`)

| Quantity | Value |
|---|---|
| Paired strains (C+SUP) | **292** |
| …that also have a genome link | **270** |
| Features tested | 16,332 |
| Species | 17 |
| **Secreted** features (q<0.05, ≥2×, ≥10% detect) | **6,724** |
| **Uniquely secreted** (Tau≥0.85) | **130** (all have MS2) |

- Top uniquely-secreted features are **lipid-like** (recurring neutral masses ~782/787/799/
  1161 Da as multiple adducts & charge states, RT ~1–2 min) — consistent with **secreted
  glycolipid / polyol-lipid biosurfactants**. A polar tail (m/z 160–300, RT 4–5) is a
  candidate **siderophore / peptide** set (e.g. rhodotorulic acid).
- Species structure is strong: 44/130 top in *R. paludigena*, then *sphaerocarpa,
  taiwanensis, dairenensis, pacifica* — ideal for presence/absence comparative genomics.

## Robustness & honest caveats

- **Null model:** the *effect-size-only* gate (median log2FC ≥ 1) yields 6,732 features vs a
  permutation-null mean of **2,443** → empirical FDR ≈ **0.36** for effect size *alone*.
  This is precisely why the reported **secreted** set *also* requires the paired-Wilcoxon
  **BH q < 0.05** (a valid sign-flip test that controls FDR at 5%) and a detection floor.
  The **uniquely-secreted** set (Tau≥0.85) is the high-confidence, actionable subset.
- **Compositional normalization caveat:** TIC normalization makes abundances relative, so a
  few dominant supernatant compounds can mildly deflate others' ratios. The paired design
  mitigates this; a non-compositional re-check is listed as future work.
- **Class assignments are mass/RT hypotheses**, not structures — MS2 annotation (SIRIUS/GNPS)
  is required before committing.
- Sensitivity table shows the unique-set size is threshold-sensitive on Tau
  (25 at 0.95 → 431 at 0.75 with the loose log2FC), so we report the mid, defensible cut.

## Candidate genes & required data

See **`outputs/candidate_gene_families.md`** for the full mapping. In brief, screen genomes
for: glycolipid/polyol-lipid clusters (glycosyltransferase + acyltransferase + FAS/elongase
+ MFS exporter + CYP450), NRPS/NRPS-like siderophore genes (rhodotorulic acid `SidA/SidD`),
carotenoid/apocarotenoid genes (`crtE/crtYB/crtI` + CCD), PKS, and tailoring/export enzymes.
Data needed: assembled + annotated genomes (have most; strain↔genome key already in
metadata), **antiSMASH** BGC calls, **OrthoFinder** presence/absence, a phylogenomic tree
for correction, **SIRIUS/GNPS** MS2 annotation, and a **presence-absence / metabolite-GWAS**
of BGC/ortholog presence vs per-strain supernatant abundance across the 270 strains.

## Outputs

| File | Contents |
|---|---|
| `outputs/secretion_scores_all_features.csv.gz` | All 16,332 features with scores, q, Tau, flags, chemistry |
| `outputs/secreted_features.csv` | 6,724 secreted features |
| `outputs/uniquely_secreted_features.csv` | 130 uniquely-secreted candidates (ranked) |
| `outputs/sensitivity_thresholds.csv` | threshold sweep |
| `outputs/null_distribution.csv` | permutation null counts |
| `outputs/secretion_summary.json`, `outputs/numbers.json` | headline values |
| `outputs/figures/volcano_secretion.png` | secretion volcano |
| `outputs/figures/secretion_score_hist.png` | score distribution |
| `outputs/figures/top_uniquely_secreted.png` | top-25 candidates |
| `outputs/candidate_gene_families.md` | gene-family shortlist + genome-data plan |

## Open questions / next steps

1. Run SIRIUS/GNPS on the top uniquely-secreted MS2 spectra → firm compound classes.
2. antiSMASH + OrthoFinder on the 270 genomes; build strain↔genome↔metabolite join.
3. Presence-absence / metabolite-GWAS of BGCs vs per-strain supernatant abundance.
4. Non-compositional normalization re-check of the secretion scores (robustness).
