---
topic: secreted-metabolites
description: Which metabolites microbes export to the supernatant vs retain intracellularly, and how secretion is measured from paired cell/supernatant MS data.
created: 2026-07-02
last_updated: 2026-07-02
status: active
---

# Secreted metabolites (supernatant-vs-cell)

## F-001: A paired supernatant-vs-cell contrast identifies ~130 uniquely & highly secreted Rhodotorula metabolites
**Status:** preliminary
**Claim:** Using paired cell (`C_*`) and supernatant (`SUP_*`) MS2 samples across 292
*Rhodotorula* strains, 6,724 of 16,332 features are significantly supernatant-enriched
(median log2(SUP/C)≥1, paired-Wilcoxon BH q<0.05, ≥10% SUP detection), and 130 of these are
also species-specific (cross-species Tau≥0.85). The top uniquely-secreted features are
lipid-like (recurring neutral masses ~782/787/799/1161 Da as multiple adducts at RT~1–2 min),
consistent with secreted glycolipid/polyol-lipid biosurfactants, and are strongly
species-structured (44/130 top in *R. paludigena*).
**Implications:** Provides a ranked, MS2-backed candidate list for structural annotation and
a direct genotype-phenotype handle: 270 of the paired strains also have a genome link, so
secreted-feature abundance can be associated with BGC/ortholog presence across strains.
**Tags:** metabolomics, secretion, biosurfactant, glycolipid, rhodotorula, specificity, bgc

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-02 | secreted_products v1 | Rhodotorula MS2 aligned features (paired C/SUP, 590 samples) | Rhodotorula_MS2_pheno_explore | 6,724 secreted; 130 uniquely secreted; lipid-like, species-structured | supports |

### Open Questions
- What are the actual structures? (SIRIUS/GNPS annotation of the top MS2 spectra pending.)
- Do the ~782/787/799 Da compounds map to a single glycolipid BGC present only in producer species?
- Does the effect survive a non-compositional (non-TIC) normalization?

## F-002: Effect-size alone is insufficient for secretion calls; significance + specificity gating is required
**Status:** preliminary
**Claim:** Under a sign-flip permutation null (shuffle C/SUP label within strain, 200×), the
bare median-log2FC≥1 gate has an empirical FDR ≈0.36 (observed 6,732 vs null mean 2,443).
Requiring the paired-Wilcoxon BH q<0.05 (a valid permutation test) plus cross-species Tau
specificity yields a defensible high-confidence set.
**Implications:** Any future secretion/enrichment screen on this data should not threshold on
fold-change alone; pair it with a paired significance test and report the null.
**Tags:** methodology, null-model, permutation, fdr, metabolomics

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-07-02 | secreted_products v1 | Rhodotorula MS2 (paired C/SUP) | Rhodotorula_MS2_pheno_explore | effect-size-only empirical FDR≈0.36; significance+Tau gating needed | supports |

### Open Questions
- Would a fully paired permutation that also recomputes Wilcoxon (not just sign-flip on log2FC) change the empirical FDR estimate materially?
