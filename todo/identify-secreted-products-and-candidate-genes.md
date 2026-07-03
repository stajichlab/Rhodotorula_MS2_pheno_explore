# Identify uniquely & highly secreted products → candidate biosynthetic gene types

| Field | Value |
|-------|-------|
| **Date** | 2026-07-02 |
| **Author** | Jason Stajich |
| **Priority** | high |
| **Status** | open |
| **Category** | analysis |
| **Related analyses** | phase0–phase3 pipeline (`analysis/`, `results/`) |
| **Related data** | `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` (paired `C_*` cell / `SUP_*` supernatant columns); `input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz` |

## Description

Identify MS2 metabolite features that are **uniquely and highly secreted** — i.e.
enriched in supernatant (`SUP_*`) relative to paired cell (`C_*`) samples, and present
in a restricted set of strains/species (uniqueness) at high abundance — then map those
metabolite classes to **candidate biosynthetic gene families** to screen for in the
*Rhodotorula* genome datasets.

## Motivation

The existing phase 0–3 pipeline correlates ~7,300 filtered features to color phenotype
but collapses the cell-vs-supernatant distinction. The raw table actually contains
**paired `C_*` (cell) and `SUP_*` (supernatant) samples — ~295 of each** — so secretion
is directly measurable as a within-sample SUP−C contrast. This is the missing axis needed
to answer "what is secreted" and to prioritize which biosynthetic gene classes explain
strain/species differences.

## Proposed Approach

1. **Confirm pairing**: match `C_<id>` ↔ `SUP_<id>` columns by strain id; report how many
   form complete pairs and which are singletons.
2. **Secretion score**: per feature, compute a paired SUP-vs-C enrichment (e.g. log2
   ratio / paired Wilcoxon across strains) on appropriately normalized areas. Flag
   "secreted" = significantly SUP-enriched.
3. **Uniqueness**: for secreted features, measure specificity across strains/species
   (e.g. Gini/entropy of abundance, or "present in ≤N species"). "Uniquely secreted" =
   high secretion score AND high specificity.
4. **Highly secreted**: rank by absolute supernatant abundance among the unique+secreted set.
5. **Annotate chemistry**: pull m/z, RT, adduct, parent_mass, `has_ms2` for the top hits;
   attempt class assignment (carotenoids/terpenoids, polyketides, NRPS peptides,
   glycolipids/biosurfactants, siderophores, etc.).
6. **Map to gene families**: for each metabolite class, list candidate biosynthetic gene
   types to search genomes for — e.g. carotenoid pathway (`crtYB`, `crtI`, `crtE`, `CAR`
   cluster), terpene synthases, PKS (type I/III), NRPS / NRPS-like, hybrid PKS-NRPS,
   fatty-acid/glycolipid biosynthesis, transporters/ABC exporters (secretion), CYP450
   tailoring enzymes. Cross-reference with antiSMASH BGC predictions on the genomes.

## Acceptance Criteria

- [ ] Ranked table of uniquely + highly secreted features (feature id, SUP/C enrichment, specificity, abundance, m/z, RT, MS2 availability).
- [ ] Chemical-class hypotheses for the top features.
- [ ] Candidate biosynthetic gene-family shortlist per class, framed as a genome-screening plan (incl. antiSMASH / orthology search).

## Notes

- Watch batch effects: `Library Plate` was a significant covariate in phase 0; secretion
  contrasts should stay within paired strains to cancel much of it.
- Normalization choice matters for ratios — revisit total-area normalization vs. a
  secretion-appropriate scheme before computing SUP/C.
- This item is the original driver for initializing mycelium in this repo (2026-07-02).
