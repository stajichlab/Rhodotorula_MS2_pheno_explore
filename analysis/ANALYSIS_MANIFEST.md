# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

## secreted_products

- **Path:** `analysis/secreted_products/` (`SECRETED_PRODUCTS.md`, `run.sh`, `scripts/01_secretion_analysis.py`)
- **Date:** 2026-07-02 · **Status:** complete (v1)
- **Question:** Which MS2 metabolites are uniquely & highly secreted (supernatant vs cell), and what biosynthetic gene families should we screen for in the genomes?
- **Inputs:** `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` (paired `C_*`/`SUP_*` columns), `input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz`
- **Method:** paired SUP-vs-C log2FC + Wilcoxon signed-rank + BH-FDR; cross-species Tau specificity; sign-flip permutation null; threshold sensitivity sweep.
- **Headline:** 292 paired strains (270 with genome link); 6,724 secreted, **130 uniquely secreted** features (all with MS2); top hits are lipid-like/glycolipid, strongly species-structured.
- **Key outputs:** `outputs/uniquely_secreted_features.csv`, `outputs/candidate_gene_families.md`, `outputs/figures/*.png`, `outputs/numbers.json`
- **Lineage:** uses the same raw MS2 table as the phase 0–3 color-phenotype pipeline, but exploits the C/SUP pairing those phases collapsed.
