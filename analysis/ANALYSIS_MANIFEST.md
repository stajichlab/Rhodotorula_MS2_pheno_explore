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

## secreted_products/sirius_annotation

- **Path:** `analysis/secreted_products/sirius_annotation/` (`run.sh`, `scripts/00_select_targets.py`, `scripts/01_export_mgf.py`, `scripts/02_run_sirius.sh`)
- **Date:** 2026-07-03 · **Status:** in-progress (v1, MGF export complete, SIRIUS formula/fingerprint/structure/canopus running)
- **Question:** What are the actual molecular formulas/structures/compound classes behind the 130 uniquely-secreted MS2 features (`secreted_products` open next-step #1)?
- **Inputs:** `analysis/secreted_products/outputs/uniquely_secreted_features.csv` (130 targets); `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` (for `adduct_rep_file`/per-sample peak areas); raw `mzML/C_*.mzML` and `ExFab_Supernatant/SUP_*.mzML` (parent dir, not in this repo).
- **Method:** positional-index join of targets to the aligned table → resolve one representative raw mzML per feature (`adduct_rep_file`, or max-peak-area fallback for 60/130 blank cases) → scan that mzML with `pyteomics` for the MS2 spectrum matching m/z (15 ppm) + RT (±0.15→0.5 min), keep highest-TIC match → MGF → SIRIUS `formula`/`fingerprint`/`structure`/`canopus`.
- **Headline:** 117/130 targets (90%) matched an MS2 spectrum and were exported (13 `no_ms2_match`, all in the ±0.15 window — widening to ±0.5 min recovered only 1 more); DDA duty cycle is the expected cause of the gap, not a script bug.
- **Key outputs:** `outputs/sirius_targets.csv`, `outputs/sirius_targets.mgf` (117 spectra), `outputs/mgf_export_summary.csv` (per-feature match audit), `outputs/sirius_project/` (SIRIUS results).
- **Caveats:** CSI:FingerID structure search and CANOPUS need `sirius login` (web service, needs internet) — not yet done as of v1; `formula`/`fingerprint` work offline, so molecular-formula/fragmentation-tree calls are unaffected by the login gap.
- **Lineage:** child of `secreted_products` — consumes its 130 uniquely-secreted feature list directly.
