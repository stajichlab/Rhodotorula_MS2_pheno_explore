# Summary Statistics: rhodotorula_auc_copper

<!-- Generated: 2026-08-11 -->
<!-- Script: manual (python3.12, ad hoc, csv/statistics stdlib) -->

## Overview

| Property | Value |
|----------|-------|
| Rows | 275 |
| Columns | 8 |
| File size | 29 KB (gzip-compressed) |
| Date range | N/A (not temporal) |
| Format | CSV, gzip-compressed |

## Column summaries

| Column | Type | Non-null | Unique | Min | Max | Mean | Top values |
|--------|------|----------|--------|-----|-----|------|------------|
| SAMPLE_NAME | categorical | 275 | 274 | — | — | — | 1 value duplicated (TFCN_17-332M-1, ×2) |
| SPECIES | categorical | 275 | 17 | — | — | — | R. mucilaginosa (209), R. paludigena (10), R. toruloides (9), R. dairenensis (7), R. diobovata (7) |
| mean_auc_rate | numeric | 275 | — | 0.806417566 | 29.83526913 | 22.2255 (median 23.2033, sd 4.1226) | — |
| MS2_SAMPLE_Cell | categorical | 275 | 273 | — | — | — | 2 rows = "No MS2 Data"; 2 rows share C_190 (see quality flags) |
| MS2_SAMPLE_Supernatant | categorical | 275 | 273 | — | — | — | Same pattern as MS2_SAMPLE_Cell (paired) |
| Location | categorical | 205 | — | — | — | — | Free text; not further tabulated |
| Medium | categorical | 1 | 1 | — | — | — | Effectively empty |
| Strain ID | numeric | 275 | 275 | 1 | (max not computed; all unique) | — | — |

## Missing data summary

| Column | Missing count | Missing % | Pattern / notes |
|--------|---------------|-----------|-----------------|
| Location | 70 | 25.5% | No obvious pattern checked yet |
| Medium | 274 | 99.6% | Effectively unusable as recorded |

## Quality flags

- `SAMPLE_NAME = "TFCN_17-332M-1"` occurs twice (rows 164–165 of the source CSV),
  both mapped to `MS2_SAMPLE_Cell/Supernatant = C_190/SUP_190` but with different
  `mean_auc_rate` (22.484 vs 23.761) and `Strain ID` (190 vs 304) — **excluded from
  downstream use per 2026-08-11 decision**, see `.living/decisions.md`.
- 2 rows (`TFCN_122C-2`, `TFCN_138A-1`) have `MS2_SAMPLE_Cell`/`MS2_SAMPLE_Supernatant`
  = `"No MS2 Data"` — not linkable to the metabolome dataset; exclude from any join.
- `Medium` column carries essentially no information (1/275 populated) — do not use.
- Exact ID overlap between this file's 273 usable `MS2_SAMPLE_Cell` values and the
  `C_*` sample columns in `Rhodotorula_MS2_aligned_features_ms2.csv.gz` was **not**
  verified at ingestion time — verify before joining in the future association analysis.

## Notes

- After excluding the 2 "No MS2 Data" rows and the 2 C_190/SUP_190 duplicate rows,
  271 rows remain usable for a metabolite-composition ~ copper-AUC association analysis.
- Species distribution is heavily imbalanced (R. mucilaginosa = 209/275, 76%) —
  any cross-species analysis will need to account for this, consistent with the
  species-confound handling already established for the color-phenotype work
  (see `.living/decisions.md`, entries on species stratification).
