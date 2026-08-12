# Provenance: rhodotorula_auc_copper

## Source

**Type**: Primary

**Origin**:
- Lab: Christian Ona — copper stress-response growth-curve assay (mean_auc_rate), lab/contact not yet recorded; fill in once confirmed.

**Citation / accession**: N/A

## Acquisition details

**Date acquired**: 2026-08-11 (file placed in `input_data/`, dated in filename)

**Obtained by**: Jason Stajich

**Method**: File received/placed directly as `input_data/Rhodotorula_AUC_copper.20260811.csv.gz`. Not downloaded via script.

**Checksum**: SHA256 4e09cfa3e1b550e978742206105cab8fa848d5932b39ebd3a928603c6fcd9660 (of the gzip-compressed file as placed 2026-08-11)

## Access restrictions

**Restriction level**: none (assumed; same handling as other input_data/ phenotype files — revisit if a DUA applies)

**Details**: None.

## Known issues

- `SAMPLE_NAME = "TFCN_17-332M-1"` is duplicated (2 of 275 rows), both mapping to
  `MS2_SAMPLE_Cell/Supernatant = C_190/SUP_190` but with conflicting `mean_auc_rate`
  (22.484 vs 23.761) and `Strain ID` (190 vs 304). Root cause not yet identified —
  possibly two distinct isolates mistakenly sharing an MS2 sample ID, or a
  data-entry duplication. **Decision (2026-08-11, user directive): exclude both
  rows for C_190/SUP_190 from any downstream analysis until reconciled with the
  source.** See `.living/decisions.md`.
- 2 rows (`TFCN_122C-2`, `TFCN_138A-1`, both R. mucilaginosa) have
  `MS2_SAMPLE_Cell`/`MS2_SAMPLE_Supernatant` = literal string `"No MS2 Data"` —
  no paired metabolomics data exists for these strains; exclude from any join to
  the MS2 feature table.
- `Medium` column is essentially empty (1/275 non-null) — assay medium not
  reliably recorded in this file; confirm with source before reporting.
- `Location` missing for 70/275 (25.5%) rows.
- Unit/definition of `mean_auc_rate` (e.g., growth-curve OD/time integration
  method, copper concentration(s) used, control normalization) not yet confirmed
  against an assay protocol document — treat as a relative sensitivity/resistance
  score pending confirmation, not an absolute/calibrated unit.
- Row-level overlap between this file's MS2_SAMPLE_Cell/Supernatant values and the
  actual `C_*`/`SUP_*` columns in `Rhodotorula_MS2_aligned_features_ms2.csv.gz` has
  not been verified — check for unmatched IDs before joining in the future
  metabolite-association analysis.

## Contact

**Primary contact**: Christian Ona (copper AUC assay data)

**Backup contact**: Jason Stajich (jasonst@ucr.edu)

## Version history

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-08-11 | Initial ingestion. File was received as plain-text CSV misnamed with a `.gz` extension; corrected to a true gzip file by the source before final placement (same content, verified row-for-row identical). |
