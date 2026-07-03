## What was worked on
- Set up `analysis/secreted_products/sirius_annotation/`: SIRIUS+CSI:FingerID MS2
  structural-annotation pipeline for the 130 uniquely-secreted features from
  `secreted_products`. Confirmed SIRIUS 5.8.1 is installed on the cluster
  (`module load sirius`, plus a conda env). Wrote scripts to (1) join targets to
  the aligned feature table's `adduct_rep_file` (with a max-peak-area fallback
  for 60/130 blank cases), (2) extract the matching MS2 spectrum directly from
  raw mzML via `pyteomics` (no MGF/GNPS export existed anywhere in the tree),
  (3) run SIRIUS `formula`/`fingerprint`/`structure`/`canopus`.

## Current state
- Branch: `main`, prior commit `fb6fb87` pushed; new sirius_annotation files not
  yet committed.
- MGF export (`scripts/01_export_mgf.py`, run under `/usr/bin/python3.12`) was
  running in the background as of end of session — check
  `outputs/mgf_export_summary.csv` for completion and match-rate before running
  `scripts/02_run_sirius.sh`.
- SIRIUS `structure`/`canopus` steps need `sirius login` (not yet done) plus
  outbound internet from whichever node runs it; `formula`/`fingerprint` work
  offline.

## Next steps
1. Confirm MGF export completed and check `no_ms2_match` rate in the summary.
2. Run `scripts/02_run_sirius.sh` (or full `run.sh`).
3. `sirius login` if compound-class/structure annotation is wanted.
4. Feed firm compound classes back into `SECRETED_PRODUCTS.md` next-steps and
   the candidate gene-family mapping.
