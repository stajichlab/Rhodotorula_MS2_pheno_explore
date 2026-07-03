# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-07-02] Mycelium scripts need Python 3.11+; cluster default is 3.9

**Category**: tooling / environment
**What happened**: `init_repo.py` (and other mycelium scripts) fail with `ImportError: cannot import name 'UTC' from 'datetime'` under the cluster default `python3` (miniconda 3.9). `/usr/bin/python3.12` runs them fine.
**Why it matters**: Every future mycelium script invocation in this repo must use a ≥3.11 interpreter or it errors before doing anything.
**Resolution**: Use `/usr/bin/python3.12` for all mycelium scripts. Documented in `CLAUDE.md` Quick Orientation tooling note. Neither the 3.9 nor the 3.12 interpreter has `pyyaml`; mycelium scripts degrade gracefully without it (text-fallback parsing), but any YAML you need to validate can't be parsed by either interpreter as-is.
**Tags**: mycelium, python-version, environment, tooling
**mitigation_type**: convention
**structural_mitigation_candidate**: a repo-local wrapper (e.g. `bin/myc` that pins `/usr/bin/python3.12`) would make the version requirement structural rather than remembered.

### [2026-07-02] init_repo.py auto-install of core packs silently no-ops; ACTIVE_CONVENTIONS.yaml stub is malformed

**Category**: tooling / bug
**What happened**: (1) `init_repo.py` computes the network dir as `skills/network/conventions` (one level too shallow — should be `<mycelium>/network/conventions`), so core convention packs are NOT auto-installed; it prints a warning but proceeds. (2) After installing packs manually with `install_convention.py`, the resulting `ACTIVE_CONVENTIONS.yaml` keeps the init template's `active_conventions: []` stub with the pack entries appended as top-level list items — malformed YAML.
**Why it matters**: A fresh `init` leaves the repo with zero conventions unless you notice the warning; the malformed registry can break any tool that YAML-parses it.
**Resolution**: Installed robust-analysis, report-generator, idea-generator, bioinformatics manually via `install_convention.py --network-dir <mycelium>/network/conventions`; rewrote `ACTIVE_CONVENTIONS.yaml` as a proper list under the `active_conventions:` key.
**Tags**: mycelium, init, bug, conventions, yaml
**mitigation_type**: ambient-awareness
**structural_mitigation_candidate**: a post-init assertion that `ACTIVE_CONVENTIONS.yaml` parses and lists ≥3 core packs would catch both failure modes.

### [2026-07-02] scilintr install on this cluster: use `uv tool install`, not pip

**Category**: tooling / environment
**What happened**: `scilintr` requires Python ≥3.11. The cluster's pip is bound to conda
Python 3.9 (rejects scilintr as "requires-python >=3.11"), and `/usr/bin/python3.12` has **no
pip module** at all. `uv tool install scilintr` (uv is at `~/.local/bin/uv`) resolved and
installed it cleanly to `~/.local/bin/scilintr`; a `python3.12 -m venv` + ensurepip also works.
**Why it matters**: The analyze convention mandates running scilintr after code edits; without
a working path it gets skipped.
**Resolution**: `uv tool install scilintr` → run `~/.local/bin/scilintr <file>`. The analysis
script now lints clean (rc=0) using ANALYSIS_OK waivers for intentional asserts/filters.
**Tags**: scilintr, uv, python-version, environment, tooling, analyze-convention
**mitigation_type**: convention
**structural_mitigation_candidate**: add the `uv tool install scilintr` step and the
`/usr/bin/python3.12`-for-mycelium note to ENVIRONMENTS_INSTALLATIONS.md so the lint path is
discoverable without re-deriving it.

### [2026-07-02] The paired C/SUP structure is the key latent axis in this MS2 dataset

**Category**: data / domain
**What happened**: The raw MS2 table encodes each strain twice — `C_*` (cell pellet) and
`SUP_*` (supernatant), ~295 each, pairable by `ATTRIBUTE_ID_1` (not `Standardized Strain`,
which is 50% null). `ATTRIBUTE_species` is the usable species label (17 species; `Species`
is 54% null). Genome linkage lives in `db_sra_run_list`/`db_biosample_list`/`db_bioproject_list`.
**Why it matters**: This pairing (collapsed by the phase 0–3 pipeline) is what makes secretion,
and genotype↔metabolite association, directly measurable; the right ID columns are non-obvious.
**Resolution**: Use `ATTRIBUTE_ID_1` for pairing, `ATTRIBUTE_species` for species, the db_*
columns for genome availability. Parse the metadata TSV with pandas (embedded newlines break awk).
**Raw data location** (parent dir `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/`):
cell mzML `mzML/C_*.mzML` (299, ~12 GB); supernatant mzML `ExFab_Supernatant/SUP_*.mzML`
(304, ~13 GB); MZmine workspace `feature_extractMS2/`. No MGF/GNPS/SIRIUS export exists yet —
MS2 must be exported (from MZmine or mzML) before SIRIUS/GNPS annotation. Genomes are NOT in
this tree (linked via SRA/BioSample accessions in the metadata; user holds them separately).
**Tags**: metabolomics, metadata, pairing, secretion, rhodotorula, data-structure, mzml, raw-data
**mitigation_type**: ambient-awareness
**structural_mitigation_candidate**: document these canonical key columns in DATA_DICTIONARY.md.

### [2026-07-03] SIRIUS is installed on the cluster and usable both as a module and a conda env

**Category**: tooling / environment
**What happened**: `module load sirius` (5.8.1) works standalone (auto-loads `java` +
`cplex-studio`). A conda env at
`/bigdata/stajichlab/jstajich/projects/Metabolomics_Workshop/sirius-ms` has the same version.
CSI:FingerID structure search and CANOPUS compound-class prediction call Boecker-lab web
services and require `sirius login` first (currently not logged in) plus outbound internet
from the node running it; `formula`/`fingerprint` (molecular formula + fragmentation tree)
work fully offline.
**Why it matters**: The `feature_index` column in `secretion_scores_all_features.csv.gz` /
`uniquely_secreted_features.csv` (analysis/secreted_products) is the **0-based positional row
index** in the aligned feature table, not the MZmine `row ID` column — confirmed by joining
feature_index 195/3043 to row IDs 205/3663 by position. Any script joining these two files
must join on position, not on an ID column.
**Resolution**: `analysis/secreted_products/sirius_annotation/scripts/00_select_targets.py`
does this positional join once and carries `row_id`/`adduct_rep_file` forward.
**Raw MS2 export gap**: no MGF/GNPS/SIRIUS export exists anywhere in this project tree; MS2
must be pulled directly from `mzML/C_*.mzML` / `ExFab_Supernatant/SUP_*.mzML` per feature
(see [[ms2-extraction-from-mzml]] decision). `pyteomics` 5.x needs Python ≥3.10 (`X | Y`
union-type syntax) — install/run it under `/usr/bin/python3.12`, not the cluster default
py3.9 conda, and needs `lxml` + `psims` as companion packages.
**Tags**: sirius, tooling, environment, feature-index, mzml, ms2, metabolomics
**mitigation_type**: ambient-awareness
