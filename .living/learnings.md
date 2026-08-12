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

### [2026-07-03] Running SIRIUS interactively OOM-killed twice; the cause was shared cgroup memory, not SIRIUS's own heap size

**Category**: tooling / environment
**What happened**: Running `sirius formula fingerprint structure canopus` on 117 spectra
directly in the interactive SLURM session (job 26025365, `--mem=16G`) was OOM-killed twice
(confirmed via `dmesg`: `oom-kill:constraint=CONSTRAINT_MEMCG ... task=java`), even after
capping the JVM heap with `JAVA_OPTS="-Xmx10G"` and lowering `--cores 2 --instance-buffer 1`.
`free -h` misleadingly reports the whole physical host's RAM (503G), not the cgroup limit —
SIRIUS's own launcher sets `-XX:MaxRAMPercentage=65` against that host-visible value, so on a
small SLURM allocation the JVM heap sizing is wrong from the start. Worse: the interactive
job's 16G cgroup was already ~9G consumed by *unrelated* concurrently-running jobs in the same
session (`conda-build` ~6.2G, `funannotate` java, `raxml-ng`), leaving only ~6-7G of real
headroom no matter how SIRIUS's own flags were tuned.
**Why it matters**: Any memory-hungry tool run inside a shared interactive HPCC session can be
killed by *other unrelated jobs'* memory use in the same cgroup — this is invisible from `free
-h` / `ps` alone; you have to check `sacct`/`scontrol show job`/`dmesg` to see the real ceiling.
**Resolution**: Submit as a dedicated `sbatch` job with its own `--mem` allocation instead of
running inside the shared interactive session (see `analysis/secreted_products/
sirius_annotation/scripts/03_sirius.sbatch`, `--mem=32G`). Also see [[hpcc-sbatch-script-paths]]
memory for the `$0`/spool-dir path bug hit while building that same sbatch script.
**Tags**: hpcc, slurm, oom, memory, cgroup, sirius, sbatch, tooling
**mitigation_type**: ambient-awareness

### [2026-08-11] The pooled metabolite-color correlations in docs/FEATURE_ANALYSIS.md were a species confound (Simpson's paradox), not a real signal

**Category**: statistics / metabolomics
**What happened**: Requested to build confidence in the strong metabolite-color
correlations reported in `docs/FEATURE_ANALYSIS.md` (e.g. Feature 2755, ρ=0.735,
q=2.4e-93 for brightness). Tracing the actual code (not just the docs) found Species was
never regressed out despite Phase 0 detecting a significant species effect (see the
paired [[joint-covariate-correction]] decision for the code-level cause). Rebuilding the
correlation with Species + Library Plate + C/SUP sample_type jointly regressed out
(`analysis/phenotype_metabolite_association/`) collapsed Feature 2755 from ρ=0.735 to
ρ=0.024 (permutation p=0.69); **0 of the original 12,269 Tier-1 hits survived** the
corrected model, a permutation null, or holdout replication (2/75 top corrected
candidates replicated out-of-sample, chance level). Restricting to *R. mucilaginosa*
alone (n=415, 210 strains — the only species with enough strains to test this, and with
real a*/b* phenotype spread: CV 13% and 35% respectively) did not recover a stronger
within-species signal either (max |ρ|=0.196, 0 FDR hits).
**Why it matters**: A pooled correlation across multiple species/strains with real
between-species phenotype and metabolite differences will produce large, highly
"significant" correlations even with zero true within-species relationship — the
metabolomics analogue of population stratification in GWAS (this project's own
`docs/GWAS_EXPERT_EVALUATION.md` flags the same failure mode on the genotype side).
Species (or any strong grouping variable) must be verified as actually controlled for by
reading the code path taken, not just the intended design (`phase0_decision.json` said
`stratified_with_plate`; the code silently didn't act on it).
**Resolution**: Use `analysis/phenotype_metabolite_association/` outputs, not
`docs/FEATURE_ANALYSIS.md`'s ρ/q-values, for any future claim about specific metabolite
features driving color phenotype. That doc now carries a caveats section pointing here.
**Tags**: metabolomics, statistics, confounding, simpsons-paradox, species, rhodotorula, correlation, validation
**mitigation_type**: process-fix

### [2026-08-11] `Species` column in phase1_phenotype_data.csv.gz is NaN for ~55% of rows (almost all SUP_* samples)

**Category**: data-quality
**What happened**: While building the species-corrected re-analysis, found
`phase1_phenotype_data.csv.gz`'s `Species` column (used throughout the legacy Phase 0-3
pipeline, including Phase 0's species-confound F-test) is NaN for 321/590 samples —
nearly every `SUP_*` (supernatant) row — because it was populated only from cell-pellet
metadata upstream. `ATTRIBUTE_species` in
`input_data/MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz` is populated
for both C_* and SUP_* (only 30/590 NaN) and was used instead. Also found one strain ID
collision: `17-332Y-1` maps to two distinct C/SUP filename pairs (C_165/SUP_165,
C_269/SUP_269), with SUP_165 and SUP_269 carrying byte-identical phenotype values — an
upstream metadata bug, not a real replicate; those 4 rows are dropped.
**Why it matters**: This means Phase 0's own species-effect significance test
(F=32.65, p=1.1e-16) was itself computed on a `Species` column that was blank for most
supernatant samples — worth keeping in mind if that number is ever cited going forward.
**Resolution**: `analysis/phenotype_metabolite_association/scripts/01_prepare_data.py`
uses `ATTRIBUTE_species` and drops the collision strain + any remaining missing-covariate
rows (40/590 dropped total, 550 retained), with counts logged rather than silently coerced.
**Tags**: metabolomics, data-quality, species, metadata, rhodotorula, join
**mitigation_type**: process-fix

### [2026-08-11] Strain-ID-to-MS2-sample collisions recur across independent phenotype files — check for them on every new phenotype ingestion

**Category**: data-quality
**What happened**: The newly ingested `rhodotorula_auc_copper` dataset
(`input_data/Rhodotorula_AUC_copper.20260811.csv.gz`) has the same shape of bug already
seen in the color-phenotype file: `SAMPLE_NAME = "TFCN_17-332M-1"` appears twice, both
mapped to the same `MS2_SAMPLE_Cell`/`MS2_SAMPLE_Supernatant` pair (`C_190`/`SUP_190`) but
with two different `mean_auc_rate` values and two different `Strain ID` values (190 vs
304) — directly analogous to the `17-332Y-1` → `{C_165/SUP_165, C_269/SUP_269}` collision
found earlier in `phase1_phenotype_data.csv.gz` (see entry above). Two independent
phenotype files from this project both had a strain incorrectly sharing/duplicating an
MS2 sample ID.
**Why it matters**: This is not a one-off fluke in a single file — it looks like a
recurring upstream metadata issue (possibly in how strain IDs get assigned to MS2 sample
runs). Any new phenotype file joined to the MS2 feature table via `C_*`/`SUP_*` IDs should
be checked for duplicate/conflicting sample-ID mappings before use, not assumed clean.
**Resolution**: For `rhodotorula_auc_copper`, both `C_190`/`SUP_190` rows are excluded
pending reconciliation (user decision, see `.living/decisions.md`), rather than guessing
which value is correct. Recommend running the same duplicate-ID check the color-phenotype
pipeline already does (`01_prepare_data.py`) on every future phenotype ingestion that
joins via `MS2_SAMPLE_Cell`/`MS2_SAMPLE_Supernatant`.
**Tags**: metabolomics, data-quality, metadata, rhodotorula, join, duplicate-records, copper
**mitigation_type**: process-fix

### [2026-08-11] SIRIUS 5.8.1's CSI:FingerID/CANOPUS calls fail with a 404 loop, not a login error -- the installed CLI is a major version behind the web API

**Category**: tooling / environment
**What happened**: After the user fixed the `sirius login` issue (confirmed working via
`sirius login --show`, active academic subscription), both SIRIUS jobs (pathway-targeted
candidates and the secreted_products rerun) hung in an infinite exponential-backoff retry
loop: `WARNING: Request to Server failed! ... Bad HTTP Response Code: 404 ... GET:
https://academic.csi.bright-giant.com/v2.6/api/fingerid/data?predictor=1`. This looks
superficially like another login problem but isn't -- it's a **client/server API version
mismatch**: the only `module load sirius` available (5.8.1) is calling an old
`/v2.6/api/...` path that the current Bright Giant server no longer serves. A web search
confirmed the current upstream release is SIRIUS 6.3.12 (5.8.1 -> 6.3.12 is a major
version jump). Both jobs were killed (`scancel`) rather than left to burn their 2h/20h
walltime allocations on a loop that would never succeed.
**Why it matters**: A 404 in SIRIUS's web-service calls after login is confirmed working
should prompt checking the installed CLI version against the current upstream release,
not more login debugging -- the error message doesn't distinguish the two failure modes.
**Resolution**: User installed `sirius/6.3.12` as a new module. Hit two follow-on snags,
both resolved: (1) the module file initially pointed at
`/opt/linux/rocky/8.x/x86_64/pkgs/sirius/6.3.12/bin`, but the package was actually
unpacked under a sibling `6.0/` directory -- a module-file path bug, not a broken
install (confirmed the `6.0/bin/sirius --version` binary correctly self-reports as
"SIRIUS 6.3.12"); user fixed the module file's path. (2) SIRIUS keeps a
**separate workspace per major version** (`~/.sirius-5.8` vs. `~/.sirius-6.3`), so the
already-working 5.8.1 login does not carry over to 6.3.12 -- a fresh `sirius login` is
required under the new module before rerunning either job.
**Tags**: sirius, tooling, environment, hpcc, module, api-version, csi-fingerid
**mitigation_type**: ambient-awareness

### [2026-08-11] SIRIUS 6.x restructured its CLI subcommand chain: `structures` now depends on `canopus`, not the other way around

**Category**: tooling / environment
**What happened**: After switching to `sirius/6.3.12` (see the prior 5.8.1-vs-6.3.12
learning), the existing pipeline command (`formula ... fingerprint structure canopus
write-summaries`, unchanged since the SIRIUS 5.x era) failed immediately with
`picocli.CommandLine$UnmatchedArgumentException: Unmatched argument at index N:
'structure'`. Checking each subcommand's own `--help` (`sirius formula --help`,
`sirius fingerprint --help`, `sirius canopus --help`) revealed the valid chain
changed: `structures` (also accepts `structure-db-search`/`structure` as aliases,
but only in the right position) is listed as a child `Command:` of `canopus`, not of
`fingerprint` -- the correct order in 6.x is `formula -> fingerprint -> canopus ->
structures -> write-summaries`, not the 5.x order (`formula -> fingerprint ->
structure -> canopus -> write-summaries`).
**Why it matters**: A parse-time `UnmatchedArgumentException` on a subcommand name
that clearly exists (confirmed via `--help`) usually means a *chaining order*
problem, not a missing/renamed subcommand -- worth checking each subcommand's own
`Commands:` section in its `--help` output before assuming the flag/subcommand was
removed.
**Caution for next time**: diagnosing this involved several rapid, live `sirius`
CLI invocations against the shared `~/.sirius-6.3` workspace in quick succession;
this appears to have raced on the auth/refresh-token file and broke the just-confirmed
login (`sirius login --show` went from showing an active subscription to "Not logged
in" with no login command run in between). Prefer `--help` (no auth touched) over
live dry-runs when debugging CLI syntax against a shared, authenticated workspace.
**Resolution**: Fixed the subcommand order in all 4 copies of the pipeline command
(`analysis/{secreted_products,pathway_targeted_association}/sirius_annotation/scripts/{02_run_sirius.sh,03_sirius.sbatch}`
-- the interactive and sbatch variants in both analyses). Login needs to be redone
before either job can be resubmitted.
**Tags**: sirius, tooling, environment, cli, api-version, csi-fingerid, canopus
**mitigation_type**: ambient-awareness
