# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-07-02] Mycelium init run in default mode on an existing repo (not --restructure)

**Context**: User asked to "initialize mycelium" in a repo that already had substantial
working structure (`input_data/`, `scripts/`, `analysis/`, `results/`) but no `.living/`.
**Decision**: Ran `init_repo.py` in **default mode**, then manually installed convention
packs and generated CLAUDE.md/INDEX.
**Alternatives considered**: (a) `--restructure` mode — rejected: it is an unimplemented
no-op (only prints an audit, moves nothing) so it would not have created `.living/`;
(b) not initializing until files were reorganized — rejected: unnecessary, default mode
is non-destructive.
**Rationale**: Read `init_repo.py` first and confirmed default mode guards every write
(`mkdir(exist_ok=True)`, `if not path.exists()`), so it only adds missing scaffolding and
cannot move or overwrite existing work. Safe on a live repo.
**Consequences**: Canonical mycelium dirs (`data/`, `algorithms/`, `reference_material/`)
now exist empty alongside the legacy dirs; documented the split in CLAUDE.md to prevent
future confusion. Existing pipeline untouched.
**Tags**: mycelium, init, repo-structure

### [2026-07-02] Installed bioinformatics domain convention pack

**Context**: init offers optional domain packs; project is genomics + metabolomics.
**Decision**: Installed `bioinformatics` (in addition to the 3 core packs). Did not
install `image-analysis` despite color-phenotype images being upstream.
**Alternatives considered**: also adding `image-analysis` — deferred: color phenotypes
arrive as pre-computed CIELab summaries, not raw images, so the pack adds little now.
**Rationale**: Upcoming work (secreted-metabolite feature tables → candidate biosynthetic
gene families) is squarely bioinformatics.
**Consequences**: `bioinformatics/analysis-conventions.md` now layers on robust-analysis
for feature-table and sequence work. Revisit image-analysis if raw plate images enter scope.
**Tags**: mycelium, conventions, bioinformatics

### [2026-07-02] Define "secretion" as a paired within-strain SUP-vs-C contrast

**Context**: The raw MS2 table has paired cell (`C_*`) and supernatant (`SUP_*`) samples per
strain; the earlier phase 0–3 pipeline pooled/collapsed them.
**Decision**: Measure secretion as a per-strain paired `log2(SUP/C)` on TIC-normalized
abundances, summarized across 292 paired strains with a paired Wilcoxon signed-rank + BH-FDR;
call features "secreted" at median log2FC≥1, q<0.05, ≥10% SUP detection.
**Alternatives considered**: (a) supernatant abundance alone (ignores that abundant
compounds may also be intracellular) — rejected; (b) two-sample SUP-vs-C across all samples
(loses pairing, exposed to strain/plate confounding) — rejected; the paired design cancels
much strain/plate nuisance.
**Rationale**: Pairing is the cleanest causal contrast for "is this exported?" and the raw
data supports it for 292/305 strains.
**Consequences**: Requires strains to have both sample types (292 kept). Compositional (TIC)
normalization caveat noted; non-compositional re-check queued as future work.
**Tags**: metabolomics, secretion, paired-design, normalization

### [2026-07-02] Uniqueness via cross-species Tau specificity (threshold 0.85)

**Context**: "Uniquely secreted" needs an operational specificity measure across 17 species.
**Decision**: Use the Tau specificity index on mean supernatant abundance per species;
uniquely-secreted = secreted AND Tau≥0.85. Report a sensitivity sweep (0.75/0.85/0.95).
**Alternatives considered**: entropy/Gini, or "present in ≤N species" hard cutoff — Tau is
continuous, standard, and comparable across features.
**Rationale**: Tau is well established for specificity and pairs naturally with a sweep.
**Consequences**: Unique-set size is Tau-sensitive (25→431 across the sweep); 130 at the mid
cut is reported as the defensible headline.
**Tags**: specificity, tau, sensitivity, metabolomics

### [2026-07-03] MS2 for SIRIUS is extracted directly from raw mzML per-feature, not via MZmine re-export

**Context**: No MGF/GNPS/SIRIUS export exists for this project (see [[raw-data-location]]
learning). Needed MS2 spectra for the 130 uniquely-secreted target features to run SIRIUS.
**Decision**: For each target, use the aligned-feature table's `adduct_rep_file` column
(MZmine's chosen representative raw file for that adduct row) to pick one mzML file, then
scan it directly with `pyteomics.mzml` for the MS2 spectrum matching the feature's m/z
(15 ppm) and RT (±0.15 min, widened to ±0.5 min on a miss), keeping the highest-TIC match.
Where `adduct_rep_file` was blank (60/130 targets), fell back to the SUP_* sample with the
highest per-sample peak area for that row.
**Alternatives considered**: (a) reopening the MZmine project to re-run its GNPS/SIRIUS
export module — rejected, no `.mzbatch`/project file was found on disk, only the flat
aligned CSV and raw mzML remain; (b) `sirius lcms-align` direct-from-mzML feature finding —
rejected for v1, since it re-detects features/alignment from scratch rather than reusing
the already-curated 130-feature target list and its QC (Tau, FDR, detection rate).
**Rationale**: `adduct_rep_file` is MZmine's own provenance field for "which raw run best
represents this row" — the most defensible single file to search first; falling back to
peak-area ranking is the natural analogue when that field is absent.
**Consequences**: A minority of targets may still fail to match (DDA duty cycle means MS2
isn't guaranteed to be triggered in the exact chosen file/window) — these are recorded as
`no_ms2_match` in `mgf_export_summary.csv` rather than silently dropped, so the true annotation
coverage is auditable.
**Tags**: metabolomics, sirius, ms2, mgf-export, mzml, secretion

### [2026-07-03] Moved SIRIUS execution from the interactive session to a dedicated sbatch job

**Context**: SIRIUS `formula`/`fingerprint`/`structure`/`canopus` on 117 spectra was
OOM-killed twice running directly in the interactive SLURM session (see [[oom-cgroup-sirius]]
learning), because that session's 16G memory cgroup was already ~9G consumed by unrelated
concurrent jobs.
**Decision**: Added `scripts/03_sirius.sbatch` (stajichlab partition/account, `--mem=32G`,
`--cpus-per-task=8`) as the canonical way to run the compute-heavy SIRIUS step, submitted with
`sbatch --chdir=<sirius_annotation dir>`. `run.sh` still documents the interactive path for
small/offline-only (`formula`+`fingerprint` without CSI:FingerID/CANOPUS) test runs.
**Alternatives considered**: (a) keep shrinking SIRIUS's own `-Xmx`/`--instance-buffer` —
rejected, headroom in the shared interactive cgroup was fundamentally too small regardless of
SIRIUS-side tuning; (b) ask the user to free up or resize the interactive session — deferred,
user chose the sbatch route when asked.
**Rationale**: A dedicated batch allocation isolates SIRIUS's memory from unrelated jobs
sharing the interactive session's cgroup, and 32G comfortably covers the ~10G heap that OOM'd
under 16G shared.
**Consequences**: Adds a SLURM dependency to the reproduce path (`run.sh` alone no longer
finishes the SIRIUS step end-to-end; `03_sirius.sbatch` must be submitted separately and waited
on). Documented in the analysis manifest.
**Tags**: hpcc, slurm, sbatch, sirius, memory, tooling

### [2026-08-11] Corrected metabolite-color correlation via joint rank-space OLS residualization instead of patching the legacy sequential-covariate loop

**Context**: `docs/FEATURE_ANALYSIS.md`'s headline metabolite-color associations
(e.g. Feature 2755, ρ=0.735) come from `scripts/phase2_correlation_analysis.py`, which
never actually took its own `'pooled' in decision['strategy']` branch (the strategy was
`stratified_with_plate`), so Species was never regressed out despite Phase 0 flagging a
significant species effect. Its `spearman_partial_corr()` helper also reassigned (rather
than accumulated) the residual on each covariate loop iteration, so with >1 covariate only
the last one would actually be controlled for, and Library Plate was regressed as a
continuous variable via `pd.factorize` rather than one-hot encoded.
**Decision**: Built a new analysis (`analysis/phenotype_metabolite_association/`) using a
single joint OLS projection (rank-transform features/phenotype, regress out one design
matrix containing one-hot Species + Library Plate + C/SUP sample_type at once, vectorized
across all 7,341 features via one matrix multiply) rather than patching the loop-based
original in place.
**Alternatives considered**: (a) patch `phase2_correlation_analysis.py` in place — rejected,
the legacy script is referenced by `docs/FEATURE_ANALYSIS.md`/`analysis/README.md` as a
specific historical run and mixing "corrected" and "original" results in the same script
would make the caveat harder to audit; (b) `statsmodels.MixedLM` with strain as a random
effect — deferred to when real phenotype replicates exist (see analysis doc's "Framework
for incorporating future phenotype replicates" section), since with only 1-2 samples/strain
a permutation-based approach that explicitly respects the C/SUP pairing is more transparent
than a mixed model with too few replicates per group to estimate variance components well.
**Rationale**: A single joint projection is both the statistically correct way to do a
partial correlation with multiple covariates (unlike the sequential-reassignment bug) and
trivially vectorizable across all features, keeping runtime under 20s for the full
22,023-test scan plus permutation/holdout checks.
**Consequences**: 0/12,269 original Tier-1 hits survive; see [[species-confound-simpsons-paradox]]
learning for the result and `analysis/phenotype_metabolite_association/PHENOTYPE_METABOLITE_ASSOCIATION.md`
for full method/results. `docs/FEATURE_ANALYSIS.md` now carries a caveats section pointing here.
**Tags**: metabolomics, statistics, partial-correlation, confounding, rhodotorula, correction

### [2026-08-11] Repeated the PLS component-count grid search inside the permutation null, rather than permutation-testing a single fixed model

**Context**: Building a multivariate/module-level test (`analysis/phenotype_metabolite_association/scripts/06_multivariate_module_test.py`)
to check for a joint many-feature signal after the univariate scan found nothing. A
naive design would pick the best PLS `n_components` once (by CV) on the real data, then
permutation-test only that fixed model.
**Decision**: The entire "grid-search over n_components, take the best CV R²" procedure
is re-run inside each of the 300 strain-block permutations, so the null distribution
already reflects the same winner's-curse-prone selection the observed statistic went
through.
**Alternatives considered**: (a) fix n_components a priori (e.g. always 5) — rejected,
an arbitrary fixed choice could either miss a real signal at the wrong component count
or (if tuned post hoc on the real data) bias the comparison in the observed statistic's
favor; (b) permutation-test each grid value separately with its own null and Bonferroni-
correct — rejected as unnecessarily conservative when the actual downstream question is
just "best achievable CV R² vs. chance," which the nested design answers directly.
**Rationale**: The univariate "top nominal" permutation tests in `03_permutation_null.py`
already demonstrated the winner's-curse problem (75/92 candidates "significant" purely
from being selected as the smallest p-values out of 7,341 tests); the same logic applies
to model-selection over components, so the selection procedure itself must be inside the
permutation loop to get a calibrated p-value.
**Consequences**: ~5 min runtime for 300 permutations (grid of 3 component counts × 4
folds × 2 fits — real + permuted — per permutation) instead of near-instant for a single
fixed model; judged worth it for a calibrated result. Observed CV R²=-0.09 vs. null p=0.56
— no joint signal detected.
**Tags**: metabolomics, statistics, pls-regression, cross-validation, permutation, model-selection, rhodotorula

### [2026-08-11] Pathway-targeted matching searches the full 16,332-row raw feature table, not the Phase-1-filtered 7,341

**Context**: Building `analysis/pathway_targeted_association/` (in response to
`phenotype_metabolite_association`'s power analysis showing the untargeted scan can't
detect true ρ<0.34), the natural first move was to match target pigment masses against
the 7,341 features already used throughout that analysis (`features_cleaned.csv.gz`).
That found only 1 hit (a low-confidence MAA).
**Decision**: Match against the full, unfiltered 16,332-row raw aligned table
(`input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz`) instead, computing a
standard total-ion-signal normalization independently for whatever the full-table
search turns up rather than reusing Phase 1's QC-survivor-only normalization
denominator.
**Alternatives considered**: (a) stick with the 7,341-feature survivor set for
consistency with the rest of the project's pipeline -- rejected once the full-table
search immediately found two additional torularhodin candidates (a diagnostic
Rhodotorula pigment) that the survivor-only search had missed entirely; (b) loosen
Phase 1's filters project-wide and redo the untargeted analysis -- rejected as
overkill and reintroducing the original multiple-testing burden this analysis exists
to avoid; a small, chemistry-motivated target list doesn't need the generic
prevalence/CV screen a genome-wide scan needs.
**Rationale**: Phase 1's blanket >=10%-of-590-samples prevalence and CV>=0.1 filters
are reasonable defaults for an untargeted scan but will systematically discard a real
pigment produced by only a subset of species/strains -- exactly the profile a
biologically meaningful, strain-specific carotenoid should have. A targeted analysis
that only searches the pre-filtered subset defeats its own purpose.
**Consequences**: Matching found 3 candidates (2 torularhodin, 1 shinorine) instead of
1; one of them (torularhodin, raw_position=12635) is the strongest, best-replicated
metabolite-color signal found anywhere in this project's phenotype-metabolite work.
Documented as a lesson in `01_match_targets.py`'s docstring for future targeted work
on this dataset.
**Tags**: metabolomics, targeted-analysis, feature-filtering, carotenoid, methodology, rhodotorula

### [2026-08-11] Set up SIRIUS structural annotation for the 3 pathway-targeted candidates, and reran the secreted_products SIRIUS batch with login fixed

**Context**: User confirmed `sirius login` is now working (verified independently via
`sirius login --show`, active academic subscription). Two SIRIUS runs were pending: the
secreted_products batch (117 spectra, twice failed at CSI:FingerID/CANOPUS on missing
login, second attempt also Slurm-cancelled at 8h walltime ~3h into CANOPUS), and the new
pathway-targeted candidates from F-006 (torularhodin x2, shinorine x1) which had never
been run through SIRIUS at all.
**Decision**: Built `analysis/pathway_targeted_association/sirius_annotation/` (scripts
00-03) mirroring `secreted_products/sirius_annotation/`'s pattern (target selection ->
MGF export via mzML scan -> sbatch SIRIUS run), with one deliberate difference: the
representative-file fallback searches both `C_*` and `SUP_*` peak-area columns (not
SUP-only, since these candidates aren't compartment-restricted the way secreted_products'
targets are). Also bumped `secreted_products/sirius_annotation/scripts/03_sirius.sbatch`'s
walltime from 8h to 20h before resubmitting, since the prior cancellation was a walltime
problem, not a resource-size problem.
**Alternatives considered**: (a) combine both target sets into one MGF/one sirius
invocation to save JVM startup overhead -- rejected, keeps the two analyses' outputs
cleanly separated under their own `sirius_annotation/outputs/sirius_project/`, and the
overhead difference is negligible next to per-compound CSI:FingerID/CANOPUS network time;
(b) reuse secreted_products' own scripts by appending the 3 new targets to its target
list -- rejected, would conflate two logically distinct analyses' provenance.
**Rationale**: Small, self-contained analysis-specific SIRIUS folders keep provenance
and output locations obvious per-analysis, consistent with this project's existing
convention (`analysis/<name>/sirius_annotation/`).
**Consequences**: Two Slurm jobs submitted (27392492: pathway-targeted, 3 spectra,
--time=02:00:00; 27392493: secreted_products rerun, 117 spectra, --time=20:00:00). Both
running as of this entry; results not yet available -- will be logged as a finding once
complete.
**Tags**: sirius, metabolomics, carotenoid, torularhodin, slurm, tooling, methodology

### [2026-08-11] Broader SIRIUS profiling combines a widened pathway search with top statistical hits, rather than loosening ppm tolerance to hit an arbitrary count

**Context**: After F-006 showed the narrow 3-candidate pathway search's mass matches
were false identifications, the user asked for SIRIUS profiling on a broader set
(~200-500 features). Widening the pathway target list alone (37 compounds across
carotenoid/sterol/melanin/MAA/flavin classes, 30 ppm instead of 15) only reached 32
matched features -- reaching 200-500 that way would require either dozens more
compound classes or a much looser ppm tolerance, the latter of which would reintroduce
exactly the false-positive-match risk F-006 just demonstrated.
**Decision**: Combined the widened pathway search (32 features,
`scripts/04_match_targets_broad.py`) with the top-80-by-nominal-p features from each of
the three already-completed corrected statistical scans (color-phenotype pooled,
copper-AUC cell, copper-AUC supernatant; 231 unique features after union) into one
263-feature SIRIUS profiling set (`sirius_annotation_broad/scripts/00_select_targets.py`).
**Alternatives considered**: (a) loosen ppm further (e.g. 50-100 ppm) to hit the count
purely via the pathway list -- rejected, dilutes the "chemistry-motivated" rationale
into noise; (b) profile only the statistical hits, dropping the pathway angle --
rejected, the pathway search is still a distinct and legitimate real question worth
keeping separate/labeled.
**Rationale**: Profiling the top statistical hits from already-validated (permutation +
holdout-checked) scans tells us what SIRIUS thinks those specific features chemically
are, independent of any pigment-pathway hypothesis -- directly useful regardless of
what class they turn out to be, and doesn't require inventing more candidate compounds
to reach a target count.
**Consequences**: `sirius_annotation_broad/outputs/target_provenance.csv` records which
source(s) each feature came from (a feature can appear in multiple tracks). MGF export
for 263 spectra takes noticeably longer than the earlier 3- and 117-spectra runs;
submitted as its own dedicated sbatch job (`sirius_broad_profiling`, 8 cores/32G/8h).
**Tags**: sirius, metabolomics, targeted-analysis, methodology, carotenoid, statistics

### [2026-08-12] SIRIUS profiling of the full pre-FDR color-phenotype candidate set (827 features), run fresh rather than reusing overlapping prior annotations

**Context**: User asked for SIRIUS on "all the features that are candidate even before
multiple testing correction," then clarified scope to color-explaining features
specifically. Filtering `phenotype_metabolite_association`'s corrected correlations to
nominal (pre-FDR) `pval_corrected<0.05` across all 3 CIELab channels (L*/a*/b*) gives
827 unique features (203/355/317 per channel, no cross-channel overlap) — far more than
the top-80-per-track subset `sirius_annotation_broad` used for the same statistical
track. 19 of these 827 already have SIRIUS results from `sirius_annotation_broad` and
`secreted_products/sirius_annotation`.
**Decision**: Created a new sibling folder,
`phenotype_metabolite_association/sirius_annotation_color_all/`, and re-ran the full
827-feature set from scratch in one self-contained SIRIUS project rather than reusing
the 19 overlapping prior results.
**Alternatives considered**: (a) tighter uncorrected threshold (e.g. p<0.01) to shrink
the set — rejected by user, who wanted the literal pre-correction set; (b) reuse the 19
already-annotated features and only export/run the 808 new ones — rejected: the ~15
minutes of extra compute to redo 19 features is far cheaper than the complexity and
provenance risk of merging two separately-run SIRIUS projects (different SIRIUS
invocations can in principle rank formula candidates slightly differently run-to-run;
one clean project avoids that ambiguity entirely).
**Rationale**: `feature_index` in `corrected_all_correlations.csv.gz` is a positional
index into `features_cleaned.csv.gz`'s columns, not a raw-table row position — this
project's existing `sirius_annotation_broad/scripts/00_select_targets.py` already
established the correct translation (column headers of `features_cleaned.csv.gz` ARE
the raw_position values); confirmed this still holds here (translated result equaled
the untranslated `feature_index` values for the color track by coincidence, verified
independently before trusting it).
**Consequences**: SLURM job `sirius_color_all` (27419982) submitted, 8 cores/32GB, 16h
budget (scaled from the broad run's 263-feature/1h28m observed rate with margin, since
per-feature runtime scaled non-linearly between the 117- and 263-feature runs — see
learnings). Results and a merge/lineage note pending job completion. Also motivated a
(not-yet-implemented) plan for profiling the entire 16,332-feature dataset via sharded
SLURM array jobs, gated on validating that SIRIUS's shared-workspace auth token
tolerates concurrent access (two prior learnings document it breaking under repeated
serial CLI invocations) — tracked in `todo/TODO_REGISTRY.md`, not started.
**Tags**: sirius, metabolomics, color-phenotype, methodology, multiple-testing, slurm
