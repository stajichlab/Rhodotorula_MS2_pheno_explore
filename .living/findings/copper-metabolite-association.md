---
topic: copper-metabolite-association
description: Does metabolite composition (cell or supernatant) predict copper stress-response AUC, once species/plate confounds are controlled for?
created: 2026-08-11
last_updated: 2026-08-11
status: active
---

# Copper-metabolite association

## F-001: No univariate feature (cell or supernatant) survives FDR correction for copper AUC; result is underpowered, not disproven
**Status:** preliminary (multivariate PLS/RF module tests still running as of this entry)
**Claim:** Across 264 strains/17 species in each of two independent tracks (cell `C_*`,
supernatant `SUP_*`, run in parallel rather than pooled — see `.living/decisions.md`),
joint rank-space partial-Spearman correlation of ~7,340 QC'd metabolite features against
copper AUC (`mean_auc_rate`), correcting for species (+Library Plate in the cell track;
Plate is constant/unusable in the supernatant track — see F-003), finds **0/7,341 (cell)
and 0/7,333 (supernatant) features at FDR<0.05** (max corrected |ρ| = 0.224 cell, 0.214
supernatant). A simulation-based power analysis (same method as
`analysis/phenotype_metabolite_association`'s F-005) puts the minimum detectable ρ at 80%
power at **0.340 (cell) and 0.345 (supernatant)** — both well above the observed max |ρ|
of ~0.22, meaning this design can only confidently rule out associations ≥~0.34, not the
0.15–0.30 range where the actual observed effect sizes sit.
**Implications:** Directly mirrors the color-phenotype analysis's own calibrated
conclusion (F-005 there): absence of an FDR hit here is **underpowered, not disproven**
for a true modest effect. Species-stratified strain-level holdout replication of the
top-25 nominal candidates in each track also failed (0/25 replicated, same-sign+p<0.05,
in both tracks; global train/test rank calibration was -0.226 [cell] and +0.031
[supernatant], i.e. at or below chance) — consistent with the FDR/power result, not an
independent contradiction of it.
**Open questions:** (1) Do the pending PLS/Random Forest joint multivariate tests find a
signal that no single feature carries individually? (2) Would restricting to
*R. mucilaginosa* alone (the dominant species, ~200/264 strains, mirroring the
color-phenotype analysis's within-species re-run) change anything, or would raising
statistical power require a fundamentally different design (more phenotype replicates,
a smaller chemistry-motivated target list, mirroring the color analysis's successful
torularhodin pathway-targeted approach)?
**Tags:** copper, metabolomics, statistics, species-confound, power-analysis, negative-result, rhodotorula

## F-002: Copper AUC itself is strongly species-structured, confirming the confound-correction step is necessary
**Status:** preliminary
**Claim:** Kruskal-Wallis test of `mean_auc_rate` across species (10 groups with n≥3)
gives H=61.88, p=5.8e-10 in both tracks (species composition is identical since both
tracks derive from the same strain set). Species medians range from ~7.6 (R.
kratochvilovae, n=3) to ~25.2 (R. pacifica, n=2), with the dominant R. mucilaginosa
(n=201) at 23.6, well above several smaller species.
**Implications:** The species-partial correction applied throughout this analysis (and
the color-phenotype analysis before it) is doing real, necessary work, not a formality —
without it, any metabolite that happens to track species identity would show a spurious
pooled correlation with copper AUC via the same Simpson's-paradox mechanism already
confirmed for color phenotype.
**Tags:** copper, metabolomics, species-confound, kruskal-wallis, rhodotorula

## F-003: Library Plate is confounded with copper AUC in the cell track but constant (unusable as a covariate) in the supernatant track
**Status:** preliminary (data-quality note, not itself a copper-biology claim)
**Claim:** Kruskal-Wallis of `mean_auc_rate` by `Library Plate`: H=21.19, p=9.6e-5 in the
cell track (4 plates, n=50-76 each). In the supernatant track, **all 264 rows have
Library Plate = 1.0** — verified against the full, unfiltered `phase1_phenotype_data.csv.gz`
(295/295 SUP_* rows, not just this analysis's subset), so this is a genuine upstream data
property (plate metadata apparently only tracked for the cell-pellet batch), not a join
bug in this analysis.
**Implications:** Library Plate is included as a covariate for the cell track only; a
constant covariate cannot be regressed out (zero variance) and would be silently dropped
by the design-matrix rank check if included anyway, so this is handled explicitly rather
than left to fail quietly.
**Tags:** copper, metabolomics, data-quality, library-plate, rhodotorula, methodology

## F-004: PLS multivariate module test finds a weak but permutation-significant joint signal in both tracks, despite negative absolute CV R²
**Status:** preliminary
**Claim:** `07_multivariate_module_test.py` (PLS regression, GroupKFold CV, component
count selected inside a 300-permutation null): observed best CV R² is **negative** in
both tracks (cell: -0.0496 at 10 components; supernatant: -0.0673 at 10 components) —
the model does not usefully predict copper AUC out-of-fold. However, permutation testing
shows this is **significantly less negative than chance** (cell: null mean R²=-0.217,
permutation p=0.017; supernatant: null mean R²=-0.226, permutation p=0.020) — a PLS
model fit to the real feature/phenotype relationship generalizes less poorly than one
fit to label-shuffled data, in both tracks independently.
**Implications:** This is a genuinely different result from the color-phenotype
analysis's own PLS test (CV R²=-0.09, permutation p=0.56 — clearly not significant
there). It's consistent with a real but weak joint multivariate signal spread across
many features, too diffuse for any single feature to clear univariate FDR (F-001) and
far too weak to be practically predictive (R² still negative), but statistically
distinguishable from a true null. Should not be oversold: "permutation-significant
negative R²" is not the same as "the model works" — it means the *degree of badness* is
informative, not that the model is good. **Random Forest corroboration (added below,
F-005) is only partial** — cell track corroborates, supernatant does not — so this
finding should be read together with F-005, not in isolation.
**Tags:** copper, metabolomics, pls, multivariate, permutation, weak-signal, rhodotorula

## F-005: Random Forest corroborates the weak joint signal in the cell track only; supernatant does not replicate across model classes
**Status:** preliminary
**Claim:** `08_random_forest_module_test.py` (same nested grid-search-in-permutation
design as F-004, N_PERM=200): cell track observed CV R²=0.0005 (max_depth=None,
essentially flat) vs. null mean -0.062 [-0.128, 0.0003], **permutation p=0.030**
(significant). Supernatant track observed CV R²=-0.0288 (max_depth=6) vs. null mean
-0.052 [-0.116, 0.003], **permutation p=0.199** (not significant — observed value
falls within the null's own 95% range).
**Implications:** Combined with F-004 (PLS significant in both tracks), the two model
classes only **agree** on the cell track: PLS p=0.017 and RF p=0.030 both clear p<0.05
independently, which is the kind of cross-model-class corroboration this project treats
as meaningful evidence (mirroring the color-phenotype analysis's practice of requiring
both PLS and RF nulls to agree before calling a multivariate result settled). The
supernatant track is significant under PLS only (p=0.020) — a single-model-class result
that could be a PLS-specific artifact (sensitivity to a particular linear combination
RF's splits don't capture) or a false positive from running 2 tracks × 2 model classes
= 4 tests with no correction across them. **Net read: tentative evidence for a weak,
diffuse joint metabolite signal in the cell (intracellular) fraction specifically; the
supernatant multivariate result should be treated with more skepticism than the cell
result.**
**Open questions:** What features drive the cell-track PLS/RF signal (loadings/importances
not yet cross-referenced)? Does a dedicated holdout for the multivariate signal itself
(not just univariate hits) replicate it? Would narrowing to a smaller, chemistry-motivated
feature set (mirroring the successful color-phenotype torularhodin lead) sharpen this
into something identifiable?
**Tags:** copper, metabolomics, random-forest, multivariate, permutation, weak-signal, model-corroboration, rhodotorula

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-11 | copper_metabolite_association v1 | Rhodotorula copper AUC x MS2 features, 264 strains x 2 tracks | Rhodotorula_MS2_pheno_explore | 0/7341 (cell), 0/7333 (SUP) FDR hits; power floor rho=0.34; holdout 0/25 both tracks | null (calibrated as underpowered, not disproven) |
| 2026-08-11 | copper_metabolite_association v1 | same | Rhodotorula_MS2_pheno_explore | PLS CV R² negative both tracks but permutation p=0.017 (cell) / 0.020 (SUP) -- less-bad-than-chance | weak positive (refines F-001: some joint signal may exist despite no univariate hit) |
| 2026-08-11 | copper_metabolite_association v1 | same | Rhodotorula_MS2_pheno_explore | RF: cell p=0.030 (significant), supernatant p=0.199 (not significant) | partial corroboration of F-004 (cell only) |

## F-006: The cell-fraction metabolite signal and the genome-side methionine signal are not carried by the same strains
**Status:** preliminary
**Claim:** `09_cross_reference_genome_signal.py` tested the direct question this whole
analysis was built toward: do strains whose cell-fraction metabolite profile better
predicts copper AUC (out-of-fold PLS/RF prediction from F-004/F-005, in
species-corrected rank-residual space) also carry the elevated genome-wide methionine
usage found independently in `Rhodotorula_Rodeo` (`metal-resistance.md` F-002)? Joined
on normalized strain name (241/264 cell-track strains matched to Rodeo's 251-strain
genome-Met cohort; a sanity check confirmed both repos' independently-loaded copies of
`mean_auc_rate` agree exactly for the joined rows, ruling out a join error). Spearman
correlation between the out-of-fold metabolite prediction and the genome Met residual:
PLS ρ=+0.077 (permutation p=0.240), Random Forest ρ=+0.050 (permutation p=0.438) — both
null (n=241).
**Implications:** The weak cell-fraction metabolite signal (F-004/F-005) and the
genome-side Met signal (Rodeo F-002) are two independent, weak findings, not two views
of one shared underlying biology — this argues against a simple "high-Met strains also
have the distinctive metabolite profile" story. This is itself a real, informative
result (not a failed analysis): it's exactly the kind of check that keeps two
independently-discovered weak signals from being conflated into a false composite
narrative. Neither signal is individually strong enough that this null cross-check
should be read as "disproving" either one on its own.
**Open questions:** Would a more exploratory link (correlating specific metabolite
features directly against Met usage, bypassing the copper-AUC-prediction step) find
anything the AUC-mediated test above would miss? Is there a genome-side signal in a
different candidate family (not Met, e.g. Cu-oxidase-adjacent, even though those didn't
survive species correction in Rodeo's F-002) that tracks the metabolite signal instead?
**Tags:** copper, metabolomics, cross-repo, genome-metabolite-link, negative-result, rhodotorula

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-11 | copper_metabolite_association v1 | Rhodotorula copper AUC x MS2 features, 264 strains x 2 tracks | Rhodotorula_MS2_pheno_explore | 0/7341 (cell), 0/7333 (SUP) FDR hits; power floor rho=0.34; holdout 0/25 both tracks | null (calibrated as underpowered, not disproven) |
| 2026-08-11 | copper_metabolite_association v1 | same | Rhodotorula_MS2_pheno_explore | PLS CV R² negative both tracks but permutation p=0.017 (cell) / 0.020 (SUP) -- less-bad-than-chance | weak positive (refines F-001: some joint signal may exist despite no univariate hit) |
| 2026-08-11 | copper_metabolite_association v1 | same | Rhodotorula_MS2_pheno_explore | RF: cell p=0.030 (significant), supernatant p=0.199 (not significant) | partial corroboration of F-004 (cell only) |
| 2026-08-11 | copper_metabolite_association v1 | cell-track metabolite predictions x Rhodotorula_Rodeo genome Met residual, 241 common strains | Rhodotorula_MS2_pheno_explore + Rhodotorula_Rodeo (cross-repo) | PLS rho=+0.077 p=0.240; RF rho=+0.050 p=0.438 | null (F-004/F-005 metabolite signal and Rodeo F-002 genome signal are independent) |

### Open Questions
- What features drive the cell-track PLS/RF joint signal (loadings/importances not yet
  cross-referenced against each other or against F-001's top nominal univariate hits)?
- Would a chemistry-motivated targeted search (mirroring the successful torularhodin
  pathway-targeted approach for color phenotype) do better here than the untargeted scan?
- Would a dedicated holdout for the multivariate (not univariate) signal replicate the
  cell-track result out of sample?
- Would a more direct (non-AUC-mediated) metabolite-vs-Met correlation, or a check
  against other genome candidate families, find a link that F-006's test missed?
