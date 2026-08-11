---
topic: color-phenotype-metabolomics
description: Whether specific MS2 metabolite features are associated with CIELab color phenotype (L*/a*/b*) in Rhodotorula strains, and what it takes to trust such an association.
created: 2026-08-11
last_updated: 2026-08-11
status: active
---

# Color-phenotype metabolomics (L*/a*/b* vs MS2 features)

## F-001: The pooled metabolite-color correlations reported in docs/FEATURE_ANALYSIS.md are a species confound, not a real effect
**Status:** preliminary
**Claim:** `docs/FEATURE_ANALYSIS.md` (legacy Phase 0-2 pipeline) reported very strong
pooled Spearman correlations between MS2 features and color phenotype (e.g. Feature 2755,
ρ=0.735 with L*), attributing these to candidate pigment metabolites. A corrected
re-analysis (`analysis/phenotype_metabolite_association/`) that jointly regresses out
Species, Library Plate, and C/SUP sample type (the original pipeline never actually
controlled for Species despite its own Phase 0 decision recommending it, and had a
residualization bug) found **0 of the original 12,269 "Tier 1" hits survive** FDR
correction, a strain-block permutation null, or species-stratified holdout replication.
Feature 2755 specifically collapses from ρ=0.735 to ρ=0.024 (permutation p=0.69); all 6
headline features named in `FEATURE_ANALYSIS.md`, across all reported phenotypes (17
feature×phenotype pairs), fail permutation testing (p>0.05). Holdout replication of the
strongest *corrected* candidates was at chance (2/75, 2.7%).
**Implications:** No specific metabolite feature from this dataset/method can currently
be claimed to control Rhodotorula color phenotype. Any future correlation-based claim on
this dataset must (a) verify the confounding-variable control was actually applied in the
code path taken, not just intended, and (b) be checked against permutation + holdout
before being reported as a candidate.
**Tags:** metabolomics, statistics, confounding, simpsons-paradox, species, rhodotorula, correlation, color-phenotype, validation

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-11 | phenotype_metabolite_association v1 | Rhodotorula MS2 aligned features + CIELab phenotype (550 samples, 17 species, corrected) | Rhodotorula_MS2_pheno_explore | 0/12,269 original hits survive; max corrected \|ρ\|=0.19; 0/17 headline feature×phenotype pairs pass permutation; 2/75 holdout replication | contradicts (docs/FEATURE_ANALYSIS.md's original claim) |

### Open Questions
- Does a multivariate/module-level signal (PLS-DA, sparse CCA) exist even though no
  single feature clears univariate testing?
- Would a mixed-effects model with real phenotype replicates (see the analysis doc's
  "Framework for incorporating future phenotype replicates") recover a smaller but real
  effect currently swamped by measurement noise?

## F-002: Restricting to R. mucilaginosa alone (the one species with enough strains to test) does not recover a within-species signal, despite real phenotype spread
**Status:** preliminary
**Claim:** *R. mucilaginosa* (n=415 samples, 210 strains — the only species in this
dataset with enough strains for a within-species test) has genuine color-phenotype
variance to work with (a* CV=13%, range 0.6-14.4; b* CV=35%, range -0.3-8.7; L* is flatter,
CV=2.3%). Re-running the corrected model restricted to this species alone (Species
dropped from the covariate set, Library Plate + C/SUP sample type retained) still finds
**0 features clearing FDR<0.05**, with max |ρ|=0.196 — the same effect-size ceiling as
the pooled, between-species-confounded model.
**Implications:** The lack of signal in F-001 isn't simply an artifact of pooling weakly
across many small species groups; even the one species with substantial genetic/strain
diversity and real phenotype spread shows no detectable univariate metabolite-phenotype
association at this sample size and correction stringency. If a real genetic effect on
pigmentation exists, it's either smaller than this design can detect at n≈210 strains, or
not capturable as a single-feature linear/rank correlation.
**Tags:** metabolomics, rhodotorula-mucilaginosa, within-species, color-phenotype, genetics, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-11 | phenotype_metabolite_association v1 (script 05) | R. mucilaginosa subset (415 samples, 210 strains) | Rhodotorula_MS2_pheno_explore | 0 FDR hits; max \|ρ\|=0.196; phenotype spread confirmed non-trivial for a*/b* | supports (F-001's conclusion that the original signal was a species artifact) |

### Open Questions
- Is 210 strains enough power to detect a realistic single-gene pigmentation effect
  given typical Mendelian/QTL effect sizes, or is this an underpowered test?
- Would stratifying R. mucilaginosa further (e.g. by phylogenomic clade, as
  `docs/PHASE3_STRATIFIED_ANALYSIS_SUMMARY.md` suggested) reveal sub-structure that a
  single pooled within-species model still averages away?

## F-003: No joint (multivariate) MS2 signal for color phenotype either, and R. mucilaginosa's own holdout replicates at chance
**Status:** preliminary
**Claim:** Two further checks close out the obvious remaining alternatives to F-001/F-002.
(a) A PLS-regression module-level test (all 7,341 residualized features jointly
predicting residualized L*/a*/b*, cross-validated R² under strain-grouped k-fold, with
the component-count model-selection step itself repeated inside a strain-block
permutation null to avoid winner's-curse) finds no joint signal: observed CV R²=-0.09,
permutation null 95% range=[-0.14, -0.03], p=0.56. (b) A dedicated strain-level 80/20
holdout restricted to *R. mucilaginosa* alone (not pooled across species, unlike F-001's
holdout) replicates the species' own top nominal candidates at 3/75 (4.0%) — chance
level, with train/test effect-size calibration for a* actually negative (Spearman
ρ=-0.46).
**Implications:** Four independent checks (FDR, permutation, holdout, multivariate) now
agree on the null result. This is a reasonably thorough case that no single- or
joint-feature MS2 signal for color phenotype is detectable in this dataset at this
sample size, rather than a fragile conclusion resting on one test. Further work should
either bring in a different data type (transcriptomics, genotype/QTL) or wait for
additional phenotype replicates to reduce measurement noise (see the analysis doc's
mixed-model framework), rather than continuing to re-slice this same MS2 dataset.
**Tags:** metabolomics, multivariate, pls-regression, cross-validation, rhodotorula-mucilaginosa, color-phenotype, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-11 | phenotype_metabolite_association v1 (scripts 06, 07) | Rhodotorula MS2 (550 samples corrected; R. mucilaginosa 415-sample subset for holdout) | Rhodotorula_MS2_pheno_explore | PLS CV R²=-0.09, permutation p=0.56; mucilaginosa-only holdout replication 3/75 (4.0%) | supports (F-001/F-002's negative conclusion) |

### Open Questions
- What sample size would be needed to detect a plausible single-gene effect size with
  this design (a formal power analysis hasn't been run)?

## F-004: A non-linear multivariate model (Random Forest) finds no signal either — the null is robust to model class
**Status:** preliminary
**Claim:** F-003 asked whether a non-linear model might recover a joint signal that
linear PLS missed. Re-running the same permutation-tested module-level design with a
Random Forest (150 trees, `max_features='sqrt'`, `max_depth` grid-searched inside the
permutation loop, same strain-grouped CV) gives the same answer: observed CV R²=-0.026,
permutation null 95% range=[-0.067, -0.008], p=0.25.
**Implications:** The null result in F-001/F-002/F-003 is not an artifact of PLS's
linearity assumption — a tree-based model that can capture arbitrary feature
interactions and non-linear dose-response relationships still finds nothing
distinguishable from chance. Five independent checks (FDR, permutation, pooled holdout,
linear multivariate, non-linear multivariate) now converge on the same negative result.
**Tags:** metabolomics, random-forest, multivariate, cross-validation, non-linear, color-phenotype, negative-result

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-11 | phenotype_metabolite_association v1 (script 08) | Rhodotorula MS2 (550 samples, corrected) | Rhodotorula_MS2_pheno_explore | RF CV R²=-0.026, permutation p=0.25 | supports (F-001/F-002/F-003's negative conclusion) |

### Open Questions
- What sample size would be needed to detect a plausible single-gene effect size with
  this design (a formal power analysis hasn't been run) — this remains the key
  unanswered question after five converging negative checks.
