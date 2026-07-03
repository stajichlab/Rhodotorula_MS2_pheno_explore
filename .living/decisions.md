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
