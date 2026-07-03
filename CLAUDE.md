# CLAUDE.md — Mycelium Living Repository

This repository is a **mycelium-enabled living repository**. It carries its own memory and grows smarter over time through structured traces of every action.

## Project Context

*Rhodotorula* MS2 metabolomics ↔ phenotype exploration. Core data: an aligned MS2
feature table (16,332 features across paired cell `C_*` and supernatant `SUP_*`
samples, ~295 of each) plus color phenotypes (CIELab L*/a*/b*) and strain metadata
for multiple *Rhodotorula* species. Downstream goal: link secreted/cell metabolites
to phenotypes and to candidate biosynthetic gene families in the genome datasets.

> **Legacy layout note.** This repo predates mycelium and keeps its own working
> directories alongside the canonical mycelium ones. Do **not** assume the mycelium
> dirs are authoritative for existing work:
> - Raw inputs live in `input_data/` (not `data/raw/`).
> - Phased pipeline scripts live in `scripts/`; phase outputs in `analysis/` and `results/`.
> - The canonical mycelium `data/`, `algorithms/`, `reference_material/` dirs were
>   created empty by init — new mycelium-managed work should use them, but old work
>   stays where it is. See `DATA_DESCRIPTION.md` and `analysis/README.md` for the
>   existing phase 0–3 pipeline.

## Quick Orientation

1. **Read `.living/INDEX.md` first** — A one-screen knowledge map: tag clusters, most-recent entries, and a tag → entry-ID inverted index. Drill into the underlying file (`learnings.md`, `decisions.md`, `conventions.md`) only when a row matches your task. The SessionStart hook regenerates this index every fresh session — trust it.
2. **For targeted lookup:** `python3 <mycelium>/skills/core/scripts/recall_lessons.py --living-dir .living/ --tag <tag>` — fetches only the matching entries instead of pulling the whole file. Also accepts `--id L-42`, `--since YYYY-MM-DD`.
3. **Read `ENVIRONMENTS_INSTALLATIONS.md`** — Environment setup, dependencies, and installation gotchas.
4. **Read the relevant manifest** — Each top-level directory has a descriptive manifest (`ANALYSIS_MANIFEST.md`, `DATA_MANIFEST.md`, `ALGORITHM_MANIFEST.md`, `REFERENCE_MANIFEST.md`).

> **Tooling note.** The mycelium scripts require Python ≥ 3.11 (`datetime.UTC`). The
> cluster default `python3` is 3.9 — use `/usr/bin/python3.12` for mycelium scripts.
> The mycelium install lives at `/rhome/jstajich/.claude/plugins/marketplaces/mycelium`.

## Installed Convention Packs

<!-- Check .living/conventions/ACTIVE_CONVENTIONS.yaml for the full list -->

### Core (auto-installed)

- **robust-analysis** — Defensive execution practices: strict error handling, data validation checks, parameter sensitivity sweeps, null hypothesis testing, adversarial self-challenge. See `.living/conventions/robust-analysis/analysis-conventions.md` for the entry point.
- **report-generator** — Structured LaTeX PDF report generation. See `.living/conventions/report-generator/analysis-conventions.md` for the workflow.
- **idea-generator** — Persona-based creative ideation for new analysis directions. See `.living/conventions/idea-generator/analysis-conventions.md` for the entry point.

### Domain (opt-in)

- **bioinformatics** — Domain-specific methodology for omics/sequence/feature-table
  analysis. See `.living/conventions/bioinformatics/analysis-conventions.md`. Layers on
  top of robust-analysis. Relevant here for metabolite feature filtering, secretion
  contrasts, and linking metabolites to biosynthetic gene families.

## Repository Structure (mycelium canonical)

```
├── algorithms/         — Reusable computational methods (see ALGORITHM_MANIFEST.md)
├── analysis/           — Analytical work (see ANALYSIS_MANIFEST.md)
├── data/               — Data assets: raw (immutable), processed, metadata (see DATA_MANIFEST.md)
├── reference_material/ — External references (see REFERENCE_MANIFEST.md)
├── todo/               — Future work items and ideas (see todo/TODO_REGISTRY.md)
└── .living/            — Repository memory layer
    ├── INDEX.md                    — Auto-regenerated knowledge map (read first)
    ├── decisions.md
    ├── learnings.md
    ├── conventions.md              — Repo-specific overrides
    ├── conventions/                — Installed convention packs
    │   └── ACTIVE_CONVENTIONS.yaml
    ├── generated-conventions/      — Conventions crystallized from learnings
    ├── log/                        — Append-only event log (LOG_REGISTRY.md)
    ├── findings/                   — Scientific findings by topic
    └── outputs/                    — Derived reports and transfer logs
```

## Workflow

### Before Starting Work

1. **Open `.living/INDEX.md`** — its tag clusters and recent-entries lists tell you which past learnings/decisions are likely relevant.
2. **Drill in selectively** — fetch just the entries with `recall_lessons.py --tag <tag>` / `--id L-42`, or read the whole file only if you need broader context.
3. Read the manifest for the area you'll be working in (e.g., `ANALYSIS_MANIFEST.md`).
4. Check `.living/conventions.md` for any repo-specific overrides.
5. If a domain convention is active, read its conventions in `.living/conventions/[domain]/`.

### While Working

- Follow analysis conventions: every analysis gets its own folder with UPPER_SNAKE_CASE.md documentation, scripts, outputs, reports.
- Follow statistical conventions: report effect sizes, confidence intervals, document assumptions.
- **Follow robust-analysis conventions**: fail loudly on unexpected data, assert shapes/types/ranges, log row counts at every step, run sensitivity analyses, test null hypotheses via permutation/bootstrap.
- **Do not subset data** without explicit user confirmation and justification.
- Every analysis must have a `run.sh` or `run.py` that reproduces final outputs.
- **Register reportable values** with `register_value.py` at the site that computes them.
- **Lint analysis code with `scilintr`** after writing/editing; drive findings to zero (fix or structured waiver).

### After Every Significant Action (Post-Action Hook Protocol)

1. **Update manifests** with new/changed entries.
2. **Update documentation** (UPPER_SNAKE_CASE.md in the affected subfolder).
3. **Log decisions** to `.living/decisions.md` for non-obvious choices.
4. **Log learnings** to `.living/learnings.md` (set `mitigation_type`).
5. **Log findings** to `.living/findings/{topic}.md` for scientific results.
6. **Log todos** to `todo/TODO_REGISTRY.md`.
7. **Validate structure**: `validate_structure.py`.
8. **Crystallize conventions** when 3+ learnings share a pattern.
9. **Convention feedback**: note whether pack practices helped or had gaps.
10. **Session summary**: write `.claude/last-session.md`.

### Automated Enforcement

Mycelium hooks are installed in `.claude/settings.local.json` (SessionStart → health,
PostToolUse → post-action / activity / read / data trackers, Stop → stop-check /
data-lineage). They enforce the post-action protocol after analysis work; read-only and
config-only sessions are never interrupted.

## Data Conventions

- `data/raw/` is **IMMUTABLE** — never modify original files. (Legacy raw inputs in
  `input_data/` are likewise immutable.)
- Every dataset has metadata in `data/metadata/[dataset-name]/` and a manifest entry.
- Large files are gitignored with download instructions documented.

## Key Files

| File | Purpose |
|------|---------|
| `ENVIRONMENTS_INSTALLATIONS.md` | How to set up the environment |
| `DATA_DESCRIPTION.md` / `analysis/README.md` | Existing phase 0–3 pipeline docs |
| `.living/INDEX.md` | Auto-regenerated knowledge map — entry point for `.living/` |
| `.living/decisions.md` | Why choices were made |
| `.living/learnings.md` | Accumulated insights and gotchas |
| `.living/conventions.md` | Repo-specific convention overrides |
| `.living/conventions/ACTIVE_CONVENTIONS.yaml` | Registry of installed convention packs |
| `todo/TODO_REGISTRY.md` | Master list of future work items |
| `*/*_MANIFEST.md` | Registry of contents in each directory |
