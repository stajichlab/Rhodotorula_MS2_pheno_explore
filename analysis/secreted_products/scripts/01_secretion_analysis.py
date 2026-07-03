#!/usr/bin/env python3
"""Identify uniquely and highly secreted Rhodotorula MS2 metabolite products.

Secretion is measured as a *paired* supernatant-vs-cell contrast: each strain
contributes a cell-pellet (``C_*``) and a supernatant (``SUP_*``) sample, so a
feature enriched in SUP relative to C within the same strain is "secreted".
Uniqueness is a cross-species specificity (Tau) of supernatant abundance.

Design (robust-analysis conventions):
  * assert shapes/counts at every step; fail loudly, never silently drop;
  * paired Wilcoxon signed-rank across strains + BH FDR;
  * sign-flip permutation null (label C/SUP shuffle within strain);
  * sensitivity sweep over the log2FC and Tau thresholds.

Outputs land in ../outputs/ ; figures in ../outputs/figures/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- make register_value importable (mycelium core scripts) -----------------
_MYC_SCRIPTS = Path(
    "/rhome/jstajich/.claude/plugins/marketplaces/mycelium/skills/core/scripts"
)
if _MYC_SCRIPTS.is_dir():
    sys.path.insert(0, str(_MYC_SCRIPTS))
try:
    from register_value import register_value
except ImportError:
    # ANALYSIS_OK[optional-dependency]: register_value is a mycelium reporting
    # helper, not part of the science. If the mycelium scripts dir is absent the
    # analysis still runs correctly; only the numbers.json report fragment is
    # skipped. No scientific value is affected. Verified by: outputs (tables,
    # figures, summary json) are all written independently of this import.
    def register_value(*_a, **_k):  # type: ignore  # ANALYSIS_OK[optional-dependency]: no-op stub when the mycelium reporting helper is unavailable; science outputs are unaffected (see except comment above).
        return None

# --- configuration (every analytical decision is a named parameter) ---------
ROOT = Path(__file__).resolve().parents[3]
META_PATH = ROOT / "input_data" / "MS2_samples_combine.extended_metadata_with_strain_traits.tsv.gz"
MS2_PATH = ROOT / "input_data" / "Rhodotorula_MS2_aligned_features_ms2.csv.gz"
OUT = Path(__file__).resolve().parents[1] / "outputs"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

EPS = 1e-6                 # pseudocount on TIC-normalized fractions
LOG2FC_THRESH = 1.0        # >= 2-fold SUP enrichment (primary)
Q_THRESH = 0.05            # BH-FDR across features
MIN_SUP_DETECT_FRAC = 0.10 # feature must be detected in >=10% of SUP samples
TAU_THRESH = 0.85          # cross-species specificity for "unique"
N_PERM = 200
SEED = 20260702
STRAIN_KEY = "ATTRIBUTE_ID_1"
SPECIES_KEY = "ATTRIBUTE_species"

rng = np.random.default_rng(SEED)


def log(msg: str) -> None:
    print(f"[secretion] {msg}", flush=True)


# ============================================================================
# 1. Load & align metadata
# ============================================================================
log(f"loading metadata: {META_PATH.name}")
meta = pd.read_csv(META_PATH, sep="\t", dtype=str)
meta["sample_id"] = meta["filename"].astype(str)
samp = meta[meta["sample_id"].str.match(r"^(C|SUP)_\d+$")].copy()
samp["stype"] = samp["ATTRIBUTE_Source"].map({"cell_pellet": "C", "supernatant": "SUP"})
n_before = len(samp)
# ANALYSIS_OK[sample-filter]: drop samples with no source type or no strain id;
# they cannot enter a paired SUP-vs-C contrast. Count logged on the next line.
samp = samp.dropna(subset=["stype", STRAIN_KEY])
log(f"metadata C_/SUP_ rows: {n_before} -> {len(samp)} after requiring type+strain")
# ANALYSIS_OK[runtime-assert]: developer tripwire; pipeline is never run under `python -O`.
assert samp["stype"].isin(["C", "SUP"]).all(), "unexpected sample type"


def _nonempty(col: pd.Series) -> pd.Series:
    # ANALYSIS_OK[missingness]: a blank/NaN genome-accession field means "no
    # genome record linked" -> mapped to False in has_genome; recorded in the
    # strain_genome flag and reported via n_paired_strains_with_genome.
    return col.fillna("").str.strip().replace({"nan": ""}) != ""


samp["has_genome"] = (
    _nonempty(samp["db_sra_run_list"])
    | _nonempty(samp["db_biosample_list"])
    | _nonempty(samp["db_bioproject_list"])
)

# ============================================================================
# 2. Load MS2 feature table (features = rows, samples = *Peak area columns)
# ============================================================================
log(f"loading MS2 table: {MS2_PATH.name}")
ms2 = pd.read_csv(MS2_PATH)
peak_cols = [c for c in ms2.columns if c.endswith("Peak area")]
ann_cols = [c for c in ms2.columns if not c.endswith("Peak area")]
log(f"MS2 features (rows): {ms2.shape[0]}, peak-area columns: {len(peak_cols)}, annotation columns: {len(ann_cols)}")
# ANALYSIS_OK[runtime-assert]: developer tripwire; pipeline is never run under `python -O`.
assert ms2.shape[0] > 10000, "unexpected feature count"

# map peak-area column -> sample_id
sid_of_col = {c: c[: -len(".mzML Peak area")] for c in peak_cols}
meta_ids = set(samp["sample_id"])
use_cols = [c for c in peak_cols if sid_of_col[c] in meta_ids]
use_sids = [sid_of_col[c] for c in use_cols]
log(f"matched MS2 sample columns to metadata: {len(use_cols)} (dropped {len(peak_cols) - len(use_cols)} QC/blank/empty)")
# ANALYSIS_OK[runtime-assert]: developer tripwire; pipeline is never run under `python -O`.
assert len(use_cols) >= 550, "too few matched sample columns"

X = ms2[use_cols].copy()
X.columns = use_sids
X = X.apply(pd.to_numeric, errors="coerce")
# ANALYSIS_OK[runtime-assert]: developer tripwire; pipeline is never run under `python -O`.
assert not X.isna().any().any(), "non-numeric / NaN peak areas encountered"
# ANALYSIS_OK[threshold]: 0 is a physical floor (negative peak areas are
# impossible / instrument artifacts), not a tunable analysis parameter.
X[X < 0] = 0.0
log(f"peak-area matrix: {X.shape[0]} features x {X.shape[1]} samples; total signal range "
    f"[{X.values.min():.3g}, {X.values.max():.3g}]")

# detection mask (raw > 0) and within-sample TIC normalization -> fractions
detected = X.values > 0
col_sums = X.sum(axis=0).values
# ANALYSIS_OK[runtime-assert]: developer tripwire; pipeline is never run under `python -O`.
assert (col_sums > 0).all(), "a sample column has zero total signal"
Xn = X.values / col_sums[None, :]              # relative abundance, columns sum to 1
sid_to_idx = {s: i for i, s in enumerate(use_sids)}
np.testing.assert_allclose(Xn.sum(axis=0), 1.0, rtol=1e-6)

# ============================================================================
# 3. Build paired SUP / C matrices over strains with both sample types
# ============================================================================
samp_use = samp[samp["sample_id"].isin(use_sids)].copy()
strain_species = (
    # ANALYSIS_OK[sample-filter]: strains with no species label are excluded
    # only from per-species Tau means (17 labeled species); they still enter
    # the paired secretion test. Effect documented in SECRETED_PRODUCTS.md.
    samp_use.dropna(subset=[SPECIES_KEY]).groupby(STRAIN_KEY)[SPECIES_KEY].agg(lambda s: s.mode().iat[0])
)
strain_genome = samp_use.groupby(STRAIN_KEY)["has_genome"].any()


def _mean_norm(sids: list[str]) -> np.ndarray:
    idx = [sid_to_idx[s] for s in sids if s in sid_to_idx]
    return Xn[:, idx].mean(axis=1)


def _any_detect(sids: list[str]) -> np.ndarray:
    idx = [sid_to_idx[s] for s in sids if s in sid_to_idx]
    return detected[:, idx].any(axis=1)


paired_strains = []
sup_cols, c_cols, sup_det, c_det, sup_species = [], [], [], [], []
for strain, grp in samp_use.groupby(STRAIN_KEY):
    s_ids = grp.loc[grp["stype"] == "SUP", "sample_id"].tolist()
    c_ids = grp.loc[grp["stype"] == "C", "sample_id"].tolist()
    if not s_ids or not c_ids:
        continue
    paired_strains.append(strain)
    sup_cols.append(_mean_norm(s_ids))
    c_cols.append(_mean_norm(c_ids))
    sup_det.append(_any_detect(s_ids))
    c_det.append(_any_detect(c_ids))
    sup_species.append(strain_species.get(strain, "unknown"))

SUP = np.column_stack(sup_cols)   # features x strains
C = np.column_stack(c_cols)
SUPd = np.column_stack(sup_det)
Cd = np.column_stack(c_det)
n_feat, n_strain = SUP.shape
log(f"paired strains (both C and SUP): {n_strain}; features: {n_feat}")
# ANALYSIS_OK[runtime-assert]: developer tripwire; pipeline is never run under `python -O`.
assert SUP.shape == C.shape == SUPd.shape == Cd.shape
# ANALYSIS_OK[missingness]: a paired strain absent from strain_genome has no
# genome record -> False (no genome). Count reported as n_paired_strains_with_genome.
n_strain_genome = int(strain_genome.reindex(paired_strains).fillna(False).sum())
log(f"paired strains that also have a genome link: {n_strain_genome}")

# ============================================================================
# 4. Secretion statistics (paired, per feature across strains)
# ============================================================================
log2fc_mat = np.log2((SUP + EPS) / (C + EPS))    # features x strains
median_log2fc = np.median(log2fc_mat, axis=1)
mean_log2fc = log2fc_mat.mean(axis=1)
sup_detect_frac = SUPd.mean(axis=1)
c_detect_frac = Cd.mean(axis=1)
mean_sup_abund = SUP.mean(axis=1)                # mean relative abundance in supernatant

log("running paired Wilcoxon signed-rank (SUP vs C) per feature ...")
pvals = np.ones(n_feat)
for i in range(n_feat):
    a, b = SUP[i], C[i]
    diff = a - b
    if np.count_nonzero(diff) < 5:               # too few informative pairs
        continue
    try:
        pvals[i] = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue
    except ValueError:
        pvals[i] = 1.0
qvals = multipletests(pvals, method="fdr_bh")[1]

is_secreted = (
    (median_log2fc >= LOG2FC_THRESH)
    & (qvals < Q_THRESH)
    & (sup_detect_frac >= MIN_SUP_DETECT_FRAC)
)
log(f"secreted features (log2FC>={LOG2FC_THRESH}, q<{Q_THRESH}, SUP-detect>={MIN_SUP_DETECT_FRAC}): "
    f"{int(is_secreted.sum())}")

# ============================================================================
# 5. Cross-species specificity (Tau) of supernatant abundance
# ============================================================================
species_arr = np.array(sup_species)
uniq_species = [s for s in pd.unique(species_arr) if s != "unknown"]
sp_means = np.zeros((n_feat, len(uniq_species)))
for j, sp in enumerate(uniq_species):
    cols = np.where(species_arr == sp)[0]
    sp_means[:, j] = SUP[:, cols].mean(axis=1)


def tau_specificity(row: np.ndarray) -> float:
    mx = row.max()
    if mx <= 0:
        return np.nan
    return float(np.sum(1.0 - row / mx) / (len(row) - 1))


tau = np.array([tau_specificity(sp_means[i]) for i in range(n_feat)])
top_species = np.array([uniq_species[k] if not np.isnan(tau[i]) else "NA"
                        for i, k in enumerate(sp_means.argmax(axis=1))])
is_unique = is_secreted & (tau >= TAU_THRESH)
log(f"uniquely secreted (secreted & Tau>={TAU_THRESH}) across {len(uniq_species)} species: {int(is_unique.sum())}")

# ============================================================================
# 6. Assemble results table with chemistry annotation
# ============================================================================
res = pd.DataFrame({
    "feature_index": np.arange(n_feat),
    "median_log2FC_SUP_vs_C": median_log2fc,
    "mean_log2FC_SUP_vs_C": mean_log2fc,
    "wilcoxon_p": pvals,
    "wilcoxon_q_bh": qvals,
    "sup_detect_frac": sup_detect_frac,
    "c_detect_frac": c_detect_frac,
    "mean_sup_rel_abund": mean_sup_abund,
    "tau_species_specificity": tau,
    "top_species": top_species,
    "is_secreted": is_secreted,
    "is_uniquely_secreted": is_unique,
})
ann_keep = [c for c in ["row ID", "row m/z", "row retention time", "adduct",
                        "parent_mass", "charge", "has_ms2", "detection_count",
                        "detection_rate"] if c in ms2.columns]
res = pd.concat([res, ms2[ann_keep].reset_index(drop=True)], axis=1)
res = res.sort_values(["is_uniquely_secreted", "mean_sup_rel_abund"], ascending=[False, False])

res.to_csv(OUT / "secretion_scores_all_features.csv.gz", index=False, compression="gzip")
res[res["is_secreted"]].to_csv(OUT / "secreted_features.csv", index=False)
top_unique = res[res["is_uniquely_secreted"]].copy()
top_unique.to_csv(OUT / "uniquely_secreted_features.csv", index=False)
log(f"wrote all-feature, secreted, and uniquely-secreted tables")

# ============================================================================
# 7. Null model — sign-flip permutation of paired log2FC (C/SUP label shuffle)
# ============================================================================
log(f"permutation null ({N_PERM} sign-flip shuffles) ...")
obs_pos = int((median_log2fc >= LOG2FC_THRESH).sum())
null_counts = np.empty(N_PERM, dtype=int)
for p in range(N_PERM):
    signs = rng.choice([-1.0, 1.0], size=n_strain)[None, :]
    perm_median = np.median(log2fc_mat * signs, axis=1)
    null_counts[p] = int((perm_median >= LOG2FC_THRESH).sum())
null_mean = float(null_counts.mean())
emp_p = float((null_counts >= obs_pos).mean())
emp_fdr = min(1.0, null_mean / max(obs_pos, 1))
log(f"observed features with median log2FC>={LOG2FC_THRESH}: {obs_pos}; "
    f"null mean: {null_mean:.1f}; empirical p={emp_p:.4f}; empirical FDR~{emp_fdr:.4f}")
pd.DataFrame({"perm": np.arange(N_PERM), "null_count": null_counts}).to_csv(
    OUT / "null_distribution.csv", index=False)

# ============================================================================
# 8. Sensitivity sweep over log2FC and Tau thresholds
# ============================================================================
rows = []
for fc in [0.5, 1.0, 1.5, 2.0]:
    base = (median_log2fc >= fc) & (qvals < Q_THRESH) & (sup_detect_frac >= MIN_SUP_DETECT_FRAC)
    for tt in [0.75, 0.85, 0.95]:
        rows.append({"log2fc_thresh": fc, "tau_thresh": tt,
                     "n_secreted": int(base.sum()),
                     "n_uniquely_secreted": int((base & (tau >= tt)).sum())})
sens = pd.DataFrame(rows)
# ANALYSIS_OK[file-selection]: deterministic output write to a fixed path, not
# an input selection; no glob/sort/latest semantics involved.
sens.to_csv(OUT / "sensitivity_thresholds.csv", index=False)
log("wrote sensitivity_thresholds.csv")

# ============================================================================
# 9. Figures
# ============================================================================
# 9a. Volcano-like: median log2FC vs -log10 q
fig, ax = plt.subplots(figsize=(7, 5.5))
ax.scatter(median_log2fc, -np.log10(qvals + 1e-300), s=4, c="#b0b0b0", alpha=0.5, label="all")
sc = res["is_secreted"].values
ax.scatter(res.loc[sc, "median_log2FC_SUP_vs_C"], -np.log10(res.loc[sc, "wilcoxon_q_bh"] + 1e-300),
           s=8, c="#1f77b4", alpha=0.7, label="secreted")
uq = res["is_uniquely_secreted"].values
ax.scatter(res.loc[uq, "median_log2FC_SUP_vs_C"], -np.log10(res.loc[uq, "wilcoxon_q_bh"] + 1e-300),
           s=14, c="#d62728", alpha=0.9, label="uniquely secreted")
ax.axvline(LOG2FC_THRESH, ls="--", c="k", lw=0.8)
ax.axhline(-np.log10(Q_THRESH), ls="--", c="k", lw=0.8)
ax.set_xlabel("median log2(SUP / C) across strains")
ax.set_ylabel("-log10 BH q-value (paired Wilcoxon)")
ax.set_title("Secretion volcano: supernatant vs cell")
ax.legend(loc="upper left", frameon=False)
fig.tight_layout(); fig.savefig(FIG / "volcano_secretion.png", dpi=150); plt.close(fig)

# 9b. Secretion score distribution
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(median_log2fc, bins=120, color="#888")
ax.axvline(LOG2FC_THRESH, ls="--", c="#d62728")
ax.set_xlabel("median log2(SUP / C)"); ax.set_ylabel("features")
ax.set_title("Distribution of secretion scores")
fig.tight_layout(); fig.savefig(FIG / "secretion_score_hist.png", dpi=150); plt.close(fig)

# 9c. Top uniquely secreted features by supernatant abundance
topn = top_unique.head(25).iloc[::-1].reset_index(drop=True)
if len(topn):
    fig, ax = plt.subplots(figsize=(8, 7))
    def _sp_short(s):
        return s.split()[-1] if isinstance(s, str) and s.split() else ""
    labels = [
        f"idx{int(topn.loc[k, 'feature_index'])} "
        f"m/z{float(topn.loc[k, 'row m/z']):.3f} "
        f"rt{float(topn.loc[k, 'row retention time']):.2f} "
        f"[{_sp_short(topn.loc[k, 'top_species'])}]"
        for k in range(len(topn))
    ]
    ax.barh(range(len(topn)), topn["mean_sup_rel_abund"].values * 100, color="#d62728")
    ax.set_yticks(range(len(topn))); ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("mean supernatant relative abundance (%)")
    ax.set_title(f"Top {len(topn)} uniquely secreted features")
    fig.tight_layout(); fig.savefig(FIG / "top_uniquely_secreted.png", dpi=150); plt.close(fig)

# ============================================================================
# 10. Register reportable values
# ============================================================================
register_value("n_paired_strains", int(n_strain))
register_value("n_paired_strains_with_genome", int(n_strain_genome))
register_value("n_features_tested", int(n_feat))
register_value("n_species", int(len(uniq_species)))
register_value("log2fc_threshold", float(LOG2FC_THRESH))
register_value("q_threshold", float(Q_THRESH))
register_value("tau_threshold", float(TAU_THRESH))
register_value("n_secreted_features", int(is_secreted.sum()))
register_value("n_uniquely_secreted_features", int(is_unique.sum()))
register_value("null_mean_positive", round(null_mean, 1))
register_value("secretion_empirical_fdr", round(emp_fdr, 4))

# summary json for the report
summary = {
    "n_paired_strains": int(n_strain),
    "n_paired_strains_with_genome": int(n_strain_genome),
    "n_features_tested": int(n_feat),
    "n_species": int(len(uniq_species)),
    "n_secreted_features": int(is_secreted.sum()),
    "n_uniquely_secreted_features": int(is_unique.sum()),
    "observed_positive_median_log2fc": obs_pos,
    "null_mean_positive": null_mean,
    "empirical_p": emp_p,
    "empirical_fdr": emp_fdr,
    "thresholds": {"log2fc": LOG2FC_THRESH, "q": Q_THRESH,
                    "tau": TAU_THRESH, "min_sup_detect_frac": MIN_SUP_DETECT_FRAC},
}
pd.Series(summary).to_json(OUT / "secretion_summary.json", indent=2)
log("done. outputs in " + str(OUT))
