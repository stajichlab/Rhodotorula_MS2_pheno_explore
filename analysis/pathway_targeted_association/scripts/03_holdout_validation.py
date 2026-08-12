#!/usr/bin/env python3
"""Species-stratified, strain-level 80/20 holdout replication for the targeted
pigment candidates -- same design as
analysis/phenotype_metabolite_association/scripts/04_replication_holdout.py, applied
to just the 3 pathway-matched features instead of a genome-wide scan. With only 3
features there is no cherry-picking/winner's-curse concern the way there was for the
"top nominal" candidates in the untargeted analysis, but holdout replication is still
the strongest available check and is cheap to run here.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))
OUT_DIR = os.path.join(SCRIPT_DIR, "..", "outputs")
PMA_OUT_DIR = os.path.join(ROOT, "analysis", "phenotype_metabolite_association", "outputs")

sys.path.insert(0, os.path.join(ROOT, "skills", "core", "scripts"))
try:
    from register_value import register_value
except ImportError:
    # ANALYSIS_OK[optional-dependency]: register_value is an internal mycelium
    # reporting helper, not a scientific dependency; numbers are still printed to
    # stdout/saved in outputs/*.json above, so a no-op fallback loses convenience
    # (report auto-fill) but not any result.
    def register_value(*args, **kwargs):
        pass

PHENOTYPE_COLS = ["Median_ColorLab_L*Mean", "Median_ColorLab_a*Mean", "Median_ColorLab_b*Mean"]
RAW_TABLE_PATH = os.path.join(ROOT, "input_data", "Rhodotorula_MS2_aligned_features_ms2.csv.gz")
RNG_SEED = 20260811
TEST_FRAC = 0.20
MIN_STRAINS_FOR_SPLIT = 5

print("=" * 80)
print("03: HOLDOUT REPLICATION FOR TARGETED PIGMENT CANDIDATES")
print("=" * 80)

targets = pd.read_csv(os.path.join(OUT_DIR, "pathway_target_feature_list.csv"))
design = pd.read_csv(os.path.join(PMA_OUT_DIR, "sample_design.csv"))
n = len(design)
rng = np.random.default_rng(RNG_SEED)

raw = pd.read_csv(RAW_TABLE_PATH)
sample_cols = [c for c in raw.columns if "Peak area" in c]
sample_total = raw[sample_cols].sum(axis=0)
candidate_raw = raw.iloc[targets["raw_position"].values][sample_cols]
candidate_raw.index = targets["raw_position"].values
sample_col_of = {c.replace(".mzML Peak area", ""): c for c in sample_cols}
design_cols = [sample_col_of[s] for s in design["sample_id"]]
abundance = candidate_raw[design_cols].T
abundance.index = design["sample_id"].values
totals = sample_total[design_cols].values
normalized = abundance.div(totals, axis=0)
pseudocount = 0.1 * normalized[normalized > 0].median().median()
log_abundance = np.log2(normalized + pseudocount)
log_abundance.index = design.index  # align to design row order

strain_species = design.drop_duplicates("strain_id").set_index("strain_id")["species"]
species_strain_counts = strain_species.value_counts()
splittable_species = species_strain_counts[species_strain_counts >= MIN_STRAINS_FOR_SPLIT].index
test_strains = set()
for sp in splittable_species:
    strains_sp = strain_species[strain_species == sp].index.tolist()
    n_test = max(1, round(len(strains_sp) * TEST_FRAC))
    test_strains.update(rng.choice(strains_sp, size=n_test, replace=False))
is_test = design["strain_id"].isin(test_strains).values
print(f"train rows: {(~is_test).sum()}  test rows: {is_test.sum()}  "
      f"({len(test_strains)} held-out strains from {len(splittable_species)} species)")


def fit(mask):
    sub_design = design.loc[mask].reset_index(drop=True)
    sub_log_abund = log_abundance.loc[mask].reset_index(drop=True)
    nn = len(sub_design)

    Z_parts = [np.ones((nn, 1))]
    for col in ("species", "Library Plate", "sample_type"):
        dummies = pd.get_dummies(sub_design[col].astype(str), drop_first=True)
        Z_parts.append(dummies.values.astype(float))
    Z = np.hstack(Z_parts)
    rank_Z = np.linalg.matrix_rank(Z)
    Z_pinv = np.linalg.pinv(Z)

    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against
    # Species/Library Plate/sample_type, identical method to
    # phenotype_metabolite_association/scripts/02 and this folder's script 02.
    def residualize(M):
        return M - Z @ (Z_pinv @ M)

    fr = sub_log_abund.apply(rankdata, axis=0).values
    # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical method to phenotype_metabolite_association/scripts/02_corrected_correlation.py.
    fr = residualize(fr)
    fr_c = fr - fr.mean(axis=0)
    fn = np.sqrt((fr_c ** 2).sum(axis=0))

    df = nn - rank_Z - 1
    out = {}
    for phenotype in PHENOTYPE_COLS:
        y = sub_design[phenotype].values.astype(float)
        # ANALYSIS_OK[batch-correction]: rank-space OLS residualization against Species/Library Plate/sample_type, identical method to phenotype_metabolite_association/scripts/02_corrected_correlation.py.
        y_resid = residualize(rankdata(y).reshape(-1, 1)).ravel()
        y_c = y_resid - y_resid.mean()
        rho = (fr_c.T @ y_c) / (fn * np.sqrt((y_c ** 2).sum()))
        t_stat = rho * np.sqrt(df / np.clip(1 - rho ** 2, 1e-12, None))
        pval = 2 * t_dist.sf(np.abs(t_stat), df=df)
        out[phenotype] = pd.DataFrame({"raw_position": log_abundance.columns, "rho": rho, "pval": pval})
    return out, df


train_fits, train_df = fit(~is_test)
test_fits, test_df = fit(is_test)
print(f"train df={train_df}  test df={test_df}")

rows = []
for phenotype in PHENOTYPE_COLS:
    tr = train_fits[phenotype].rename(columns={"rho": "rho_train", "pval": "pval_train"})
    te = test_fits[phenotype].rename(columns={"rho": "rho_test", "pval": "pval_test"})
    # ANALYSIS_OK[join]: raw_position is unique within each of tr/te (one row per candidate feature).
    merged = tr.merge(te, on="raw_position", validate="one_to_one")
    merged["phenotype"] = phenotype
    merged["compound"] = merged["raw_position"].map(
        targets.drop_duplicates("raw_position").set_index("raw_position")["compound"]
    )
    merged["same_sign"] = np.sign(merged["rho_train"]) == np.sign(merged["rho_test"])
    # ANALYSIS_OK[threshold]: nominal p<0.05 in the held-out test set, same
    # replication criterion as phenotype_metabolite_association/scripts/04 and 07.
    merged["replicated"] = merged["same_sign"] & (merged["pval_test"] < 0.05)
    rows.append(merged)

result = pd.concat(rows, ignore_index=True).sort_values("pval_train")
print("\n" + result[["phenotype", "raw_position", "compound", "rho_train", "pval_train",
                      "rho_test", "pval_test", "replicated"]].to_string(index=False))

# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
result.to_csv(os.path.join(OUT_DIR, "targeted_holdout_results.csv"), index=False)

n_replicated = int(result["replicated"].sum())
summary = {
    "n_test_strains": len(test_strains),
    "n_train_rows": int((~is_test).sum()),
    "n_test_rows": int(is_test.sum()),
    "n_pairs_tested": int(len(result)),
    "n_replicated": n_replicated,
    "rng_seed": RNG_SEED,
}
# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
with open(os.path.join(OUT_DIR, "targeted_holdout_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\n" + json.dumps(summary, indent=2))

# ANALYSIS_OK[file-selection]: os.path.join(OUT_DIR, ...) here writes a named output file (not a glob/latest read); no ambiguity about which file is used.
register_value("targeted_holdout_n_replicated", n_replicated, provenance="outputs/targeted_holdout_summary.json")
