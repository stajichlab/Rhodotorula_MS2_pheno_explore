#!/usr/bin/env python3
"""Broader companion to build_target_list.py, for the wider SIRIUS profiling pass
requested after F-006's narrow (17-compound, 15 ppm) search turned out to have matched
none of its hypothesized identities. This does NOT replace or overwrite
pigment_pathway_targets.csv (F-006's original target list stays reproducible as-is);
it writes a separate pigment_pathway_targets_broad.csv.

Widening strategy (still chemistry-motivated, not a blind mass scan):
  - Same core carotenoid pathway (phytoene -> torularhodin) as before, PLUS common
    carotenoid ESTERS (Rhodotorula xanthophylls are frequently fatty-acid-esterified
    in vivo) with C14:0/C16:0/C16:1/C18:0/C18:1/C18:2 acyl groups on the two candidates
    that have a free hydroxyl in principle (torularhodin's carboxylic acid can form an
    ester too) -- esterification shifts mass substantially, which is exactly why the
    narrow search (unesterified masses only) could have missed real carotenoid-derived
    features.
  - Ergosterol pathway (ergosterol, ergosterol peroxide, zymosterol, lanosterol,
    episterol) -- the dominant fungal sterol pathway, plausible for any yeast MS2 dataset
    and a reasonable "what else is coloured/lipid-related" complement to carotenoids.
  - Broader melanin/quinone pathway: adds pyomelanin/homogentisate-pathway compounds
    (homogentisic acid, 4-maleylacetoacetate) alongside the existing DOPA-melanin set,
    since basidiomycete yeasts can use either melanin route.
  - Broader MAA set: adds palythinol, asterina-330, usujirene, mycosporine-serine.
  - Flavins (riboflavin, FAD, FMN) -- some yeasts produce visible flavin-derived
    yellow pigmentation; cheap to include and mechanistically distinct from the above.
  - Adducts widened from 4 to include the additional types actually observed in the raw
    table's `adduct` column ([M+NH3+H]+, [M+ACN+H]+, [M+K]+, [M+2H]2+ already included,
    [M+3H]3+ for the larger compounds).
  - ppm tolerance widened from 15 to 30 (still tighter than a truly blind scan; matches
    matching_scripts/02_match_targets_broad.py's own default).

See build_target_list.py for the atomic-mass constants and formula_mass() logic (kept
here as an independent copy rather than imported, so this script has no import-order
dependency on the other one and stays runnable standalone).
"""
import csv
import os

ATOMIC_MASS = {"C": 12.0000000, "H": 1.00782503207, "N": 14.0030740048, "O": 15.9949146196}
PROTON_MASS = 1.007276467
ELECTRON_MASS = 0.00054858
NA_MASS = 22.98976928
K_MASS = 38.96370649
NH3_MASS = 14.0030740048 + 3 * 1.00782503207
H2O_MASS = 2 * 1.00782503207 + 15.9949146196
ACN_MASS = 2 * 12.0 + 3 * 1.00782503207 + 14.0030740048  # CH3CN

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def formula_mass(formula: dict) -> float:
    return sum(ATOMIC_MASS[el] * n for el, n in formula.items())


def add_formula(base: dict, extra: dict) -> dict:
    out = dict(base)
    for el, n in extra.items():
        out[el] = out.get(el, 0) + n
    return out


# fatty acyl groups for carotenoid ester enumeration: name -> formula of the ACID
# (esterification: alcohol -OH + acid -COOH -> ester + H2O, i.e. add acid formula
# minus H2O)
FATTY_ACIDS = {
    "C14:0(myristoyl)": {"C": 14, "H": 28, "O": 2},
    "C16:0(palmitoyl)": {"C": 16, "H": 32, "O": 2},
    "C16:1(palmitoleoyl)": {"C": 16, "H": 30, "O": 2},
    "C18:0(stearoyl)": {"C": 18, "H": 36, "O": 2},
    "C18:1(oleoyl)": {"C": 18, "H": 34, "O": 2},
    "C18:2(linoleoyl)": {"C": 18, "H": 32, "O": 2},
}

COMPOUNDS = [
    # --- carotenoid pathway core (same as narrow list), high confidence ---
    ("phytoene", {"C": 40, "H": 64}, "carotenoid", "PSY product (step 1)", "high", ""),
    ("phytofluene", {"C": 40, "H": 62}, "carotenoid", "desaturation step 2", "high", ""),
    ("zeta-carotene", {"C": 40, "H": 60}, "carotenoid", "desaturation step 3", "high", ""),
    ("neurosporene", {"C": 40, "H": 58}, "carotenoid", "desaturation step 4", "high", ""),
    ("lycopene", {"C": 40, "H": 56}, "carotenoid", "desaturation step 5 (linear)", "high", "isomeric with gamma-/beta-carotene"),
    ("gamma-carotene", {"C": 40, "H": 56}, "carotenoid", "monocyclization", "high", "isomeric with lycopene/beta-carotene"),
    ("beta-carotene", {"C": 40, "H": 56}, "carotenoid", "bicyclization", "high", "isomeric with lycopene/gamma-carotene"),
    ("torulene", {"C": 40, "H": 54}, "carotenoid", "extra desaturation", "high", "diagnostic Rhodotorula pigment"),
    ("torularhodin", {"C": 40, "H": 52, "O": 2}, "carotenoid", "torulene oxidation", "high", "diagnostic Rhodotorula pigment; refuted for raw_position 11564/12635 by SIRIUS, still valid to search elsewhere"),
    # --- ergosterol pathway, medium-high confidence (dominant fungal sterol) ---
    ("ergosterol", {"C": 28, "H": 44, "O": 1}, "sterol", "major fungal membrane sterol", "high", ""),
    ("ergosterol_peroxide", {"C": 28, "H": 44, "O": 3}, "sterol", "oxidized ergosterol", "medium", ""),
    ("zymosterol", {"C": 27, "H": 44, "O": 1}, "sterol", "ergosterol pathway intermediate", "medium", ""),
    ("lanosterol", {"C": 30, "H": 50, "O": 1}, "sterol", "ergosterol pathway intermediate", "medium", ""),
    ("episterol", {"C": 28, "H": 46, "O": 1}, "sterol", "ergosterol pathway intermediate", "medium", ""),
    # --- DOPA-melanin precursors (same as narrow list) ---
    ("L-DOPA", {"C": 9, "H": 11, "N": 1, "O": 4}, "melanin_precursor", "tyrosinase substrate", "medium", ""),
    ("5,6-dihydroxyindole", {"C": 8, "H": 7, "N": 1, "O": 2}, "melanin_precursor", "DHI", "medium", ""),
    ("DHICA", {"C": 9, "H": 7, "N": 1, "O": 4}, "melanin_precursor", "5,6-dihydroxyindole-2-carboxylic acid", "medium", ""),
    # --- pyomelanin / homogentisate pathway, medium confidence ---
    ("homogentisic_acid", {"C": 8, "H": 8, "O": 4}, "melanin_precursor", "pyomelanin precursor (HGA pathway)", "medium", ""),
    ("4-maleylacetoacetate", {"C": 8, "H": 8, "O": 6}, "melanin_precursor", "HGA pathway intermediate", "low", ""),
    # --- MAAs, low confidence (broadened set) ---
    ("mycosporine-glycine", {"C": 11, "H": 17, "N": 1, "O": 6}, "maa", "MAA core", "low", ""),
    ("shinorine", {"C": 15, "H": 23, "N": 2, "O": 8}, "maa", "MAA, glycine+serine", "low", ""),
    ("porphyra-334", {"C": 14, "H": 22, "N": 2, "O": 8}, "maa", "MAA, glycine+threonine", "low", ""),
    ("palythine", {"C": 10, "H": 16, "N": 2, "O": 5}, "maa", "MAA core", "low", ""),
    ("palythinol", {"C": 11, "H": 18, "N": 2, "O": 6}, "maa", "hydroxylated palythine", "low", ""),
    ("asterina-330", {"C": 13, "H": 20, "N": 2, "O": 6}, "maa", "MAA, alanine substituted", "low", ""),
    ("usujirene", {"C": 13, "H": 20, "N": 2, "O": 6}, "maa", "MAA, threonine-derived isomer", "low", ""),
    ("mycosporine-serine", {"C": 12, "H": 19, "N": 1, "O": 8}, "maa", "MAA core, serine-substituted", "low", ""),
    # --- flavins, medium confidence (visible yellow pigments) ---
    ("riboflavin", {"C": 17, "H": 20, "N": 4, "O": 6}, "flavin", "vitamin B2, yellow pigment", "medium", ""),
    ("FMN", {"C": 17, "H": 21, "N": 4, "O": 9, "S": 0}, "flavin", "flavin mononucleotide", "medium", ""),
]
# FMN actually contains P not S -- fix formula explicitly (kept as a visible correction
# rather than silently right the first time, since P wasn't in ATOMIC_MASS yet).
ATOMIC_MASS["P"] = 30.97376199
ATOMIC_MASS["S"] = 31.97207100
COMPOUNDS[-1] = ("FMN", {"C": 17, "H": 21, "N": 4, "O": 9, "P": 1}, "flavin", "flavin mononucleotide", "medium", "")
COMPOUNDS.append(("FAD", {"C": 27, "H": 33, "N": 9, "O": 15, "P": 2}, "flavin", "flavin adenine dinucleotide", "medium", ""))

# --- carotenoid fatty-acid esters: torularhodin (has a -COOH, forms esters too via
# acyl on a putative hydroxylated derivative is speculative; more defensibly, xanthophyll
# esters need a free -OH which none of the core list has except via hydroxylation -- so
# instead enumerate esters on a hydroxy-torulene hypothetical (C40H54O, +O relative to
# torulene, a plausible hydroxylation product) as the most defensible ester scaffold. ---
COMPOUNDS.append(("hydroxy-torulene", {"C": 40, "H": 54, "O": 1}, "carotenoid", "hypothetical hydroxylated torulene", "low", "not a confirmed natural product; included only as an ester-formation scaffold"))
hydroxy_torulene_formula = {"C": 40, "H": 54, "O": 1}
for acyl_name, acid_formula in FATTY_ACIDS.items():
    ester_formula = add_formula(hydroxy_torulene_formula, acid_formula)
    ester_formula["O"] -= 1  # lose one O and 2 H as H2O on esterification
    ester_formula["H"] -= 2
    COMPOUNDS.append((
        f"hydroxy-torulene-{acyl_name}-ester", ester_formula, "carotenoid",
        "hypothetical xanthophyll ester", "low",
        "speculative ester scaffold, not a confirmed natural product",
    ))

ADDUCTS = [
    ("[M+H]1+", lambda m: m + PROTON_MASS, 1),
    ("[M+Na]1+", lambda m: m + NA_MASS - ELECTRON_MASS, 1),
    ("[M+K]1+", lambda m: m + K_MASS - ELECTRON_MASS, 1),
    ("[M+2H]2+", lambda m: (m + 2 * PROTON_MASS) / 2, 2),
    ("[M+3H]3+", lambda m: (m + 3 * PROTON_MASS) / 3, 3),
    ("[M-H2O+H]1+", lambda m: m - H2O_MASS + PROTON_MASS, 1),
    ("[M+NH3+H]1+", lambda m: m + NH3_MASS + PROTON_MASS, 1),
    ("[M+ACN+H]1+", lambda m: m + ACN_MASS + PROTON_MASS, 1),
]

PPM_TOLERANCE_DEFAULT = 30  # documented here for reference; actual matching script owns the constant

rows = []
for name, formula, category, step, confidence, notes in COMPOUNDS:
    formula_str = "".join(f"{el}{n}" for el, n in formula.items() if n)
    m = formula_mass(formula)
    for adduct_label, adduct_fn, charge in ADDUCTS:
        mz = adduct_fn(m)
        rows.append({
            "compound": name, "formula": formula_str, "monoisotopic_mass": round(m, 5),
            "category": category, "pathway_step": step, "confidence": confidence,
            "adduct": adduct_label, "charge": charge, "expected_mz": round(mz, 5),
            "notes": notes,
        })

out_path = os.path.join(OUT_DIR, "pigment_pathway_targets_broad.csv")
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"wrote {len(rows)} compound x adduct target rows -> {out_path}")
print(f"  {len(COMPOUNDS)} compounds x {len(ADDUCTS)} adducts")
for cat in sorted(set(c[2] for c in COMPOUNDS)):
    n = sum(1 for c in COMPOUNDS if c[2] == cat)
    print(f"  {cat}: {n} compounds")
