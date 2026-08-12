#!/usr/bin/env python3
"""Build the pigment-pathway candidate compound list used by
analysis/pathway_targeted_association/. Masses are computed from molecular
formulas using exact monoisotopic atomic masses (not hand-typed literature
values), so a formula typo is easy to catch and the precision is consistent.

Compound set and provenance:
  - Carotenoid pathway (phytoene -> torularhodin): the canonical fungal/yeast
    carotenogenesis pathway, well-documented for Rhodotorula/Rhodosporidium/
    Sporidiobolus and already the basis of docs/CAROTENOID_GENE_BRIGHTNESS_ANALYSIS.md's
    gene-family work in this project. High confidence these are the right
    compound class; formulas are standard (e.g. Britton, Liaaen-Jensen &
    Pfander, "Carotenoids Handbook").
  - DOPA-melanin precursors: some Rhodotorula/Rhodosporidium produce
    tyrosinase-mediated DOPA-melanin under stress; included per user's request
    to broaden scope beyond carotenoids. Lower confidence than carotenoids for
    THIS species/dataset -- flagged as such in the `confidence` column.
  - Mycosporine-like amino acids (MAAs): UV-protective compounds, well
    documented in algae/cyanobacteria/lichenized fungi; much less established
    in Rhodotorula specifically. Included per user's request; lowest
    confidence tier, flagged accordingly.

This is a hypothesis-generation list, not a validated standards library --
matches below should be treated as candidates for RT/MS2 follow-up, not
confirmed identifications.
"""
import csv
import os

# Exact atomic monoisotopic masses (IUPAC 2021 values)
ATOMIC_MASS = {"C": 12.0000000, "H": 1.00782503207, "N": 14.0030740048, "O": 15.9949146196}
PROTON_MASS = 1.007276467
ELECTRON_MASS = 0.00054858
NA_MASS = 22.98976928
K_MASS = 38.96370649
NH3_MASS = 14.0030740048 + 3 * 1.00782503207
H2O_MASS = 2 * 1.00782503207 + 15.9949146196

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def formula_mass(formula: dict) -> float:
    return sum(ATOMIC_MASS[el] * n for el, n in formula.items())


# name, formula (as element:count dict), category, pathway_step, confidence, notes
COMPOUNDS = [
    # --- Carotenoid pathway (phytoene -> torularhodin), high confidence ---
    ("phytoene", {"C": 40, "H": 64}, "carotenoid", "PSY product (step 1)", "high",
     "first committed carotenoid; colorless"),
    ("phytofluene", {"C": 40, "H": 62}, "carotenoid", "desaturation step 2", "high", ""),
    ("zeta-carotene", {"C": 40, "H": 60}, "carotenoid", "desaturation step 3", "high", ""),
    ("neurosporene", {"C": 40, "H": 58}, "carotenoid", "desaturation step 4", "high", ""),
    ("lycopene", {"C": 40, "H": 56}, "carotenoid", "desaturation step 5 (linear)", "high",
     "isomeric with gamma-/beta-carotene; cannot distinguish by mass alone"),
    ("gamma-carotene", {"C": 40, "H": 56}, "carotenoid", "monocyclization (CrtYB)", "high",
     "isomeric with lycopene/beta-carotene; cannot distinguish by mass alone"),
    ("beta-carotene", {"C": 40, "H": 56}, "carotenoid", "bicyclization (CrtYB)", "high",
     "isomeric with lycopene/gamma-carotene; cannot distinguish by mass alone; major yellow-orange pigment"),
    ("torulene", {"C": 40, "H": 54}, "carotenoid", "beta-carotene + extra desaturation", "high",
     "diagnostic Rhodotorula/Sporidiobolus pigment (orange-red)"),
    ("torularhodin", {"C": 40, "H": 52, "O": 2}, "carotenoid", "torulene oxidation (terminal)", "high",
     "diagnostic Rhodotorula/Sporidiobolus pigment (red); carboxylic acid, may prefer negative mode / [M-H]-"),
    # --- DOPA-melanin precursors, medium/lower confidence for this species ---
    ("L-DOPA", {"C": 9, "H": 11, "N": 1, "O": 4}, "melanin_precursor", "tyrosinase substrate", "medium", ""),
    ("dopaquinone", {"C": 9, "H": 9, "N": 1, "O": 4}, "melanin_precursor", "tyrosinase product", "low",
     "highly reactive/unstable, unlikely to be detected intact"),
    ("5,6-dihydroxyindole", {"C": 8, "H": 7, "N": 1, "O": 2}, "melanin_precursor", "DHI, downstream of DOPA", "medium", ""),
    ("DHICA", {"C": 9, "H": 7, "N": 1, "O": 4}, "melanin_precursor", "5,6-dihydroxyindole-2-carboxylic acid", "medium", ""),
    # --- Mycosporine-like amino acids, low confidence for this species ---
    ("mycosporine-glycine", {"C": 11, "H": 17, "N": 1, "O": 6}, "maa", "MAA core", "low", ""),
    ("shinorine", {"C": 15, "H": 23, "N": 2, "O": 8}, "maa", "MAA, glycine+serine substituted", "low", ""),
    ("porphyra-334", {"C": 14, "H": 22, "N": 2, "O": 8}, "maa", "MAA, glycine+threonine substituted", "low", ""),
    ("palythine", {"C": 10, "H": 16, "N": 2, "O": 5}, "maa", "MAA core", "low", ""),
]

ADDUCTS = [
    # (adduct_label, delta_mass_fn given monoisotopic M and target integer charge)
    ("[M+H]1+", lambda m: m + PROTON_MASS, 1),
    ("[M+Na]1+", lambda m: m + NA_MASS - ELECTRON_MASS, 1),
    ("[M+2H]2+", lambda m: (m + 2 * PROTON_MASS) / 2, 2),
    ("[M-H2O+H]1+", lambda m: m - H2O_MASS + PROTON_MASS, 1),
]

rows = []
for name, formula, category, step, confidence, notes in COMPOUNDS:
    formula_str = "".join(f"{el}{n}" for el, n in formula.items())
    m = formula_mass(formula)
    for adduct_label, adduct_fn, charge in ADDUCTS:
        mz = adduct_fn(m)
        rows.append({
            "compound": name, "formula": formula_str, "monoisotopic_mass": round(m, 5),
            "category": category, "pathway_step": step, "confidence": confidence,
            "adduct": adduct_label, "charge": charge, "expected_mz": round(mz, 5),
            "notes": notes,
        })

out_path = os.path.join(OUT_DIR, "pigment_pathway_targets.csv")
with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"wrote {len(rows)} compound x adduct target rows -> {out_path}")
print(f"  {len(COMPOUNDS)} compounds x {len(ADDUCTS)} adducts")
for cat in ("carotenoid", "melanin_precursor", "maa"):
    n = sum(1 for c in COMPOUNDS if c[2] == cat)
    print(f"  {cat}: {n} compounds")
