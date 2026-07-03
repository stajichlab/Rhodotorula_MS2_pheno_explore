# Species-Specific Metabolite Features — Candidate Identity Report

**Date:** 2026-07-02
**Source analysis:** `scripts/species_specific_features.py` → `results/species_specific/`
**Input:** `results/phase0/phase0_ms2_aligned.csv.gz` (16,332 features × 590 samples, 13 species with ≥3 samples), annotated against `input_data/Rhodotorula_MS2_aligned_features_ms2.csv.gz` (adduct/isotope/MS2 metadata)

## 1. Background

Of 16,332 aligned LC-MS features tested across 13 *Rhodotorula* species (≥3 samples each), only **two** passed strict species-specificity criteria (≥80% detection in the target species, ≤20% max detection in any other species, Fisher's exact FDR q < 0.05). This report characterizes those two features using precursor mass, charge/adduct, retention time, and detection pattern, and gives reasoned candidate compound classes. **No individual MS2 fragment ion list (m/z + intensity pairs) was available for annotation** — only precursor-level metadata (`has_ms2 = True`, scan counts) — so identities below are hypotheses to prioritize follow-up work, not confirmed structures.

## 2. Feature 1 — Row ID 20534 (*Rhodotorula sp.* clade XI–specific)

| Property | Value |
|---|---|
| m/z | 1360.1685 |
| Charge / adduct | 2+ / [M+2H]²⁺ |
| Isotope peak | M+1 (not monoisotopic) |
| Implied neutral mass | ≈2718.3 Da (M+1); monoisotopic ≈2717.3 Da |
| Retention time | 2.63 min (near-void — highly polar/hydrophilic) |
| MS2 | present, but only 9 scans; Gaussian similarity 0.96; noise score 0.22 |
| Detection | 6/590 samples total (1.0% overall) |
| Detected in | *R. sp.* clade XI: 2 cell-pellet (mean area 82,883) + 2 supernatant (mean area 2,459,679); *R. sp.* clade I: 1 supernatant (645,442). Zero elsewhere (15/17 other species). |

**Interpretation.** This is an unusually large (~2.7 kDa), highly polar, doubly-charged, almost exclusively **secreted (supernatant)** species detected in only 4 samples of one undescribed clade (plus one trace hit in a sister clade). Candidate classes, most to least plausible for a basidiomycete yeast exudate at this mass/polarity:

1. **Extracellular polysaccharide/exopolysaccharide oligomer or capsular glycan fragment** (medium–high plausibility) — *Rhodotorula* species are known producers of capsular/extracellular polysaccharides (mannans, glucuronoxylomannan-like glycans). ~2717 Da is consistent with an oligosaccharide of roughly 15–17 hexose/pentose units, matches the very early (near-void) retention time, and multiply-hydroxylated glycans can carry 2+ charge in positive mode.
2. **Extended glycolipid/glycosylated-ceramide conjugate** (medium) — possible but this mass is larger than typical yeast glycosphingolipids (usually <1500 Da), so would require an unusually long glycan chain.
3. **Cyclic glycopeptide or lipopeptide siderophore** (medium–low) — fungal siderophores (ferrichrome/coprogen types) are normally 500–800 Da; a ~20–24 residue peptide could reach this size but has no strong precedent in *Rhodotorula*.

**Caveat — likely artifact risk is non-trivial.** With only 6/590 detections, a thin MS2 scan count (9), and restriction to a single unnamed clade with just 4 replicate detections, this could alternatively reflect an in-source dimer/adduct misassignment (true monomer misannotated as 2+), a batch-specific media/culture contaminant, or a feature-picking artifact. Confidence that this is a genuine, clade-specific secondary metabolite is **low–medium** pending isotope-envelope confirmation of the 2+ charge state and independent replication.

## 3. Feature 2 — Row ID 47208 (*Rhodotorula taiwanensis*–enriched)

| Property | Value |
|---|---|
| m/z | 445.3146 |
| Charge / adduct | 1+ / [M−2H₂O+H]⁺ (non-default; software-assigned) |
| Parent (neutral) mass | 480.3281 Da |
| Isotope peak | M+0 (monoisotopic) |
| Retention time | 4.63 min (mid-polarity) |
| MS2 | present, 7 scans; Gaussian similarity 0.96; noise score 0.45 |
| Detection | 99/590 samples overall (16%) — a "leaky" background feature genus-wide, but concentrated in *R. taiwanensis* |
| *R. taiwanensis* detection | 10/12 samples (83%) — cell pellet (n=5, mean area 429,136) and supernatant (n=5, mean area 937,958), both several-fold above background |
| Other species | ≤18% detection rate each, background peak areas ~35,000–110,000 (e.g. *R. mucilaginosa* SUP, 74 samples, mean 84,750) |

**Interpretation.** This feature is present at low background levels across most of the genus but is both more prevalent **and** several-fold more abundant specifically in *R. taiwanensis*, in both cell pellet and supernatant — consistent with a true biosynthetic up-regulation in this species rather than simple presence/absence. The neutral mass (~480 Da) losing two waters readily to ionize is chemically diagnostic of a polyhydroxylated, ring-containing skeleton. Candidate classes, ranked:

1. **Oxidized/hydroxylated sterol or ergosterol-pathway intermediate** (medium–high plausibility) — facile loss of two H₂O from the protonated ion is classic behavior for di-/tri-hydroxylated sterols or triterpenoids (e.g. hydroxy-ergosterol-type intermediates). A formula near C₂₉–₃₀H₄₈O₃–₄ (MW ~460–480) fits, and *Rhodotorula* has a well-characterized, active ergosterol/carotenoid pathway — the strongest chemotype match for this genus.
2. **Oxylipin / long-chain hydroxy-fatty-acid derivative** (medium) — a di-hydroxy C28–30 fatty acid or wax-ester-type species could show similar dehydration, though less typical for yeast.
3. **Glycosphingolipid aglycone / ceramide fragment** (low–medium) — phytoceramides dehydrate readily but usually show a single dominant water loss and different RT/adduct behavior.
4. **Carotenoid-pathway apocarotenoid intermediate** (low) — *Rhodotorula* carotenoids (torulene, torularhodin) fall in a similar mass range (~400–540 Da) but are typically more conjugated/nonpolar than this feature's moderate RT suggests.

## 4. What would be needed to confirm either identity

- The actual MS2 fragment ion list (m/z + intensity) for diagnostic neutral losses (sugar oxonium ions, sterol ring fragments, fatty-acyl losses)
- Sub-5-ppm accurate mass with formula/RDBE prediction
- Isotope-envelope confirmation of charge state (especially important for Feature 1's 2+ assignment)
- Independent biological replicates to rule out batch/contaminant artifacts (particularly Feature 1, n=6 detections)
- Spectral library match (GNPS/MassBank) or authentic standard comparison

All class assignments above should be treated as **prioritization hypotheses**, not confirmed identifications.

## 5. Suggested next steps

1. Extract the raw MS2 fragment spectra for row IDs 20534 and 47208 from the source `.mzML`/MGF files (`SUP_302.mzML` scan 1601 and `SUP_155.mzML` scan 2999, respectively) and run GNPS/SIRIUS/CSI:FingerID formula and spectral-library matching.
2. For Feature 1, verify the 2+ charge assignment via isotope spacing (0.5 Th between isotope peaks) and check for a monomer-mass match in the same run (in-source dimer check).
3. For Feature 2, compare against ergosterol-pathway intermediate standards (e.g. ergosterol, ergosterol peroxide, hydroxy-ergosterol species) if available, given the strong genus-level prior for this pathway.
4. Re-run detection on raw (unfiltered) data specifically split by sample type (C vs SUP) per the existing `SUMMARY.md` "Suggested Next Steps," since both features are secretion-associated.
