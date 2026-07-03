# Candidate biosynthetic gene families & the data needed to test them

Purpose: translate the **uniquely & highly secreted** MS2 features (see
`uniquely_secreted_features.csv`) into concrete gene families to screen for in the
*Rhodotorula* genome datasets, and specify exactly what data each test needs.

The logic is genotype↔phenotype: a metabolite that is **secreted** (enriched in
supernatant vs cell) and **species/strain-specific** should be explained by a
biosynthetic gene / gene cluster (and an exporter) that is **present in the producer
strains and absent/divergent in non-producers**. With 270 paired strains that have
both a genome link and paired metabolomics, this is directly testable as a
presence/absence (or SNP) association across strains.

---

## What the chemistry points to

Signals from the 130 uniquely-secreted features (all have MS2 spectra):

| Observation | Interpretation |
|---|---|
| Recurring neutral masses **~782, 787, 799, 1161, 1356 Da** seen as several adducts/charge states ([M+H]+, [M+2H]2+, [M+NH4]+) at **RT ~1–2.6 min** | A handful of real, moderately polar, mid-MW **lipid-like / glycolipid** compounds — the classic profile of secreted **polyol-ester / glycolipid biosurfactants** that *Rhodotorula* are known to excrete |
| A tail of **small polar masses (m/z 160–300, RT 4–5 min)** | Candidate secreted **siderophores / amino-acid-derived** metabolites (e.g. rhodotorulic acid, a hallmark *Rhodotorula* dihydroxamate siderophore, ~344 Da) and other peptide/organic-acid products |
| Strong **species structure** (44/130 top in *R. paludigena*, others in *sphaerocarpa, taiwanensis, dairenensis, pacifica*) | Species-restricted BGCs — the ideal case for presence/absence comparative genomics |

> These are hypotheses from mass/RT/adduct alone. **Structural annotation of the MS2
> spectra is the required next step** (see data list) before committing to a class.

---

## Candidate gene families to search the genomes for

Ranked by fit to the secreted, lipid-like, species-specific signal.

1. **Glycolipid / polyol-lipid biosurfactant clusters** (top priority for the ~780–800 Da RT~2 features)
   - Glycosyltransferases (GT-1 / GT family) — sugar–lipid linkage
   - Acyltransferases (BAHD / MBOAT / `Mac1`,`Mac2`-like) — fatty-acyl decoration
   - Fatty-acid synthase (FAS1/FAS2) and elongases/desaturases (`OLE1`) — the acyl backbone
   - Dedicated **MFS exporter** (`Mmf1`-like) — explains *secretion*
   - CYP450 fatty-acid hydroxylase (ω/subterminal) — adds hydroxyls seen in polyol lipids
   - *Reference model:* the Ustilaginomycete **MEL cluster** (`Emt1/Mac1/Mac2/Mmf1/Cyp1`) and *Rhodotorula* polyol-ester-of-fatty-acid (PEFA) literature.

2. **NRPS / NRPS-like — siderophores & secreted peptides**
   - NRPS (`SidD`-like) + L-ornithine-N5-monooxygenase (`SidA`-like) + N5-hydroxyornithine acetyltransferase → **rhodotorulic acid** pathway
   - Stand-alone NRPS-like (single-module) adenylation enzymes
   - Siderophore ABC/MFS transporters and esterases

3. **Terpenoid / carotenoid & apocarotenoid genes** (pigment axis; mostly cell-associated but oxidized apocarotenoids can be secreted)
   - `crtE` (GGPP synthase), `crtYB` (bifunctional phytoene synthase/lycopene cyclase), `crtI` (phytoene desaturase), `crtS`/CYP (torulene/torularhodin oxidation)
   - Carotenoid cleavage dioxygenases (CCD) → secreted apocarotenoids
   - (These tie the secretion analysis back to the existing color-phenotype GWAS.)

4. **Polyketide synthases (PKS)**
   - Type I iterative PKS and type III PKS + tailoring (KR/DH/ER domains, O-methyltransferases)

5. **Tailoring / diversifying enzymes that create species-specificity**
   - CYP450 monooxygenases, SAM-dependent methyltransferases, acetyl/acyltransferases, FAD/NAD(P) oxidoreductases, sulfotransferases, glycosyltransferases

6. **Export / secretion machinery (why it's in the supernatant at all)**
   - ABC transporters (PDR/MDR subfamilies), MFS transporters, lipid-transfer proteins

---

## Data needed to investigate candidates

You said you have genomes for most of the metabolite strains — here is what turns that
into a test, in priority order.

**A. To even map chemistry → genes (prerequisite)**
- **MS2 structural annotation** of the target features: run **SIRIUS + CSI:FingerID** (the
  `sirius-ms` env already exists in the lab) and/or **GNPS molecular networking** to get
  molecular formula + compound class per feature. Without this, gene mapping stays at the
  class-hypothesis level.
- **Raw data is on disk** at `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/`:
  - cell extracts: `mzML/C_*.mzML` (299 files, ~12 GB)
  - supernatant:   `ExFab_Supernatant/SUP_*.mzML` (304 files, ~13 GB)
  - MZmine feature-extraction workspace: `feature_extractMS2/` (source of the aligned table used here)
- **No MS2 export exists yet** (no MGF / GNPS / SIRIUS artifacts found). Next action is to
  export MS2 spectra for the target `feature_index` rows — either from the original MZmine
  project (GNPS-FBMN / SIRIUS export) or directly from the mzML — then feed SIRIUS/GNPS. The
  README notes a GNPS2 rerun by the Petras lab may also produce these.

**B. Genome assets (you have most of this)**
- **Assembled genomes** for the ~270 strains that have both a genome link and paired
  metabolomics (the metadata already carries `db_biosample_list` / `db_sra_run_list` /
  `db_bioproject_list` per strain — that is the strain↔genome key).
- **Structural + functional annotation** (proteomes + GFF): `funannotate` (lab already
  uses it) or equivalent; needed for orthology and BGC calls.
- A **strain ↔ genome-accession ↔ metabolite-abundance** join table (buildable now from
  the metadata + `secretion_scores_all_features.csv.gz`).

**C. Prediction & comparative-genomics layers**
- **antiSMASH (fungal mode)** on every genome → per-strain catalog of PKS/NRPS/terpene/
  other BGCs. This is the primary hypothesis generator for clusters 1–4 above.
- **OrthoFinder** (or similar) across the strains → gene-family presence/absence matrix.
- **Phylogenomic tree** (the `nf_phyling` skill / BUSCO markers) → needed to control for
  relatedness when associating gene presence with metabolite production.

**D. The association test (the actual experiment)**
- For each target secreted feature: correlate **per-strain supernatant abundance**
  (continuous, from this analysis) against **gene/ortholog/BGC presence-absence** (or SNPs)
  across the 270 strains — a **metabolite-GWAS / presence-absence association**, ideally
  phylogeny-corrected. The repo already contains GWAS scaffolding (`scripts/02_gwas_*`)
  that can be repurposed with the secreted-feature abundances as phenotypes.
- Strongest evidence = a BGC/ortholog whose presence tracks the producer species
  (e.g. *R. paludigena*-specific clusters for the 44 features topping in that species),
  and that includes both a synthase and an exporter.

**E. Confirmation (optional, later)**
- Knockout / heterologous expression of the top candidate cluster, or targeted LC-MS of a
  producer vs non-producer, to close the loop.

---

## Immediate, no-new-data next step

Everything in **A–D** except SIRIUS annotation can start now:
1. Build the strain↔genome↔metabolite join table from existing metadata.
2. Run antiSMASH + OrthoFinder on the genomes you already have.
3. Take the top ~10 uniquely-secreted features (this analysis) and test presence/absence
   association of BGCs/orthologs against their per-strain supernatant abundance.
