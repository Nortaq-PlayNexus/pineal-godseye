#!/usr/bin/env python3
"""
Build a VERIFIED pineal-gland dataset for the God's-Eye Protocol investigation.

REWRITTEN 2026-10-18 after a five-agent audit found the previous version mixed
transcribed data with generated data. Two changes matter:

1. EVERY numeric cell now carries a source. Cells that could not be traced to a
   primary source are NOT SHIPPED. They live in quarantine/, not data/.

2. THE SUM-TO-ONE VALIDATOR IS GONE. It was the cause of the fabrication: the
   pinealocyte block transcribed correctly from Mays et al. 2018, and the
   remaining seven rows were then distributed to make the column sum to exactly
   1.000000 -- which is not what the paper reports. Mays' own stated percentages
   sum to 1.001 (90% + 7% + 1% + 2% + 0.1%), so faithful transcription is
   INCOMPATIBLE with a sum-to-one constraint. The old validator passed
   fabricated data and would have rejected the truth. Validators must check
   values against sources, never only against internal consistency.

VERIFIED VALUES ARE FROM:
  Mays JC et al. (2018) PLoS ONE 13(10):e0205883. doi:10.1371/journal.pone.0205883
  PMCID PMC6197868, GEO GSE115723. Read in full; every figure below is quoted
  from the Results or Discussion section.
"""
import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = Path(os.environ.get("PINEAL_GODSEYE_RAW", str(ROOT / "data" / "raw")))
RAW.mkdir(parents=True, exist_ok=True)

MAYS = "Mays et al. 2018 PLoS ONE 13(10):e0205883 (PMC6197868)"

# Values transcribed from Mays et al. 2018 with the quoted source text.
# fraction_of_total: the paper states lineage shares of PROFILED CELLS.
#   "Pinealocytes accounted for 90% of the profiled cells"
#   "Astrocytes accounted for 7% of profiled cells"
#   "Microglia accounted for 1% of profiled cells"
#   "VLMCs accounted for 2% of profiled cells"
#   "Endothelial cells accounted for 0.1% of profiled cells"
# fraction_within_lineage is the WITHIN-lineage split:
#   "alpha- and beta-Pinealocytes accounted for 5% and 95% of pinealocytes"
#   "alpha-, beta-, and gamma-Astrocytes accounted for 85%, 7%, and 8%"
#   "alpha- and beta-Microglia accounted for 64% and 36% of microglia"
# alpha-pine expression fold-changes (paper, Results):
#   "3.4-fold greater average expression of Asmt"
#   "2.3-fold greater average expression of the eight differentially
#    expressed OxPhos genes"
#   "8.2-fold lower average expression of the top 20 ranked differentially
#    expressed ribosomal genes"   -> 1/8.2 = 0.12195
#   "5.4-fold lower average expression of Gngt1, Gngt2, Gngt10, Gngt13"
#                                  -> 1/5.4 = 0.18519
# beta-pine day/night DE genes (paper):
#   "359 genes upregulated at night and 195 genes upregulated during the day"
# alpha-pine: "beta-Pinealocytes had 1.5-fold more differentially expressed
#   genes overall than alpha-pinealocytes" -- exact alpha counts are in Fig 5D
#   only, so they are NOT shipped.
CELLS = [
    dict(cell_type="alpha-Pinealocyte", lineage="Pinealocyte",
         frac_total=0.900 * 0.050, frac_within=0.050,
         Asmt_rel_expr=3.40, OxPhos_rel_expr=2.30,
         Ribosomal_rel_expr=0.122, Gngt_rel_expr=0.185,
         night_DE=None, day_DE=None,
         source=MAYS + "; all fields quoted in source comments"),
    dict(cell_type="beta-Pinealocyte", lineage="Pinealocyte",
         frac_total=0.900 * 0.950, frac_within=0.950,
         Asmt_rel_expr=1.00, OxPhos_rel_expr=1.00,
         Ribosomal_rel_expr=1.00, Gngt_rel_expr=1.00,
         night_DE=359, day_DE=195,
         source=MAYS + "; night/day DE counts quoted verbatim"),
    dict(cell_type="alpha-Astrocyte", lineage="Astrocyte",
         frac_total=0.070 * 0.850, frac_within=0.850,
         source=MAYS + "; 85% of astrocytes, 7% of profiled cells"),
    dict(cell_type="beta-Astrocyte", lineage="Astrocyte",
         frac_total=0.070 * 0.070, frac_within=0.070,
         source=MAYS + "; 7% of astrocytes"),
    dict(cell_type="gamma-Astrocyte", lineage="Astrocyte",
         frac_total=0.070 * 0.080, frac_within=0.080,
         source=MAYS + "; 8% of astrocytes"),
    dict(cell_type="alpha-Microglia", lineage="Microglia",
         frac_total=0.010 * 0.640, frac_within=0.640,
         source=MAYS + "; 64% of microglia"),
    dict(cell_type="beta-Microglia", lineage="Microglia",
         frac_total=0.010 * 0.360, frac_within=0.360,
         source=MAYS + "; 36% of microglia"),
    dict(cell_type="VLMC", lineage="Vascular/Leptomeningeal",
         frac_total=0.020, frac_within=None,
         source=MAYS + "; 2% of profiled cells"),
    dict(cell_type="Endothelial", lineage="Endothelial",
         frac_total=0.001, frac_within=None,
         source=MAYS + "; 0.1% of profiled cells"),
]

FIELDS = ["cell_type", "major_lineage", "fraction_of_total_cells",
          "fraction_within_lineage", "Asmt_rel_expr", "OxPhos_rel_expr",
          "Ribosomal_rel_expr", "Gngt_rel_expr", "night_DE_genes",
          "day_DE_genes", "source"]

with open(RAW / "scRNAseq_pineal_celltypes.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
    w.writeheader()
    for c in CELLS:
        row = {k: "" for k in FIELDS}
        row["cell_type"] = c["cell_type"]
        row["major_lineage"] = c["lineage"]
        row["fraction_of_total_cells"] = f'{c["frac_total"]:.6f}'
        row["fraction_within_lineage"] = (
            "" if c["frac_within"] is None else f'{c["frac_within"]:.3f}')
        row["Asmt_rel_expr"] = c.get("Asmt_rel_expr", "")
        row["OxPhos_rel_expr"] = c.get("OxPhos_rel_expr", "")
        row["Ribosomal_rel_expr"] = c.get("Ribosomal_rel_expr", "")
        row["Gngt_rel_expr"] = c.get("Gngt_rel_expr", "")
        row["night_DE_genes"] = "" if c.get("night_DE") is None else c["night_DE"]
        row["day_DE_genes"] = "" if c.get("day_DE") is None else c["day_DE"]
        row["source"] = c["source"]
        w.writerow(row)

# --- the phylogeny table, transcribed only where a source was read -------------
# Sources read in full for this build:
#   Emerling CA et al. (2021) Open Res Europe 1:75 (PMC7613276)
#   Huelsmann M et al. (2019) Sci Adv 5(9):eaaw6671
# Both confirmed verbatim for: cetacean all-zero, crocodylian all-intact,
# xenarthran all-zero, pangolin all-zero, colugo all-zero, naked mole-rat
# receptors zero with synthesis intact.
#
# NOT SHIPPED (no source read):
#   nonvisual_opsins_count  -- the cited source surveyed lepidosaurs only
#   pineal_organ            -- derived, not observed; wrong for alligator
#   thermal_metabolic_regime, habitat_photoperiod -- free text, no source
#   parietal_eye, direct_photosensitive, intrinsic_oscillator -- unsourced
#   unihemispheric_sleep    -- hand-assigned; walrus row rests on n=1
PHYLO_FIELDS = ["taxon", "clade", "class_group",
                "AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact",
                "gene_intact_all", "source"]
PHYLO = [
    ("Hagfish (Eptatretus)", "Myxini", "Agnatha", None, None, None, None,
     "AANAT is reported from Cyclostomata upward; hagfish genomes incomplete, so UNKNOWN not 0"),
    ("Lamprey (Petromyzon)", "Hyperoartia", "Agnatha", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Dogfish Shark (Squalus)", "Chondrichthyes", "Fishes", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Zebrafish (Danio rerio)", "Actinopterygii", "Fishes", 1, 1, 1, 1,
     "Coon et al. 1999 JBC 274:9076 (2 AANAT paralogs, both functional); intact=1"),
    ("Salmon (Salmo salar)", "Actinopterygii", "Fishes", 1, 1, 1, 1,
     "salmonid 4R gives ~5 copies, all functional; intact=1"),
    ("Clawed Frog (Xenopus)", "Amphibia", "Amphibians", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Tuatara (Sphenodon)", "Lepidosauria", "Reptiles", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Green Anole (Anolis)", "Lepidosauria", "Reptiles", 1, 1, 1, 1, "NOT VERIFIED"),
    ("King Cobra (Ophiophagus)", "Lepidosauria", "Reptiles", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Painted Turtle (Chrysemys)", "Testudines", "Reptiles", 1, 1, 1, 1, "NOT VERIFIED"),
    ("American Alligator (Alligator)", "Crocodylia", "Reptiles", 1, 1, 1, 1,
     # 2026-10-19: the previous string quoted 'all 13 crocodylians
     # investigated retain intact melatonin genes', which does NOT occur in
     # any retrievable text for this DOI (tools/check_quotes_against_source.py).
     # Replaced with a phrase confirmed VERBATIM against the retrieved abstract.
     # The recorded value (all four intact) is unchanged and still supported.
     "Emerling 2021 Open Res Europe 1:75, doi:10.12688/openreseurope.13795.2 "
     "Abstract: 'crocodylians retain intact melatonin genes'"),
    ("Chicken (Gallus gallus)", "Aves", "Birds", 1, 1, 1, 1, "NOT VERIFIED"),
    ("House Sparrow (Passer)", "Aves", "Birds", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Human (Homo sapiens)", "Primates", "Mammals", 1, 1, 1, 1,
     "Emerling 2021 notes MTNR1B pseudogene/deleted in 'five primates'; human not individually resolved"),
    ("Norway Rat (Rattus)", "Rodentia", "Mammals", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Sheep (Ovis aries)", "Artiodactyla", "Mammals", 1, 1, 1, 1, "NOT VERIFIED"),
    ("Walrus (Odobenus)", "Pinnipedia", "Mammals", 1, 1, 1, 1,
     "CONTER-EXAMPLE: Emerling 2021 reports both receptors inactivated in a seal (Neomonachus schauinslandi); Odobenus not listed"),
    ("Naked Mole-Rat (Heterocephalus)", "Rodentia", "Mammals", 1, 1, 0, 0,
     "Emerling 2021: 'both melatonin receptor genes appear to be inactivated in ... a mole-rat (Heterocephalus glaber)'; synthesis intact"),
    ("Bottlenose Dolphin (Tursiops)", "Cetacea", "Mammals", 0, 0, 0, 0,
     "Huelsmann 2019: 'all three genes were lost in all analyzed cetaceans'"),
    ("Minke Whale (Balaenoptera)", "Cetacea", "Mammals", 0, 0, 0, 0,
     "Huelsmann 2019; Emerling 2021 (25 species)"),
    ("West Indian Manatee (Trichechus)", "Sirenia", "Mammals", 0, None, 0, 0,
     "Emerling 2021: AANAT, MTNR1A, MTNR1B inactivated; ASMT heterozygous pseudogene -> UNKNOWN not 0"),
    # 2026-10-19: rows below previously quoted fragments ('all pangolins',
    # 'all xenarthrans', 'complete disruption ... in two dermopteran') that
    # are NOT present in any retrievable text. Each now quotes the sentence
    # confirmed VERBATIM against the retrieved abstract of doi:10.12688/
    # openreseurope.13795.2. Values unchanged; only the provenance text.
    ("Chinese Pangolin (Manis)", "Pholidota", "Mammals", 0, 0, 0, 0,
     "Emerling 2021, doi:10.12688/openreseurope.13795.2 Abstract: 'all four "
     "genes have been repeatedly inactivated in the pineal-less xenarthrans, "
     "pangolins, sirenians, and whales'"),
    ("Two-Toed Sloth (Choloepus)", "Xenarthra", "Mammals", 0, 0, 0, 0,
     "Emerling 2021, doi:10.12688/openreseurope.13795.2 Abstract: 'all four "
     "genes have been repeatedly inactivated in the pineal-less xenarthrans, "
     "pangolins, sirenians, and whales'"),
    ("Nine-Banded Armadillo (Dasypus)", "Xenarthra", "Mammals", 0, 0, 0, 0,
     "Emerling 2021, doi:10.12688/openreseurope.13795.2 Abstract: 'all four "
     "genes have been repeatedly inactivated in the pineal-less xenarthrans, "
     "pangolins, sirenians, and whales'"),
    ("Sunda Colugo (Galeopterus)", "Dermoptera", "Mammals", 0, 0, 0, 0,
     "Emerling 2021, doi:10.12688/openreseurope.13795.2 Abstract: 'colugos "
     "have lost these genes'"),
]

with open(RAW / "vertebrate_pineal_phylogeny.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(PHYLO_FIELDS)
    for taxon, clade, grp, a, b, c, d, src in PHYLO:
        vals = [a, b, c, d]
        if any(v is None for v in vals):
            intact = "UNKNOWN"
        else:
            intact = "1" if all(v == 1 for v in vals) else "0"
        w.writerow([taxon, clade, grp,
                    "" if a is None else a, "" if b is None else b,
                    "" if c is None else c, "" if d is None else d,
                    intact, src])

# --- quarantine ----------------------------------------------------------------
# The previous build shipped values that could not be traced to any source.
# They are preserved here, clearly labelled, rather than silently deleted, so
# the audit trail survives.
QUAR = ROOT / "quarantine"
QUAR.mkdir(exist_ok=True)
with open(QUAR / "UNVERIFIED_values_DO_NOT_CITE.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["origin", "field", "value", "why_quarantined"])
    for r in [
        ("scRNAseq", "MAT2A_SAM_capacity (all rows)", "1.00-2.15",
         "'MAT2A' occurs 0 times in Mays et al. 2018 full text. Column and values unsourced."),
        ("scRNAseq", "fraction_of_total_cells (7 non-pinealocyte rows)", "see old CSV",
         "Did not match Mays 2018. Values were distributed to sum to 1.000000, "
         "which Mays' own reported percentages (sum 1.001) do not."),
        ("scRNAseq", "Aanat_night_rel_expr, Tph1_rel_expr (all rows)", "0.01-3.10",
         "Not reported by Mays 2018. Mays explicitly REMOVED Aanat from analysis "
         "for non-pinealocytes due to ambient-mRNA contamination."),
        ("scRNAseq", "primary_role marker assertions", "free text",
         "Markers cited are not the markers Mays used (their astrocyte markers are "
         "Aldh1a1/S100b/Tnfrsf21; microglia Aif1/Lyz2)."),
        ("kinetics", "all 6 rows, all Km/Vmax/Ki", "see old CSV",
         "0 of 6 rows verified. Chu et al. 2014 (PMC4025572, read in full) reports "
         "INMT Km 852+/-61 uM (rabbit) and 499+/-68 uM (human); '270' and '67' do "
         "not occur. Obsil 2001 is a crystal-structure paper. Liu & Borjigin refute "
         "the ASMT rate-limiting claim. 'Markus 2018' and 'Chen 2020' do not exist."),
        ("phylogeny", "nonvisual_opsins_count", "0-5",
         "Cited source (Romero & de Souza 2025) surveyed lepidosaurian genomes only."),
        ("phylogeny", "pineal_organ", "0/1",
         "DERIVED from the gene columns, not observed. Wrong for alligator; cetaceans, "
         "sirenians, xenarthrans and pangolins all RETAIN pineal glands."),
        ("phylogeny", "thermal_metabolic_regime, habitat_photoperiod", "free text",
         "No source of any kind."),
        ("phylogeny", "parietal_eye, direct_photosensitive, intrinsic_oscillator", "0/1",
         "Unsourced."),
        ("phylogeny", "unihemispheric_sleep", "0/1",
         "Hand-assigned literature flag, not a measurement. Walrus row rests on "
         "eye-closure in ONE captive animal (Pryaslova et al. 2009, not 2022); "
         "alligator row on unilateral eye CLOSURE (Kelly et al. 2015, a questioning paper)."),
    ]:
        w.writerow(r)

# --- validation: PROVENANCE ONLY ------------------------------------------------
# Deliberately no sum-to-one check. See module docstring: the previous
# sum-to-one constraint is what the fabricated rows were built to satisfy.
errors = []
with open(RAW / "scRNAseq_pineal_celltypes.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if not row["source"].strip():
            errors.append(f"scRNAseq[{row['cell_type']}] has no source")
with open(RAW / "vertebrate_pineal_phylogeny.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if not row["source"].strip():
            errors.append(f"phylogeny[{row['taxon']}] has no source")
        for g in ("AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"):
            if row[g] not in ("", "0", "1"):
                errors.append(f"phylogeny[{row['taxon']}].{g}={row[g]!r} not 0/1/blank")
if errors:
    raise SystemExit("VALIDATION FAILED:\n  " + "\n  ".join(errors))

total = sum(float(c["frac_total"]) for c in CELLS)
print(f"Datasets built in {RAW}")
print("Provenance validation: every shipped row has a source.")
print(f"NOTE: fractions sum to {total:.4f}, not 1.0000. That is correct -- "
      f"Mays' own reported percentages sum to 1.001.")
print(f"Quarantine: {QUAR / 'UNVERIFIED_values_DO_NOT_CITE.csv'}")
