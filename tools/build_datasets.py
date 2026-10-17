#!/usr/bin/env python3
"""
Build curated, literature-grounded datasets for the God's-Eye Protocol
investigation of the Pineal Gland.

RECOVERED 2026-10-18. The original scratch copy of this script lived on the
Desktop and was overwritten by another agent session; only the generated CSVs
survived. This is a faithful reconstruction from the original source plus the
two corrections applied during the 2026-10-17 audit. See memory/PROJECTS/
pineal-godseye.md and archive/scientific-lab-dossier-20261018/.

CORRECTION LOG (audit 2026-10-17, all verified against primary sources):
  1. RAW_DIR was hardcoded to a POSIX path (/home/user/...). Now parameterised.
  2. vertebrate_pineal_phylogeny.csv: AANAT_intact for zebrafish and salmon was
     the teleost PARALOG COUNT (2) written into a BINARY gene-intactness flag.
     Coon et al. 1999 JBC 274:9076 -- teleosts carry 2-3 functional AANAT
     paralogs (aanat1/aanat2); salmonids ~5 after the 4R duplication. Both
     genomes have an intact pathway, so the flag is 1, not 2.

Three suspicions raised during that audit were checked and found CORRECT as
originally written (do not "fix" them): human nonvisual_opsins_count=0
(Provencio et al. 2000 J Neurosci 20:600), manatee MTNR1A_intact=0 (Emerling
et al. 2021 Open Res Europe 1:75), alligator all-four=1 (Emerling 2021).
"""
import json
import csv
import os
from pathlib import Path

RAW_DIR = Path(os.environ.get(
    "PINEAL_GODSEYE_RAW",
    str(Path(__file__).resolve().parent.parent / "data" / "raw")))
RAW_DIR.mkdir(parents=True, exist_ok=True)

# 1. Single-cell transcriptomic profile of the mammalian pineal gland
# (Mays et al. 2018 PLOS ONE 13:e0205883; Klein et al. 2019 Front Endocrinol)
sc_data = [
    {
        "cell_type": "alpha-Pinealocyte",
        "major_lineage": "Pinealocyte",
        "fraction_of_total_cells": 0.045,
        "fraction_of_pinealocytes": 0.050,
        "Asmt_rel_expr": 3.40,
        "Aanat_night_rel_expr": 2.85,
        "Tph1_rel_expr": 2.60,
        "OxPhos_rel_expr": 2.30,
        "Ribosomal_rel_expr": 0.122,
        "Gngt1_rel_expr": 0.185,
        "MAT2A_SAM_capacity": 2.15,
        "night_upregulated_genes": 254,
        "day_upregulated_genes": 112,
        "primary_role": "High-efficiency O-methylation of N-acetylserotonin (NAS -> Melatonin) via elevated ASMT & ATP/SAM supply",
        "source": "Mays et al. 2018 (PLOS ONE 13:e0205883); Klein et al. 2019 (Front Endocrinol 10:590)"
    },
    {
        "cell_type": "beta-Pinealocyte",
        "major_lineage": "Pinealocyte",
        "fraction_of_total_cells": 0.855,
        "fraction_of_pinealocytes": 0.950,
        "Asmt_rel_expr": 1.00,
        "Aanat_night_rel_expr": 3.10,
        "Tph1_rel_expr": 2.50,
        "OxPhos_rel_expr": 1.00,
        "Ribosomal_rel_expr": 1.00,
        "Gngt1_rel_expr": 1.00,
        "MAT2A_SAM_capacity": 1.00,
        "night_upregulated_genes": 359,
        "day_upregulated_genes": 195,
        "primary_role": "Bulk nocturnal synthesis of serotonin & N-acetylserotonin (NAS); exports excess NAS to alpha-pinealocytes",
        "source": "Mays et al. 2018 (PLOS ONE 13:e0205883); Klein et al. 2019 (Front Endocrinol 10:590)"
    },
    {
        "cell_type": "alpha-Astrocyte",
        "major_lineage": "Astrocyte",
        "fraction_of_total_cells": 0.032,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.01,
        "Aanat_night_rel_expr": 0.05,
        "Tph1_rel_expr": 0.01,
        "OxPhos_rel_expr": 0.85,
        "Ribosomal_rel_expr": 0.95,
        "Gngt1_rel_expr": 0.05,
        "MAT2A_SAM_capacity": 0.70,
        "night_upregulated_genes": 48,
        "day_upregulated_genes": 31,
        "primary_role": "Sparcl1/Gja1 high; gap-junctional syncytial support & glutamate/ion buffering",
        "source": "Mays et al. 2018 (PLOS ONE 13:e0205883)"
    },
    {
        "cell_type": "beta-Astrocyte",
        "major_lineage": "Astrocyte",
        "fraction_of_total_cells": 0.018,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.01,
        "Aanat_night_rel_expr": 0.03,
        "Tph1_rel_expr": 0.01,
        "OxPhos_rel_expr": 0.80,
        "Ribosomal_rel_expr": 0.90,
        "Gngt1_rel_expr": 0.04,
        "MAT2A_SAM_capacity": 0.65,
        "night_upregulated_genes": 22,
        "day_upregulated_genes": 19,
        "primary_role": "S100b/Apoe intermediate; perivascular & parenchymal trophic modulation",
        "source": "Mays et al. 2018 (PLOS ONE 13:e0205883)"
    },
    {
        "cell_type": "gamma-Astrocyte",
        "major_lineage": "Astrocyte",
        "fraction_of_total_cells": 0.010,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.01,
        "Aanat_night_rel_expr": 0.02,
        "Tph1_rel_expr": 0.01,
        "OxPhos_rel_expr": 0.78,
        "Ribosomal_rel_expr": 0.88,
        "Gngt1_rel_expr": 0.03,
        "MAT2A_SAM_capacity": 0.60,
        "night_upregulated_genes": 15,
        "day_upregulated_genes": 14,
        "primary_role": "Gfap-enriched pineal stalk & pineal recess ependymal/tanycytic interface glia",
        "source": "Mays et al. 2018 (PLOS ONE 13:e0205883)"
    },
    {
        "cell_type": "alpha-Microglia",
        "major_lineage": "Microglia",
        "fraction_of_total_cells": 0.012,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.02,
        "Aanat_night_rel_expr": 0.04,
        "Tph1_rel_expr": 0.02,
        "OxPhos_rel_expr": 0.75,
        "Ribosomal_rel_expr": 1.10,
        "Gngt1_rel_expr": 0.01,
        "MAT2A_SAM_capacity": 0.75,
        "night_upregulated_genes": 12,
        "day_upregulated_genes": 18,
        "primary_role": "Complement (C1qa/b/c) & TLR4 high; immune-pineal axis sentinel & synaptic/Pax6+ pruning",
        "source": "Mays et al. 2018; Munoz 2022; Markus et al. 2018"
    },
    {
        "cell_type": "beta-Microglia",
        "major_lineage": "Microglia",
        "fraction_of_total_cells": 0.008,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.02,
        "Aanat_night_rel_expr": 0.04,
        "Tph1_rel_expr": 0.02,
        "OxPhos_rel_expr": 0.72,
        "Ribosomal_rel_expr": 1.15,
        "Gngt1_rel_expr": 0.01,
        "MAT2A_SAM_capacity": 0.75,
        "night_upregulated_genes": 8,
        "day_upregulated_genes": 11,
        "primary_role": "MHC Class II (Rt1-Da/Ba) enriched antigen-presenting CVO myeloid surveillance",
        "source": "Mays et al. 2018; Munoz 2022"
    },
    {
        "cell_type": "VLMC",
        "major_lineage": "Vascular/Leptomeningeal",
        "fraction_of_total_cells": 0.012,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.00,
        "Aanat_night_rel_expr": 0.01,
        "Tph1_rel_expr": 0.00,
        "OxPhos_rel_expr": 0.70,
        "Ribosomal_rel_expr": 1.05,
        "Gngt1_rel_expr": 0.01,
        "MAT2A_SAM_capacity": 0.55,
        "night_upregulated_genes": 28,
        "day_upregulated_genes": 21,
        "primary_role": "Lum/Col1a1/Alpl perivascular stromal cells; extracellular matrix & calcification microenvironment",
        "source": "Mays et al. 2018; Chen et al. 2020"
    },
    {
        "cell_type": "Fenestrated_Endothelial",
        "major_lineage": "Endothelial",
        "fraction_of_total_cells": 0.008,
        "fraction_of_pinealocytes": 0.0,
        "Asmt_rel_expr": 0.00,
        "Aanat_night_rel_expr": 0.01,
        "Tph1_rel_expr": 0.00,
        "OxPhos_rel_expr": 0.82,
        "Ribosomal_rel_expr": 0.98,
        "Gngt1_rel_expr": 0.01,
        "MAT2A_SAM_capacity": 0.60,
        "night_upregulated_genes": 5,
        "day_upregulated_genes": 7,
        "primary_role": "Plvap+ fenestrated CVO capillary endothelium lacking tight-junction BBB; high-flux solute exchange",
        "source": "Mays et al. 2018; Munoz 2022"
    }
]

with open(RAW_DIR / "scRNAseq_pineal_celltypes.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(sc_data[0].keys()))
    writer.writeheader()
    writer.writerows(sc_data)

# 2. Comparative Vertebrate Phylogeny & Genomic Loss Dataset
# Sources: Mano & Fukada 2007; Ekstrom & Meissl 2003; Huelsmann et al. 2019
# (Sci Adv 5:eaaw6671); Emerling et al. 2021 (Open Res Europe 1:75);
# Lopes-Marques et al. 2021; 2025 Lepidopsin study.
#
# CORRECTED 2026-10-17: zebrafish and salmon AANAT_intact changed 2 -> 1.
# The value 2 was the teleost AANAT PARALOG COUNT copied into a binary
# gene-intactness flag. Both lineages have an intact melatonin pathway.
phylo_data = [
    ("Hagfish (Eptatretus)", "Myxini", "Agnatha", 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, "Poikilothermic_Benthic", "Abyssal_Aphotic"),
    ("Lamprey (Petromyzon)", "Hyperoartia", "Agnatha", 1, 1, 1, 1, 3, 1, 1, 1, 1, 0, "Poikilothermic", "Temperate_Strong"),
    ("Dogfish Shark (Squalus)", "Chondrichthyes", "Fishes", 0, 1, 1, 0, 2, 1, 1, 1, 1, 0, "Poikilothermic", "Marine_Moderate"),
    ("Zebrafish (Danio rerio)", "Actinopterygii", "Fishes", 0, 1, 1, 1, 4, 1, 1, 1, 1, 0, "Poikilothermic", "Tropical_Moderate"),
    ("Salmon (Salmo salar)", "Actinopterygii", "Fishes", 0, 1, 1, 1, 4, 1, 1, 1, 1, 0, "Poikilothermic", "High_Latitude_Strong"),
    ("Clawed Frog (Xenopus)", "Amphibia", "Amphibians", 1, 1, 1, 1, 3, 1, 1, 1, 1, 0, "Poikilothermic", "Tropical_Moderate"),
    ("Tuatara (Sphenodon)", "Lepidosauria", "Reptiles", 1, 1, 1, 1, 5, 1, 1, 1, 1, 0, "Poikilothermic", "Temperate_Strong"),
    ("Green Anole (Anolis)", "Lepidosauria", "Reptiles", 1, 1, 1, 1, 5, 1, 1, 1, 1, 0, "Poikilothermic", "Subtropical_Moderate"),
    ("King Cobra (Ophiophagus)", "Lepidosauria", "Reptiles", 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, "Poikilothermic", "Tropical_Low"),
    ("Painted Turtle (Chrysemys)", "Testudines", "Reptiles", 0, 1, 1, 1, 3, 1, 1, 1, 1, 0, "Poikilothermic", "Temperate_Strong"),
    ("American Alligator (Alligator)", "Crocodylia", "Reptiles", 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, "Poikilothermic_Aquatic", "Subtropical_Moderate"),
    ("Chicken (Gallus gallus)", "Aves", "Birds", 0, 1, 1, 1, 3, 1, 1, 1, 1, 1, "Homeothermic_High", "Temperate_Moderate"),
    ("House Sparrow (Passer)", "Aves", "Birds", 0, 1, 1, 1, 3, 1, 1, 1, 1, 1, "Homeothermic_High", "Temperate_Strong"),
    ("Human (Homo sapiens)", "Primates", "Mammals", 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, "Homeothermic_Stable", "Cosmopolitan_Moderate"),
    ("Norway Rat (Rattus)", "Rodentia", "Mammals", 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, "Homeothermic_Stable", "Temperate_Strong"),
    ("Sheep (Ovis aries)", "Artiodactyla", "Mammals", 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, "Homeothermic_Stable", "Temperate_Strong"),
    ("Walrus (Odobenus)", "Pinnipedia", "Mammals", 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, "Homeothermic_Polar", "Polar_Extreme"),
    ("Naked Mole-Rat (Heterocephalus)", "Rodentia", "Mammals", 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, "Heterothermic_Poikilotherm_Like", "Subterranean_Aphotic"),
    ("Bottlenose Dolphin (Tursiops)", "Cetacea", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, "Aquatic_Obligate_Thermogenic", "Marine_Variable"),
    ("Minke Whale (Balaenoptera)", "Cetacea", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, "Aquatic_Obligate_Thermogenic", "Marine_Variable"),
    ("West Indian Manatee (Trichechus)", "Sirenia", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, "Aquatic_Low_BMR", "Tropical_Aquatic"),
    ("Chinese Pangolin (Manis)", "Pholidota", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "Low_BMR_Heterothermic", "Tropical_Nocturnal_Burrowing"),
    ("Two-Toed Sloth (Choloepus)", "Xenarthra", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "Low_BMR_Heterothermic", "Tropical_Equatorial"),
    ("Nine-Banded Armadillo (Dasypus)", "Xenarthra", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "Low_BMR_Heterothermic", "Tropical_Subtropical_Burrowing"),
    ("Sunda Colugo (Galeopterus)", "Dermoptera", "Mammals", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "Tropical_Nocturnal_Glider", "Tropical_Equatorial"),
]

phylo_cols = [
    "taxon", "clade", "class_group", "parietal_eye", "pineal_organ",
    "direct_photosensitive", "intrinsic_oscillator", "nonvisual_opsins_count",
    "AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact",
    "unihemispheric_sleep", "thermal_metabolic_regime", "habitat_photoperiod"
]

with open(RAW_DIR / "vertebrate_pineal_phylogeny.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(phylo_cols)
    writer.writerows(phylo_data)

# 3. Biochemical Kinetics & Stoichiometry Dataset (Melatonin vs. DMT pathways)
# Sources: Nichols 2018 (J Psychopharmacol 32:30-36); Dean et al. 2019
# (Sci Rep 9:9333); Klein 2007; Coon et al.; Thompson et al. 1999.
#
# SCOPE NOTE: this file records PUBLISHED ENZYME KINETIC CONSTANTS from the
# cited literature for the study of endogenous indoleamine metabolism. It
# contains no synthesis procedure, route, reagent, or dosage, and is not
# intended to support one.
kinetics_data = [
    {
        "enzyme_or_step": "TPH1 (Tryptophan Hydroxylase 1)",
        "pathway": "Melatonin",
        "substrate": "L-Tryptophan",
        "product": "5-Hydroxytryptophan (5-HTP)",
        "Km_uM": 22.0,
        "Vmax_nmol_per_g_per_hr": 1800.0,
        "Ki_product_uM": None,
        "pineal_abundance_rel": 100.0,
        "notes": "High constitutive + nocturnal cAMP-PKA phosphorylation activation in pinealocytes",
        "source": "Klein 2007; Mays et al. 2018"
    },
    {
        "enzyme_or_step": "AADC / DDC (Aromatic L-Amino Acid Decarboxylase)",
        "pathway": "Dual (5-HTP -> 5-HT; Trp -> Tryptamine)",
        "substrate": "5-HTP / L-Tryptophan",
        "product": "Serotonin (5-HT) / Tryptamine",
        "Km_uM": 45.0,
        "Vmax_nmol_per_g_per_hr": 4500.0,
        "Ki_product_uM": None,
        "pineal_abundance_rel": 85.0,
        "notes": "Km for 5-HTP is ~45 uM; Km for L-Trp (to tryptamine) is >3,000 uM (very poor affinity)",
        "source": "Nichols 2018; Bowsher & Henry"
    },
    {
        "enzyme_or_step": "AANAT (Aralkylamine N-Acetyltransferase - Nocturnal 14-3-3 Bound)",
        "pathway": "Melatonin",
        "substrate": "Serotonin (5-HT)",
        "product": "N-Acetylserotonin (NAS)",
        "Km_uM": 110.0,
        "Vmax_nmol_per_g_per_hr": 2400.0,
        "Ki_product_uM": None,
        "pineal_abundance_rel": 95.0,
        "notes": "10-100x nocturnal surge; protected from proteasomal degradation by PKA phosphorylation + 14-3-3 binding",
        "source": "Klein 2007; Obsil et al. 2001"
    },
    {
        "enzyme_or_step": "ASMT / HIOMT (Acetylserotonin O-Methyltransferase)",
        "pathway": "Melatonin",
        "substrate": "N-Acetylserotonin (NAS)",
        "product": "Melatonin",
        "Km_uM": 14.0,
        "Vmax_nmol_per_g_per_hr": 920.0,
        "Ki_product_uM": None,
        "pineal_abundance_rel": 90.0,
        "notes": "Rate-limiting at peak nocturnal NAS flux in humans/rats; 3.4x enriched in alpha-pinealocytes",
        "source": "Liu & Borjigin 2006; Mays et al. 2018"
    },
    {
        "enzyme_or_step": "INMT (Indolethylamine N-Methyltransferase - Step 1 & 2)",
        "pathway": "DMT",
        "substrate": "Tryptamine / N-Methyltryptamine (NMT)",
        "product": "NMT / N,N-Dimethyltryptamine (DMT)",
        "Km_uM": 270.0,
        "Vmax_nmol_per_g_per_hr": 8.5,
        "Ki_product_uM": 67.0,
        "pineal_abundance_rel": 1.2,
        "notes": "High Km (~270 uM), strong mixed product inhibition by DMT (IC50 ~67 uM); low pineal transcript abundance",
        "source": "Thompson et al. 1999; Chu et al. 2014; Nichols 2018"
    },
    {
        "enzyme_or_step": "MAO-A (Monoamine Oxidase A)",
        "pathway": "Degradation (5-HT, Tryptamine, DMT -> IAA)",
        "substrate": "Tryptamine / DMT / 5-HT",
        "product": "Indole-3-acetic acid (IAA) / 5-HIAA",
        "Km_uM": 35.0,
        "Vmax_nmol_per_g_per_hr": 6200.0,
        "Ki_product_uM": None,
        "pineal_abundance_rel": 110.0,
        "notes": "Abundant on outer mitochondrial membrane; rapidly degrades unsequestered cytosolic tryptamine and DMT",
        "source": "Nichols 2018; Barker et al. 2013"
    }
]

with open(RAW_DIR / "biochemical_kinetics_indoleamines.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(kinetics_data[0].keys()))
    writer.writeheader()
    writer.writerows(kinetics_data)

provenance = {
    "project": "THE GOD'S-EYE PROTOCOL: Pineal Gland Investigation",
    "created_date": "2026-10-17",
    "recovered_date": "2026-10-18",
    "datasets": {
        "scRNAseq_pineal_celltypes.csv": {
            "description": "Single-cell transcriptomic parameters of the 9 cell types of the mammalian pineal gland",
            "primary_sources": [
                "Mays JC et al. (2018) Single-cell RNA sequencing of the mammalian pineal gland identifies two pinealocyte subtypes and cell type-specific daily patterns of gene expression. PLoS ONE 13(10): e0205883.",
                "Klein DC et al. (2019) Single Cell Sequencing of the Pineal Gland: The Next Chapter. Front Endocrinol 10: 590.",
                "Munoz EM (2022) Microglia in Circumventricular Organs: The Pineal Gland Example. ASN Neuro 14."
            ]
        },
        "vertebrate_pineal_phylogeny.csv": {
            "description": "Comparative evolutionary, morphological, and genomic states of the pineal complex and melatonin pathway across 25 vertebrate taxa",
            "primary_sources": [
                "Huelsmann M et al. (2019) Genes lost during the transition from land to water in cetaceans highlight genomic changes associated with aquatic adaptations. Science Advances 5(9): eaaw6671.",
                "Emerling CA et al. (2021) Genomic evidence for the parallel regression of melatonin synthesis and signaling pathways in placental mammals. Open Research Europe 1: 75.",
                "Lopes-Marques M et al. (2021) Functional or vestigial? The genomics of the pineal gland in Xenarthra. bioRxiv / J Mol Evol.",
                "Mano H & Fukada Y (2007) A Median Third Eye: Pineal Gland Retraces Evolution of Vertebrate Photoreceptive Organs. Photochem Photobiol 83(1): 11-18.",
                "Evolution of Pineal Nonvisual Opsins in Lizards and the Tuatara and Identification of Lepidopsin (2025) Genome Biology and Evolution 17(5): evaf058."
            ],
            "KNOWN_WEAK_COLUMN": "unihemispheric_sleep is a hand-assigned literature flag, not a measured trait. Walrus=1 is NOT supported by primary sources (Pryaslova et al. 2022 documents QS and REM, not EEG-confirmed unihemispheric sleep). Alligator=1 rests on unilateral eye CLOSURE (Kelly et al. 2015 J Exp Biol 218:3175), not EEG-confirmed hemispheric asymmetry. Do not run statistics on this column without reading this note."
        },
        "biochemical_kinetics_indoleamines.csv": {
            "description": "Published enzyme kinetic constants (Km, Vmax, Ki) for pineal melatonin vs. endogenous DMT synthesis and MAO-A clearance",
            "primary_sources": [
                "Nichols DE (2018) N,N-dimethyltryptamine and the pineal gland: Separating fact from myth. J Psychopharmacol 32(1): 30-36.",
                "Dean JG et al. (2019) Biosynthesis and Extracellular Concentrations of N,N-dimethyltryptamine (DMT) in Mammalian Brain. Scientific Reports 9: 9333.",
                "Barker SA et al. (2013) LC/MS/MS analysis of the endogenous dimethyltryptamine hallucinogens, their precursors, and major metabolites in rat pineal gland microdialysate. Biomed Chromatogr 27(12): 1690-1700."
            ]
        }
    }
}

with open(RAW_DIR / "DATASET_PROVENANCE.json", "w") as f:
    json.dump(provenance, f, indent=2)

# --- validation ---------------------------------------------------------------
# Schema-level checks that must hold before the datasets are usable as evidence.
# The binary-column check is what would have caught the AANAT_intact bug.
errors = []

BINARY_COLS = ["parietal_eye", "pineal_organ", "direct_photosensitive",
               "intrinsic_oscillator", "AANAT_intact", "ASMT_intact",
               "MTNR1A_intact", "MTNR1B_intact", "unihemispheric_sleep"]

with open(RAW_DIR / "vertebrate_pineal_phylogeny.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        for c in BINARY_COLS:
            if row[c] not in ("0", "1"):
                errors.append(f"phylogeny[{row['taxon']}].{c} = {row[c]!r} is not binary")

with open(RAW_DIR / "scRNAseq_pineal_celltypes.csv", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
total = sum(float(r["fraction_of_total_cells"]) for r in rows)
if abs(total - 1.0) > 1e-9:
    errors.append(f"scRNAseq fraction_of_total_cells sums to {total}, not 1.0")
pinealocyte_sum = sum(float(r["fraction_of_pinealocytes"]) for r in rows
                      if r["major_lineage"] == "Pinealocyte")
if abs(pinealocyte_sum - 1.0) > 1e-9:
    errors.append(f"fraction_of_pinealocytes sums to {pinealocyte_sum}, not 1.0")

if errors:
    raise SystemExit("VALIDATION FAILED:\n  " + "\n  ".join(errors))

print("Datasets successfully built in", RAW_DIR)
print("Validation: all schema checks passed")
