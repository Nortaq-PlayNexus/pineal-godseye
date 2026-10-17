# pineal-godseye

The God's-Eye Protocol investigation of the pineal gland.

## Read these first

| File | What it is |
|---|---|
| AUDIT_2026-10-18.md | Five-agent audit. What was wrong, what was fixed. |
| ALPHA_BETA_FINDINGS.md | The project's positive result, from primary data |
| S9_HYPOTHESES.md | Hypothesis register. H1/H2 retired, H5 supported |
| S8_UNCERTAINTY.md | Open questions. Contains known-false claims |
| KNOWLEDGE_LANDSCAPE.md | Knowledge map. Contains known-false claims |

## Current state

The dataset was rebuilt from verified sources on 2026-10-18 after an audit
found fabricated values in it. It now ships 25 gene-intactness rows and 9
cell-type rows, every one carrying a source.

Three positive findings, all from GEO GSE115723 primary count matrices:

1. Alpha-pinealocytes are **not** a droplet artefact (they are the
   deepest-sequenced cells, not the shallowest)
2. They have **genuine pinealocyte identity** (7/9 and 8/9 independent markers)
3. Alpha and beta respond **differently** to adrenergic stimulation
   (33 genes at BH q<0.05 vs ~699 expected)

Two hypotheses were retired. H1 built a mechanism on a column that does not
exist in the source it cites. H2 misread a paper in the direction that
manufactured the pattern it claimed to find. The mechanism of the alpha/beta
division remains **unexplained**.

## Setup

`ash
python -m venv .venv
.venv/Scripts/python.exe -m pip install numpy pandas h5py
`

GEO data downloads into data/geo/ (~180 MB, excluded from git):

`ash
python tools/geo_load.py
`

## Reproducing

`ash
python tools/build_datasets.py           # rebuild + provenance check
python tools/verify_fractions.py         # fractions vs Mays 2018
python tools/regression_test.py          # validator fires on both defect classes
python tools/test_fisher_exact.py        # 5 derivable invariants
python tools/test_alpha_artifact.py      # complexity signature
python tools/test_alpha_robustness.py    # inverted control + sensitivity
python tools/test_alpha_identity.py      # marker identity, 2 samples
python tools/test_iso_response.py        # adrenergic response
`

## Standing rules

1. A validator that tests internal consistency cannot distinguish a
   transcribed number from a well-fitted fabrication. It must test against
   the source.
2. Never hardcode a summary table beside the claim it supports.
3. A schema that cannot express "unknown" forces unknown values into false
   ones. Add the third state before the data needs it.
4. Every self-correction in this project landed on code. None changed a
   scientific claim. That asymmetry is the lesson worth carrying.
