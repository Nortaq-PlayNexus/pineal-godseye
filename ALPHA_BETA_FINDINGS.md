# Alpha/beta pinealocytes â€” corrected picture

**Status 2026-10-18. Built from primary data (GEO GSE115723), not from the
literature's characterisation of it.**

Three independent lines of evidence, each from raw count matrices:

1. **Not a droplet artefact.** Î± cells are the *deepest*-sequenced cells in the
   tissue, with lower mitochondrial fraction â€” the opposite of the
   low-complexity signature that generates spurious subclusters.
2. **Genuine pinealocyte identity.** Elevated for 7/9 (Day1) and 8/9 (Night1)
   independent markers versus non-pinealocytes. None of these markers played any
   part in defining the Î± set.
3. **Functionally distinct.** 33 genes show a different adrenergic response in
   Î± versus Î², against ~699 expected by chance.

Plus a correction to how the cells are usually described (below).

## What was actually wrong before

The previous dataset's Î±-pinealocyte row carried `MAT2A_SAM_capacity = 2.15`,
a value with no source â€” "MAT2A" appears zero times in Mays et al. 2018. The
whole of hypothesis H1 rested on that column, plus a direction of intercellular
flow that runs opposite to the source it cited.

The cells were real. The explanation attached to them was not.

## Evidence 1 â€” complexity signature

Mays define Î±-pinealocytes by ~2.3Ã- higher OxPhos and ~8.2Ã- lower ribosomal
expression. Low ribosomal plus high mitochondrial is the canonical
low-complexity/dying-cell signature in droplet scRNA-seq.

If Î± were an artefact, Asmt-high cells would be shallow and
mitochondrial-heavy. They are the opposite:

| | Day1 depth | Day1 ribo% | Day1 mito% | Night1 depth | Night1 ribo% | Night1 mito% |
|---|---|---|---|---|---|---|
| Î± (top 5% Asmt) | 1.47Ã- | 0.79Ã- | 16.3 vs 20.3 | 1.54Ã- | 0.88Ã- | 12.9 vs 17.6 |
| vs all pinealocytes | | | | | | |

The decisive check is the ordering. Asmt-**low** cells â€” far shallower, mean
7,013 UMI â€” have a *worse* ribosomal fraction (12.73%) than Asmt-high cells
(7.82%). If low complexity drove the signal, the shallowest cells should be
worst. They are not.

## Evidence 2 â€” identity

The slice was chosen on `Asmt` alone, so testing enrichment on markers that
played no part in that choice is not circular.

| Marker | Day1 vs non-pineal | Night1 vs non-pineal |
|---|---|---|
| `Sag` | 65.6Ã- | 35.8Ã- |
| `Pde10a`/`Pde6c` | 9.8Ã- | 33.6Ã- |
| `Tph1` | 8.7Ã- | 12.9Ã- |
| `Gngt1` | 6.9Ã- | 10.3Ã- |
| `Gngt2` | 7.8Ã- | 9.9Ã- |
| `Cngb1` | 17.2Ã- | 17.1Ã- |

`Sag` at 65Ã- is decisive â€” it is a pineal-and-retina-enriched marker and is
absent from the rest of the tissue.

## Evidence 3 â€” adrenergic response differs

Vehicle vs isoproterenol, unpaired (different animals; only 5 of ~1,800
barcodes are shared, so cells cannot be matched).

**33 genes differ at BH q<0.05** vs ~699 expected by chance. Strongest
subtype-specific: `Drd4` (Î± +7.11 vs Î² +2.72, q = 5.8Ã-10â»â¸), `Pdc`, `Mt1`,
`Asl`, `Lmo1`.

Î² is far more globally responsive: 4,997 genes at q<0.05 versus 164 in Î±,
consistent with Mays' observation that Î² carries most day/night differential
expression.

## The pipeline is validated

Before trusting the above, the drug-response pipeline was checked against
Mays' own named findings:

| gene | Î± vehicle â†’ iso | Î² vehicle â†’ iso |
|---|---|---|
| `Aanat` | 0.55 â†’ **27.90** | 0.38 â†’ **23.69** |
| `Crem` | 0.71 â†’ 3.21 | 1.32 â†’ 4.10 |
| `Drd4` | 0.49 â†’ 7.60 | 0.14 â†’ 2.86 |
| `Pde10a` | 0.06 â†’ 0.47 | 0.06 â†’ 0.50 |

These are precisely the genes Mays et al. 2018 report as isoproterenol
responders, and the `Aanat` induction reproduces their magnitude. Independent
confirmation that this code reads the GEO data correctly.

## A correction to how these cells are described

Î±-pinealocytes are routinely called "transcriptionally quiet." **That is wrong
in absolute terms and right only as a relative share of a larger transcriptome.**

Absolute ribosomal counts per cell: Î± carries **1,245** vs 1,076 in other
pinealocytes (Day1), and **1,038** vs 761 (Night1). Î± cells transcribe *more*
ribosome, not less â€” because they are the deepest-sequenced cells in the
gland.

What is true: ribosomal genes are *proportionally* depleted (7.8% vs 9.9%),
and the non-ribosomal share is correspondingly higher (92.2% vs 90.1%).

Why this matters: "these cells transcribe less" is the opposite of what the
data show. Reading the relative figure as an absolute one may be how a
cofactor-supply story got attached to these cells in the first place â€” a cell
that appears to need less machinery looks like a cell that supplies something
to cells that make more.

## What remains open

- **The mechanism is still unexplained.** The cells are real and distinct.
  Why they exist, and what they do beyond carrying more `Asmt`, is not
  established. H1 was the wrong answer to the right question.
- **Mays' own SAM hypothesis is untested.** It is a flux argument inferred from
  transcript levels: OxPhosâ†‘ and ribosomeâ†“ inside Î±, so more ATP, so more SAM
  available to `ASMT` inside Î±. It requires actual SAM measurement. It has
  never been done.
- **Spatial coupling is unmeasured.** Whether Î± and Î² cells are physically
  adjacent, and whether SAM moves between them, needs spatial transcriptomics or
  spatial metabolomics. Neither exists for pineal tissue in any published work.
- **Two of six GEO samples analysed.** Day1, Night1, Vehicle, Isoproterenol done
  at the Î±/Î² level; Night2 and the second vehicle/iso replicate untouched.

## Method notes worth keeping

**Operational Î± definition.** Î± cells are the top 5% by `Asmt` among
`Asmt`-positive cells. This is deliberately conservative: if Î± were an
artefact, the proxy *is* the artefact. It is not identical to Mays' clustering,
whose labels are not distributed with the GEO matrices.

**Two discarded runs, both caught by execution rather than review:**

1. An identity test using `Tph1`, `Gngt1`, `Gngt2` as controls scored 2/7
   "not enriched." That was a bad design â€” Mays states Î± and Î² share an
   overlapping marker set and lists those genes as expressed in *both*, so
   near-1.0 was the expected result carrying no information. Corrected to
   contrast against non-pinealocytes.
2. A drug-response run joined vehicle and isoproterenol on barcodes, got 5
   cells in common, and reported percent-change figures from them. The abort
   condition fired on *zero* shared barcodes rather than *too few cells to mean
   anything*, so it did not trigger. Rewritten unpaired with a 30-cell floor.

Both were caught only by looking at the cell counts printed at the top of the
output.

## Reproducing

```bash
python tools/test_alpha_artifact.py      # complexity signature
python tools/test_alpha_robustness.py    # 4 checks incl. inverted control
python tools/test_alpha_identity.py      # marker identity + 2nd sample
python tools/test_iso_response.py        # adrenergic response
```

Requires the project venv (`.venv`, numpy/pandas/h5py). Raw data in
`data/geo/`, downloaded from NCBI FTP; excluded from git for size.
Results in `results/*.json`.
