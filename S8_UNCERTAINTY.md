# Section 8 — Uncertainty map

Classified per the protocol's categories. Every entry states what would resolve
it. Categories are the protocol's own; entries are specific to the pineal gland.

## Unknown because relevant observations have not been collected

| # | Question | Why unresolved | What would resolve it |
|---|---|---|---|
| U1 | What does pineal DMT do, functionally? | Only trace amounts measured; no perturbational study in humans | In vivo loss-of-function (INMT knockdown) with behavioural readout in a model organism that produces it |
| U2 | Is pineal calcification benign or pathological? | Correlational imaging + post-mortem histology only | Longitudinal MRI with concurrent histology in the same individuals |
| U3 | Do pinealocytes do anything besides make melatonin? | scRNAseq is descriptive; no perturbation | Cell-type-specific ablation in vivo with phenotype capture |
| U4 | What selective pressure drove melatonin-gene loss in 4 mammal clades? | Genomic signature present; no fitness measurement | Quantitative-genetic or comparative-fitness work; test the "unihemispheric sleep enabled it" hypothesis directly |
| U5 | Does the human pineal contribute to seasonal affective disorder? | Human pineal is small, imaging resolution-limited | Purpose-built high-field human pineal imaging across seasons |

## Unknown because measurements are inadequate

| # | Question | Why | What would resolve it |
|---|---|---|---|
| U6 | Is `unihemispheric_sleep` real in walrus? | Flag rests on eye-closure in ONE captive animal; no EEG | EEG on multiple free-ranging pinnipeds |
| U7 | Does the alligator show true hemispheric asymmetry? | Kelly et al. 2015 demonstrate unilateral eye CLOSURE only | EEG + local-field-potential in crocodilians |
| U8 | Actual AANAT paralog usage in each teleost | Dataset records a binary flag, not paralog identity | Per-species aanat1/aanat2 expression quantification |

## Unknown because theories disagree

| # | Question | Disagreement | Resolvable? |
|---|---|---|---|
| U9 | Did melatonin loss enable unihemispheric sleep, or co-occur with it? | Huelsmann 2019 says "may have been a precondition"; our data cannot test it (p = 0.63, n = 25, no phylogenetic correction) | Not with this dataset. Needs many more taxa + phylogenetic regression |
| U10 | Is AANAT or ASMT rate-limiting in humans? | AANAT is rate-limiting in rodent; Liu & Borjigin 2006 argue ASMT constrains at peak human flux | Simultaneous enzyme-activity measurement across the night |
| U11 | Are pineal volume reductions in depression/autism causal or incidental? | Meta-analyses exist; direction of causation unresolved | Longitudinal imaging with symptom tracking |

## Unknown because computationally difficult

| # | Question | Difficulty |
|---|---|---|
| U12 | Does pineal-gene loss track photoperiod unpredictability across all vertebrates? | Would need a genome-scale survey of melatonin genes vs photoperiod regime in hundreds of species. 25 taxa is far too few. |

## Unknown because experiments are impractical or too expensive

| # | Question | Barrier |
|---|---|---|
| U13 | Does the human pineal have a function beyond melatonin that matters clinically? | Human pinealectomy is rare and not ethical to study prospectively |
| U14 | Why did 4 lineages independently lose the gland? | Comparative genomics can identify the loss, not the selection on it |

## Unknown because disciplines have not been integrated

| # | Question | Missing integration |
|---|---|---|
| U15 | Is pineal calcification a cause or consequence of age-related neural change? | Chronobiology + mineralogy + neuroimaging rarely combined |
| U16 | Does the CVO immune axis (microglia, VLMC, MHC-II) respond to photoperiod change? | Immunology + chronobiology; the scRNAseq data has the cell types but no functional assay |

## Unknown because the question has not been formulated precisely

| # | Question | Problem |
|---|---|---|
| U17 | "Is the pineal a third eye?" | Not a scientific question as posed. Splits into: (a) does it transduce light? (b) does it share homology with retina? (c) is there a residual photoreceptive capacity? These need separate tests. |
| U18 | "Does the pineal produce DMT?" | Answered (yes, trace). The live question is what it does there. Original phrasing invites a yes/no that is already settled. |

## Potentially unknowable with current methods

| # | Question | Why |
|---|---|---|
| U19 | Did the gland ever serve consciousness/soul functions in any lineage? | Not testable; no measurement bears on it |

## Apparently unknown but literature may address

| # | Question | Status |
|---|---|---|
| U20 | Do pineal organoids recapitulate melatonin rhythm? | Yale 2024 reports they DO replicate circadian rhythms. Needs direct verification against that paper before treating as open. |

## Cross-cutting: what this investigation cannot do

- **No live animal data.** Every mechanistic claim here is second-hand from
  published literature. Nothing was reproduced experimentally.
- **25 taxa is a convenience sample**, not a phylogenetically balanced one.
  It over-represents mammals and lacks most non-mammalian diversity.
- **Phylogenetic non-independence is ignored** in the contingency table. Taxa
  share ancestry; treating 25 rows as 25 independent observations is a known
  methodological error, and it inflates the effective sample problem on top of
  the small n.
- **Two dataset columns are hand-assigned literature flags**
  (`unihemispheric_sleep`, `habitat_photoperiod`), not measurements.
