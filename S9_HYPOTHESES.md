# Section 9 — Hypothesis register

Each hypothesis is scored for novelty status, falsifiability, and current
support. **No hypothesis here is claimed as a discovery.** Most are explicitly
labelled as known-but-under-tested or as synthesis, per the protocol's rule
that novelty requires a targeted literature search, which I have not run for
any of these to completion.

---

## H1 — The alpha/beta pinealocyte division is a SAM-supply division, not a step-division

**Claim.** alpha-pinealocytes exist to supply S-adenosylmethionine (SAM) to
beta-pinealocytes, not to do a distinct enzymatic step. The 4.5%/95% split
inverts the enzyme abundance you'd expect if each cell type ran its own
complete pathway.

**Motivating observations.** From the dataset: alpha = 4.5% of cells,
5.0% of pinealocytes, but `Asmt_rel_expr` = 3.40 vs beta's 1.00, and
`MAT2A_SAM_capacity` = 2.15 vs 1.00. Meanwhile beta has `Ribosomal_rel_expr`
= 1.00 vs alpha's 0.122 — alpha cells are transcriptionally quiet, beta cells
are not. Mays et al. 2018 report ASMT is enriched in alpha cells.

**Reasoning.** If alpha cells were simply "more efficient at the last step,"
you would expect them to be transcriptionally active. They are the *opposite* —
0.122 ribosomal expression. A cell that is metabolically specialised but
transcriptionally quiet looks more like a cofactor supplier than a factory.

**Assumptions.** (a) `MAT2A_SAM_capacity` is a meaningful proxy for SAM supply;
(b) the `rel_expr` values are comparable across cell types; (c) transcriptional
quiescence implies a supply role.

**Alternative explanations.** (a) Alpha cells may simply be more differentiated
and therefore less transcriptionally active, with no SAM coupling at all.
(b) The 4.5% fraction may reflect cell-cycle state or dissociation artefact.
(c) The division may be an artefact of the clustering resolution in a single
scRNAseq paper.

**Predictions that would distinguish it.** If SAM coupling is real, then:
- perturbing MAT2A should selectively reduce melatonin output while leaving
  beta-cell AANAT activity unchanged;
- SAM and NAS should show *spatial* coupling (alpha-adjacent beta cells
  converting NAS preferentially) detectable by spatial transcriptomics.

**Evidence to support.** SAM microdialysis showing local depletion around
alpha cells; MAT2A knockdown reducing melatonin disproportionately to AANAT
knockdown at matched NAS levels.

**Evidence to falsify.** MAT2A knockdown with no effect on melatonin output;
or spatial transcriptomics showing alpha and beta cells are randomly
interleaved with no SAM gradient.

**Practical test.** Spatial transcriptomics on pineal tissue (the dataset has
no spatial dimension — scRNAseq discarded it). This is the single highest-value
next experiment.

**Novelty status.** Known but underexplored. The alpha/beta division is
established; the *functional interpretation* as SAM supply is a synthesis of
the dataset's own numbers, and I have NOT run a prior-art search for it.
Treat as hypothesis, not discovery.

---

## H2 — Gene loss order is constrained by pathway topology, not by drift

**Claim.** In the lineages that lost the melatonin pathway, the order of gene
loss is not random. Receptor genes (MTNR1B, MTNR1A) are lost *before* synthesis
genes (ASMT, AANAT) in every documented case, consistent with selection acting
first on the signalling side.

**Motivating observations.** Emerling 2021: in cetaceans, MTNR1B is inactivated
in the *stem* lineage, while ASMT and MTNR1A are inactivated *after* the
odontocete/mysticete split. Huelsmann 2019 reports the same asymmetry. The
dataset's phylogenetic pattern (all four genes lost together in 7 of 9 clades)
hides this ordering entirely, because it is binary.

**Reasoning.** Losing the receptor first disconnects the gland from its
downstream effects while leaving the (energetically costly) synthesis machinery
intact — a classic evolutionary sequence: stop paying for an output whose
consumers are gone, then stop paying for production.

**Assumptions.** (a) the inactivating mutations are genuinely ancestral, not
convergent; (b) absence of relaxed selection in Emerling's dN/dS table means
neutral loss rather than selective loss.

**Alternative explanations.** (a) Mutation rate differences between genes
could produce apparent ordering without selection. (b) The orderings may be
uncoincidentally similar and I have read a pattern into 4 lineages.

**Predictions.** If topology-constrained, receptor loss should also precede
synthesis loss in the *other* independently-lost clades (Xenarthra, Sirenia,
Pholidota) — not just cetaceans. If drift-driven, order should be uncorrelated
with pathway position.

**Current status.** **TESTABLE NOW AND ALREADY PARTLY ANSWERED.** Emerling
2021's own data is the test, and it reports the same stem-lineage-first pattern
in sirenians as in cetaceans. This is therefore NOT a novel hypothesis — it is
an existing empirical result I re-derived. Recording it as such is the point.

**Novelty status.** Established in the literature (Emerling 2021, Huelsmann
2019). Recorded here because it is the one place where the dataset's binary
encoding actively destroys real signal.

---

## H3 — `unihemispheric_sleep` and melatonin loss are independent

**Claim.** The association Huelsmann proposes does not hold in our data.
Melatonin-pathway loss and unihemispheric sleep co-occur in cetaceans and
sirenians, but unihemispheric sleep ALSO occurs in taxa with a fully intact
pathway (birds, possibly crocodilians, possibly pinnipeds), so loss is
neither necessary nor sufficient for the trait.

**Current support.** Odds ratio 0.46, **Fisher two-sided p = 0.63** — not
significant. n = 25. Phylogenetic independence ignored.

**Why this matters even as a null.** A null at n = 25 with 7 events cannot
distinguish "no association" from "association too small to detect." The
honest statement is: **this dataset has no power to test the hypothesis either
way.** That is itself the finding — the convenience sample is inadequate for the
question, which is why U12/U9 remain open.

**What would give it power.** Phylogenetic comparative analysis (e.g. Pagel's
lambda-corrected logistic regression) across hundreds of species with genome
data. Not feasible with 25 hand-selected taxa.

**Novelty status.** This is not a hypothesis — it is a null result, recorded as
such. It must not be presented as evidence against Huelsmann 2019.

---

## H4 — The dataset's binary encoding conceals more than it reveals

**Claim.** `vertebrate_pineal_phylogeny.csv` uses binary flags for what are
multi-state or continuous properties, and this has twice produced wrong
conclusions: (a) `AANAT_intact=2` (a paralog count in a binary flag), and
(b) `unihemispheric_sleep` conflating eye-closure observation with EEG-confirmed
hemispheric asymmetry.

**Why this is a hypothesis, not just a bug report.** If binary encoding is a
recurring failure mode, the prediction is that OTHER datasets in this estate
encode multi-state properties as binary and will fail the same way. That is a
testable, general claim with a defined search: audit the flag columns across
the project's datasets.

**Prediction.** At least one other dataset has a flag column mixing observation
types (e.g. "observed" vs "measured" vs "inferred from related species").

**Status.** Partly confirmed within this project (3 of 3 audit rounds found an
error). The cross-project extension is untested.

**Novelty status.** A synthesis of this project's own audit history. Not
claimed as novel.

---

## Register summary

| # | Hypothesis | Status | Novelty | Falsifiable now? |
|---|---|---|---|---|
| H1 | alpha/beta = SAM supply division | Untested | Known but underexplored | Yes — spatial transcriptomics |
| H2 | Loss order topology-constrained | **Already answered in literature** | Established | Already answered |
| H3 | Sleep/loss independence | **NULL, underpowered** | Not a hypothesis | No — needs more taxa |
| H4 | Binary encoding hides signal | Partly confirmed | Synthesis | Yes — cross-project audit |

**Nothing in this register is a discovery.** H2 is an existing result re-derived;
H3 is a null; H1 is the only genuinely open, testable hypothesis, and it needs
spatial data the project does not have.
