# Section 9 — Hypothesis register (REVISED 2026-10-18)

> **H1 and H2 are RETIRED.** Both were withdrawn after independent adversarial
> review found their evidence did not support them. This file now records the
> retired hypotheses alongside the surviving one, because the failures are more
> informative than the hypotheses were. See `AUDIT_2026-10-18.md`.

## Summary

| # | Hypothesis | Status | Basis |
|---|---|---|---|
| H1 | alpha/beta pinealocytes = SAM-supply division | **RETIRED** | Premise is prior art; direction inverted vs source; central datum unsourced |
| H2 | Gene-loss order topology-constrained | **RETIRED** | Misreading of the cited paper, wrong in the direction that manufactures the claimed asymmetry |
| H3 | USWS / melatonin-loss independence | **UNTESTED** | Data quarantined; design underpowered even when it existed |
| H4 | Unsourced values cluster outside validation coverage | **PARTIALLY SUPPORTED** | See correction below — downgraded 2026-10-18 after adversarial testing of my own validator |
| H5 | **The alpha/beta division is a real, functionally distinct population** | **SUPPORTED** | Three independent lines from primary GEO data |

## H5 — SUPPORTED (added 2026-10-18)

This is the project's first positive result, and it was reached by attacking
the cell type rather than by building on it.

**Claim.** Alpha-pinealocytes are a genuine pinealocyte population that differs
from beta-pinealocytes in identity, library complexity, and drug response.

**Evidence, all from GEO GSE115723 primary count matrices:**

1. **Not a droplet artefact.** Alpha cells are the deepest-sequenced cells in
   the tissue (1.47× Day1, 1.54× Night1) with lower ribosomal *and* lower
   mitochondrial fraction. The artefact hypothesis predicts the opposite. The
   ordering check is decisive: Asmt-LOW cells, which are far shallower, have a
   worse ribosomal fraction (12.73%) than Asmt-high cells (7.82%).
2. **Genuine identity.** Elevated for 7/9 (Day1) and 8/9 (Night1) markers not
   used to define the set. `Sag` at 65.6× and 35.8× is decisive.
3. **Functionally distinct.** 33 genes differ in adrenergic response between
   alpha and beta at BH q<0.05, against ~699 expected by chance. Strongest:
   `Drd4` (alpha +7.11 vs beta +2.72, q = 5.8e-08).

**Pipeline validated against the source.** `Aanat` induction (alpha
0.55 -> 27.90 CP10K), `Crem`, `Drd4`, `Pde10a` all reproduce Mays et al.
2018's reported isoproterenol response, confirming this code reads the GEO
data correctly.

**Correction attached to this finding.** Alpha cells are routinely called
"transcriptionally quiet." They are not — they carry MORE absolute ribosomal
counts than any other pinealocyte (1,245 vs 1,076 Day1; 1,038 vs 761 Night1).
The description is true only as a relative share of a larger transcriptome.
Reading it as absolute inverts the biology and may be how a cofactor-supply
story got attached to these cells.

**Scope.** H5 says the division is real and functional. It says nothing about
the mechanism. H1 was the wrong answer to a real question, and the question is
still open. Full detail in `ALPHA_BETA_FINDINGS.md`.

---

## H1 — RETIRED

**Was:** alpha-pinealocytes supply S-adenosylmethionine (SAM) to beta-pinealocytes.

**Why it died, three independent reasons:**

1. **The premise is published prior art.** Mays et al. 2018 abstract, verbatim:
   *"α-pinealocytes have transcriptomic differences that likely enhance melatonin
   formation by increasing the availability of the Asmt cofactor
   S-adenosylmethionine, resulting from increased production of a precursor of
   S-adenosylmethionine, ATP."*

2. **The distinguishing claim reverses the source.** Mays' Fig 7C moves
   **NAS beta → alpha**: *"NAS that is not converted to melatonin in
   β-pinealocytes enters the α-pinealocyte by passive diffusion through
   membranes and gap junctions."* H1 claimed SAM flows alpha → beta. The
   project's own dataset annotation said the same as Mays, and the H1 test
   never read it.

3. **The central datum does not exist.** `MAT2A_SAM_capacity = 2.15`. The
   string "MAT2A" occurs **zero times** in Mays et al. 2018 (verified against
   PMC6197868 full text). Mays' actual mechanism is a flux argument —
   OxPhos-derived ATP → SAM, *within* alpha — with no gene named.

Aggravating: the test grouped `MAT2A` (makes SAM) with `ASMT` (consumes SAM) as
a single "SAM-supply axis." ASMT enrichment marks alpha as the SAM **sink**,
which argues *against* H1. The test could not have failed on the transcriptional
axis, and printed "PREDICTION UPHELD" in every draft.

**Deeper problem nobody had raised until the audit:** alpha cells are 19× rarer
than beta but are credited with 2.15× the SAM capacity. A 5%-abundance
population cannot provision a 19×-larger consumer at equal flux. H1's
stoichiometry does not close.

**Also unresolved:** alpha-pinealocytes are defined by low ribosomal content
(8.2×) and high mitochondrial content (2.3×) — the classic low-quality/dying-cell
signature in droplet scRNA-seq. Mays profiled only n=275 alpha cells and did not
report per-cluster UMI distributions. **Whether alpha-pinealocytes are a real
cell type or a technical artefact is untested by anyone.** GEO GSE115723 permits
this check; nobody has run it.

---

## H2 — RETIRED

**Was:** receptor genes are lost before synthesis genes in melatonin-pathway
losing lineages, consistent with selection acting first on the signalling side.

**Why it died:** Emerling et al. 2021's actual stem-Cetacea pattern pairs a
**synthesis** gene with a **receptor**:

| gene | Emerling 2021 evidence | branch |
|---|---|---|
| AANAT | premature stop exon 1, all 25 species; dN/dS stem Cetacea **3.2** | **stem** |
| MTNR1A | "a whole exon deletion occurred on the stem Cetacea branch" | **stem** |
| ASMT | 1-bp insn exon 1, 21 spp., mysticetes + odontocetes | post-split |
| MTNR1B | 1-bp del exon 2, 17 spp., mysticetes + odontocetes | post-split |

There is **no receptor-before-synthesis asymmetry in cetaceans.** H2 was wrong
in precisely the direction that manufactures the pattern it claimed to find.

The Sirenia "corroboration" was also absent from the source: the only *shared*
inactivating mutation in sirenians is MTNR1B — the gene H2 called post-split in
cetaceans. AANAT and MTNR1A inactivations are species-specific (manatee,
dugong).

Emerling's whole framing is relaxed selection toward omega ≈ 1 — drift, which is
the alternative H2 positioned against, and which their dN/dS data supports.

**The deeper error:** H2 asserted exactly the inference S8 U14 had already
declared unobtainable (*"Comparative genomics can identify the loss, not the
selection on it"*). S8 and S9 contradicted each other and nothing reconciled
them. Emerling's Supplementary Tables S2–S5 contain per-species calls that could
actually test H2. They were never opened.

---

## H3 — UNTESTED (not refuted)

**Was:** the association between melatonin-pathway loss and unihemispheric sleep
is absent.

**Now:** the `unihemispheric_sleep` column is **quarantined**, so the
association cannot be computed from project data at all.

Two independent reasons it was never going to work:

1. **Power.** n=25, 7 events, actual test size 0.0267 (not 0.05). Minimum
   detectable odds ratio at 80% power is **~27×**. A 2-fold effect is
   undetectable in either direction.
2. **The column mixes mammal and non-mammal architectures.** Melatonin-**duration**
   reading is a mammalian readout rule that birds demonstrably do not use —
   Juss et al. 1993, *Proc R Soc B* 254:21-28: *"unlike mammals, birds do not
   read the duration of the melatonin signal to measure scotoperiod."* Averaging
   chicken and hamster rows is averaging two different systems.

**Correct statement:** the hypothesis is OPEN. The dataset could not test it in
either direction. Huelsmann 2019 is neither supported nor refuted here.

---

## H4 — SURVIVED, and is now the project's main finding

**Claim:** unsourced values in this project cluster systematically *outside*
the coverage of any validation rule, and inside it wherever a constraint could
be satisfied by fitting rather than by transcription.

**Evidence:**

| Finding | Location |
|---|---|
| Fractions fabricated to sum to exactly 1.000000 | scRNAseq, 7 rows |
| `MAT2A` column and value unsourced | scRNAseq |
| All 6 kinetics rows untraceable; 2 cited papers don't exist | kinetics |
| 5 morphology columns derived from gene columns, not observed | phylogeny |
| 2 free-text ecological columns with no source | phylogeny |
| `pineal_organ=0` for alligator — flatly wrong | phylogeny |

The pinealocyte block was transcribed correctly. Everything that could be
*satisfied by constraint* was generated. The sum-to-one validator added on
2026-10-17 passes every fabricated row. **The validator validated the
fabrication.**

**Prediction (untested):** if this pattern recurs, other projects that added
schema validators to internally-consistent data will show the same failure. This
project has no cross-project data, so the generalisation is unproven — and
noting it as proven is the error this project has already made once.

**Test that would settle it:** re-analyse a sibling project's data against its
sources, checking specifically whether its validated columns are also its
correct ones.

### CORRECTION 2026-10-18 — H4 downgraded after adversarial testing

The cross-project test (commit `6a72ff7`) and an adversarial test of this
project's own validator (`tools/test_provenance_adversarial.py`) together
downgrade H4. Full audit: `AUDIT_H4_VERIFICATION_2026-10-18.md`.

**What I originally claimed, and why it was too strong:**

1. *"Provenance architecture explains why the projects differ in their
   resistance to fabricated values."* — The Phase 3 result refutes this for
   this project. Every shipped row now carries a source string, and the
   validator still cannot detect:
   - a fabricated but well-formed citation
   - a citation to a nonexistent file
   - a value that contradicts its source (test 5c: set crocodylian
     `AANAT_intact=0` against Emerling 2021 and recompute the derived column —
     it passes)

   The validator is **syntactic**. It checks that a citation is present and
   non-blank. It never resolves a citation and never compares a value to its
   source.

2. *"The pattern does not recur in the comparison projects."* — Coverage was
   2 catalogued files of 22,679 under `data/`, and the high-risk category
   (hand-entered reference constants) was not audited. Two problems were found
   in the comparison project itself: `verified=True` is a hardcoded literal
   (never computed), and `bytes_expected`/`bytes_actual` are tautological on
   the genome path.

**Corrected status:**

- **H4a** (pattern does not recur) — **INCONCLUSIVE**. "No instance found in
  the material inspected" is all that was established.
- **H4b** (architecture explains the difference) — **PARTIALLY SUPPORTED**.
  The observable part holds: a value originating as a digest of real bytes has
  a receipt that fitting cannot reproduce; a hand-typed value does not. The
  causal claim does not: three projects differ in four ways at once, and
  recording provenance is demonstrably not the same as verifying it.

**What survives, and is now the project's transferable finding:** a validator
that tests internal consistency cannot distinguish a transcribed value from a
well-fitted fabrication, because satisfying the constraint is exactly how the
fabrication was built. That claim was tested directly (case 5c) and holds.

**What does not survive:** any suggestion that adding provenance fields gives a
project resistance to fabricated data. It gives traceability. Verification
requires resolving each citation and comparing the value to the source — a
step this build does not perform.

**Update, 2026-10-18 — partial support.** The audit did what it could to test
this on the primary data, and the pattern held. Of the GSE115723 alpha-pinealocyte
parameters: `Asmt` 3.4x, `OxPhos` 2.3x, ribosomal 1/8.2 = 0.122 and
`Gngt1` 1/5.4 = 0.185 are all **verified against Mays 2018**, because they were
read off the paper. `MAT2A_SAM_capacity` was **not**, because it was not. Both
categories sit in the same file, in the same column block, with the same
validator. Verified values and fabricated values were indistinguishable to every
check the project had at the time.

The distinguishing feature is now clear: **a value read from a source has a
receipt; a value generated to fit a constraint does not.** The fix is per-cell
provenance, not more validators. That is now implemented — every shipped row
carries a source string, and the build refuses rows without one.

**Still untested:** whether this pattern recurs in a sibling project. Four
agents found defects here that four prior audit passes missed; nobody has
looked at another project. That remains the open half of H4.
