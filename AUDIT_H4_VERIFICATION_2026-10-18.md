# Independent Verification Audit — 2026-10-18

Audit of the H4 cross-project finding, commits `6a72ff7` / `1390268` / `838a3b4`,
and per-cell provenance enforcement.

**Verdict up front: H4b is downgraded to PARTIALLY SUPPORTED, and one of my own
claims is refuted.** My provenance validator is syntactic, not evidentiary.
Detail below.

---

## Phase 1 — Repository integrity

| Check | Result |
|---|---|
| `git status --short` | empty (clean) |
| Tracked files | 27 |
| Tracked under `.venv/` or `data/geo/` | 0 (correctly ignored) |
| `.gitmsg` tracked? | No |
| `.gitignore` unintended effects | None found |
| Test exit codes | all 0 |

Commands run, actual exit codes:

```
tools\build_datasets.py        EXIT 0
tools\verify_fractions.py      EXIT 0   (all lineage fractions match Mays 2018)
tools\regression_test.py       EXIT 0
tools\test_fisher_exact.py     EXIT 0
```

### FINDING 1 — a defect in my own commit sequence

Commit `6a72ff7` (the H4 finding) **also tracked `.gitmsg`**, a 35-line scratch
file containing the commit message itself. `1390268` removed it; `838a3b4` added
it to `.gitignore`.

Net state is clean, but two of the three commits I reported as the H4 finding
were partly cleanup of my own error. I reported "clean tree, everything
committed" without noticing that I had created the problem in the first place.
Recorded as a finding rather than smoothed over.

---

## Phase 2 — Independent verification of the H4 cross-project claims

### pineal-godseye: hand-typed literals — CONFIRMED

`tools/build_datasets.py` contains `Asmt_rel_expr=3.40`, `MAT2A`-adjacent
values and gene flags as Python literals. No digest, no retrieval. Confirmed by
reading the source.

### genome-anomaly-lab: digests computed from downloaded bytes — CONFIRMED, with a caveat

`gal/data/acquire.py`:
- L125-133 `_sha256(path)` — genuine streaming SHA-256
- L205 `digest = _sha256(dest)` (genome path)
- L292 `digest = _sha256(dest)` (ENA path)

Independently recomputed both catalogued files. Both match:

```
GCF_009914755.1.zip        sha256 MATCHES CATALOG
ERR1624145_1.fastq.gz     sha256 MATCHES CATALOG
```

### FINDING 2 — `verified=True` is a hardcoded literal, never computed

All 8 occurrences (`acquire.py` L211, L298; `targets.py` L63, L81, L110, L127,
L131, L135, L139, L143) set `verified=True` as a constant. Nothing in the
project computes a verification result and writes that flag. The column
**asserts** verification rather than **recording** it.

### FINDING 3 — `bytes_expected` and `bytes_actual` are tautological on the genome path

`acquire.py` L207-208:

```python
size = _download_resumable(url, dest)   # no expected= passed
rec = DownloadRecord(..., bytes_expected=size, bytes_actual=size, ...)
```

`_download_resumable`'s size-mismatch guard (L155-156) only fires when an
`expected` argument is supplied. The genome path supplies none, so
`bytes_expected` is set to the observed size. Both catalogued rows show
`bytes_expected == bytes_actual`, consistent with this.

On the ENA path `expected` comes from the ENA API (`fastq_bytes`, L277) and
*is* independent — so the check is meaningful there and tautological on the
genome path. The column pair looks like verification and is not, for that row
type.

### ScientificDiscoveryLab: hashes checked at verify time — CONFIRMED, strongest of the three

`04_SHARED_ENGINE/engine/utilities/core.py` L236-246:

```python
actual_artifact_hash = sha256_file(str(artifact))
report["artifact_hash_valid"] = actual_artifact_hash == artifact_info.get("sha256")
if not report["artifact_hash_valid"]:
    report["errors"].append("result artifact SHA-256 mismatch")
report["valid"] = (report["scope_valid"] and report["result_hash_valid"]
                   and report["artifact_hash_valid"] is True
                   and report["artifact_size_valid"] is True ...)
```

This **gates** the verdict — a hash mismatch makes the whole record invalid. It
is not decorative. The same pattern appears in
`experiment_verify.py` L128-135.

`tests/unit/test_result_provenance.py` contains both a positive assertion
(L48-53) and a **negative** test asserting `artifact_hash_valid is False`
(L116) — i.e. the failure path is tested.

### FINDING 4 — I could NOT execute SDL's tests

SDL has no project venv and `pytest` is not installed for the available
interpreters (`No module named pytest`, exit 1). The code-level evidence above
stands on inspection. **I am not claiming its tests pass.** That remains
unverified by execution.

---

## Phase 3 — Adversarial tests of my own provenance validator

`tools/test_provenance_adversarial.py`. Temporary environment; canonical data
untouched.

| # | Attack | Result | Which check fired |
|---|---|---|---|
| 1 | source field deleted | CAUGHT | source |
| 2 | whitespace-only source | CAUGHT | source |
| 3 | fabricated well-formed citation | **NOT CAUGHT** | — |
| 4 | citation to nonexistent file | **NOT CAUGHT** | — |
| 5 | value contradicts a real source | caught | **internal only** |
| 5c | same, derived column recomputed | **NOT CAUGHT** | — |
| 6 | value modified, source unchanged | caught | **internal only** |
| 7 | value AND source changed consistently | caught | **internal only** |

**caught by source check: 2 · internal only: 3 · not caught: 3**

### FINDING 5 — REFUTED: my provenance claim was too strong

I wrote, in commit `6a72ff7`:

> "require a per-cell source string and validate that every cell has one"

and in `S9_HYPOTHESES.md`:

> "every shipped row must carry a source string; the build refuses rows
> without one"

Those statements are literally true, but I presented them under the heading
"provenance enforcement", and I separately framed H4b as "provenance
architecture explains why the projects differ in their resistance to
fabricated values." **That framing implies resistance to fabrication. The
validator provides none.**

Case 5c is the decisive test: set `AANAT_intact=0` for the American alligator —
contradicting Emerling 2021, which says crocodylians retain all four genes — and
recompute the derived column so internal consistency holds. The validator
passes it.

The correct, defensible claim is narrow:

> **No shipped cell lacks a source string.**

That is a completeness check on citations, not verification of them. The
validator never resolves a citation and never compares a value to its source.

Note also that attacks 5, 6 and 7 were *initially* reported as caught by the
first version of this test, and I wrote an explanation saying they were "NOT
caught" — the explanation was wrong and the numbers right. They fire on
`gene_intact_all disagrees`, an internal-consistency check. **That is the same
class of constraint the original fabricated rows were built to satisfy.** The
test caught my own misattribution, which is the only reason it is correct now.

---

## Phase 4 — Calibrated verdicts

### H4a — "The pattern does not recur in the two comparison projects"

**Operational definition.** A hand-entered numeric value exists in the two
comparison projects that is not traceable to a source, and which satisfies an
internal-consistency constraint.

**Status: INCONCLUSIVE — cannot be established from what was tested.**

Evidence *for*: genome-anomaly-lab computes digests from bytes at download;
SDL gates verdicts on digest comparison. Neither showed the hand-typed-value
pattern in what I inspected.

Evidence *against / limits*:
- **Catalog coverage is 2 files of 22,679** under `data/`. "No fabrication
  found" over 2 catalogued files says little about the other 22,677.
- **`verified=True` is hardcoded** (Finding 2), so the catalog's own claim of
  verification is unearned.
- **I did not audit either project's hand-entered reference values**, which is
  the actual high-risk category. genome-anomaly-lab has 962 test files and hard
  constants (`MAX_DP_FRAGMENT`, effect-size expectations) I did not trace.
- **SDL's tests could not be executed** (Finding 4).
- Both projects have documented histories of *estimator* defects — which are a
  different failure mode from fabricated data.

**Falsifier.** A hand-typed constant in either project with no source, found on
inspection.

**What I should have said originally:** "I found no instance in the material I
inspected." Not "the pattern does not recur."

### H4b — "Provenance architecture explains the difference in resistance"

**Status: PARTIALLY SUPPORTED — and narrower than I claimed.**

**Supported:** for the three projects, the origin of a numeric value predicts
whether it has a receipt. Values originating as digests of real bytes carry a
receipt that fitting cannot reproduce. Values typed by hand do not. That
contrast is real and is visible in the code.

**Not supported:** the causal claim that architecture *explains* resistance.
Three projects cannot establish a general rule, and my own Phase 3 result
demonstrates the limit directly — pineal-godseye now has per-cell source
strings on every row and **still cannot detect a fabricated citation or a value
contradicting its source**. Architecture that records provenance is not
architecture that verifies it.

**Additional confound.** SDL and genome-anomaly-lab differ from pineal-godseye
in at least four ways at once: value origin, provenance architecture, test count,
and audit history. Attributing the difference to architecture alone is not
warranted.

**Falsifier.** A project with digest-based provenance that nonetheless ships a
fabricated value; or pineal-godseye-style citation records that successfully
detect a fabricated citation without source resolution.

---

## Findings summary

| # | Finding | Severity |
|---|---|---|
| 1 | My commit `6a72ff7` tracked a scratch file; two commits were cleanup of my own error | low (self-inflicted, resolved) |
| 2 | genome-anomaly-lab `verified=True` hardcoded, never computed | **material to H4a** |
| 3 | genome-anomaly-lab `bytes_expected`/`bytes_actual` tautological on genome path | **material to H4a** |
| 4 | SDL tests could not be executed here | limits verification |
| 5 | **My provenance validator is syntactic, not evidentiary** | **material — refutes the framing of H4b** |

---

## Discrepancies between my previous claims and current evidence

1. **"provenance architecture explains resistance to fabricated values"** —
   overstated. Phase 3 shows my own validator does not resist fabrication.
2. **"the pattern does not recur"** — overstated. Coverage was 2 of 22,679
   files, and the hand-entered-constant category was not audited.
3. **"enforced per-cell provenance"** — accurate only in the narrow sense that
   no cell lacks a citation. Should not be described as verification.

No scientific result, dataset value, or hypothesis status changed in this audit.
Findings 2–5 affect the *framing* of H4, not any pineal result.
