# Audit: citation verification â€” is the provenance limit fixable? (2026-10-19)

Third pass on the question the first two audits raised. Both established that
the shipped validator is **syntactic, not evidentiary**: it checks a citation
string is non-blank and that derived columns agree with their inputs. It never
resolves a citation.

This pass asked whether that limit is **inherent or fixable**, tested it, and
found that the answer changed three times while I was measuring it.

---

## Finding 1 â€” my own probe's first result (8/8) was an artifact

`tools/probe_is_limit_fixable.py` v1 modelled an "upgraded validator" with an
evidence table mapping each citation to a verbatim passage. It scored **8/8**
attacks caught.

That number was worthless. The evidence entries had a `quote` field, and
`classify_upgraded()` **never read it**. The upgrade reduced to:

```python
key = next((k for k, v in evidence.items()
            if v["citation"].split()[0] in src), None)
```

â€” an allowlist of author first-words (`"Emerling"` appearing anywhere in the
source string), plus a hand-written expected value. No sentence from any paper
was consulted. That is the same syntactic check the two prior audits already
refuted, renamed.

Second defect, same file: the table covered **2 of 25 rows**, and every attack
targeted one of those two rows. Detection on the covered subset said nothing
about the other 23.

---

## Finding 2 â€” I fabricated a quotation, reproducing the original sin

While building the evidence table I wrote, from memory:

> "all 13 crocodylians investigated retain intact melatonin genes"

and attributed it to Emerling et al. 2021.

I retrieved the actual abstract (Europe PMC, DOI `10.12688/openreseurope.13795.2`).
The paper says:

> "We report that **crocodylians retain intact melatonin genes** and express
> AANAT and ASMT in their eyes"

The string I wrote **does not occur**. The number "13" does not appear anywhere
in the retrieved abstract. The paper examined **123** vertebrate species.

**The value (all four genes intact) is correct.** The *quotation* was not. And
the fabricated string was already in the shipped dataset â€” I reproduced it
independently, from memory, without checking.

That is the finding. The defect is not one stale row; it is in the process
that produces citation strings. Nothing in the build ever asked whether a
quoted passage exists.

---

## Finding 3 â€” 0 of 8 shipped quotations were verifiable

`tools/check_quotes_against_source.py` now **extracts** every single-quoted
fragment from the source column automatically (v1 hand-listed them, which is
precisely what would miss a row) and tests each against retrieved text.

Before the fix:

```
quotation claims found : 8
  VERBATIM in source   : 0
  NOT FOUND            : 7
  NO SOURCE (no DOI)   : 1
rows asserting a value with NO source at all : 11/25
```

**Caveat, stated up front:** Europe PMC serves only the abstract for this
record (no PMCID, `isOpenAccess=N`). A genuine sentence from the paper body
would also read NOT FOUND. So NOT FOUND means *not quotable from retrievable
text* â€” not *proven false*.

What holds regardless of that caveat: the dataset places these strings inside
quotation marks, asserting they are verbatim. Nothing verified that assertion,
and the audit could not confirm it.

---

## Finding 4 â€” what I changed, and what I deliberately did not

Five rows had a passage in the retrieved abstract that directly supports the
recorded value. Those now quote it, with the DOI:

| Taxon | Now quotes (VERBATIM) |
|---|---|
| American Alligator | `crocodylians retain intact melatonin genes` |
| Chinese Pangolin | `all four genes have been repeatedly inactivated in the pineal-less xenarthrans, pangolins, sirenians, and whales` |
| Two-Toed Sloth | *(same sentence)* |
| Nine-Banded Armadillo | *(same sentence)* |
| Sunda Colugo | `colugos have lost these genes` |

**Three left unfixed on purpose** â€” I could not verify them, and rewriting an
unverifiable quotation as a paraphrase would substitute one unverifiable claim
for another:

- Human â€” `'five primates'` not found
- Naked Mole-Rat â€” `'both melatonin receptor genes appear to be inactivated in â€¦'` not found
- Bottlenose Dolphin â€” Huelsmann 2019, no DOI recorded, **NO SOURCE**

After the fix: **5 VERBATIM / 2 NOT FOUND / 1 NO SOURCE.**

`tools/verify_provenance_edit_is_text_only.py` diffs the change cell by cell
against `HEAD` and asserts no column but `source` moved. Result: **0 value cells
changed.** This was a provenance edit, not a data edit. **No hypothesis status
changes.**

Originals preserved in `quarantine/UNVERIFIED_QUOTATIONS_2026-10-19.csv`; the
pre-change CSV hashed `5DBE4BBFâ€¦17A25`.
---

## Finding 5 â€” a latent harness defect, exposed by the fix

Correctly quoting a source string that contains a comma **broke**
`tools/regression_test.py`:

```
FAIL: unsourced row not caught ([])
```

The cause was in the test, not the validator. `regression_test.py:78` did:

```python
mutated[i] = ",".join(ln.split(",")[:-1])   # drop the source field
```

That assumes the source field contains no comma. Once a source string was
properly CSV-quoted *and* contained a comma, the naive split dropped only the
final comma-segment â€” leaving a row whose source was still **non-empty**. The
validator correctly reported nothing, and the test failed.

This is the **same defect class** already found and fixed in
`tools/test_provenance_adversarial.py` two turns earlier. It was never fixed
here. It stayed latent only because no shipped source string contained a comma
until this edit created one.

Fixed by round-tripping through the `csv` module. `regression_test.py` exits 0.

> A harness defect that cannot fire is indistinguishable from a harness that
> works. This one survived three commits and only surfaced when the data
> changed in a way it had never seen.

---

## Finding 6 â€” the upgrade is real, and strictly better

Probe v3 scores each case on the **target row alone** (v2 aggregated over all 25
rows, so 23 rows' coverage verdict masked whether the target was caught), and
**adds** to the current checks rather than replacing them.

| # | Case | current | v3 | mode | verdict |
|---|---|---|---|---|---|
| 1 | source deleted | CAUGHT | CAUGHT | n/a | SOURCE_CHECK |
| 2 | fabricated citation swapped in | MISSED | **CAUGHT** | n/a | UNRESOLVED_CITATION |
| 3 | value contradicts, derived recomputed | MISSED | **CAUGHT** | name | EVIDENCE_CONFLICT |
| 4 | value flipped + prose flipped, DOI kept | MISSED | **CAUGHT** | doi | EVIDENCE_CONFLICT |
| 5 | value contradicts, derived not recomputed | CAUGHT | CAUGHT | n/a | INTERNAL_CHECK |
| 6 | correct DOI, correct value, prose asserts opposite | MISSED | **MISSED** | doi | â€” |
| 7 | value contradicts (uncovered row) | MISSED | NOCOV | n/a | UNCOVERED |
| 8 | fabricated citation (uncovered row) | MISSED | NOCOV | n/a | UNCOVERED |

**Current 2/8 â†’ upgraded 5/8.**

Two corrections to my own earlier readings, both recorded rather than quietly
fixed:

- **v2 scored 4/8 because the upgrade had *replaced* the internal-consistency
  check instead of adding to it** â€” losing case 5, which the shipped validator
  caught. An upgrade that is not a strict superset is a regression. After the
  fix: 5/8, and it catches everything the current validator catches.
- **Case 4 was not the attack it claimed to be.** It flipped the value without
  recomputing the derived column, making it a duplicate of case 5. The real
  "value AND citation flipped consistently" case needs `redrive=True`; with it,
  the upgrade **does** catch it. That was a mis-specified test, not a hole.

The **MODE** column matters: the shipped dataset records **no DOIs**, so every
row resolves by author-year name only. That mapping is an assumption in the
probe file, not evidence in the data.

---

## Finding 7 â€” the one hole that remains, and it is the one that bit us

Case 6: keep the correct DOI, keep the correct value, and write prose asserting
the opposite. The citation is used as a **key**; its wording is never compared
to the passage.

That is a **misquotation** defect, not a wrong-value defect â€” and it is exactly
the defect this audit found in the shipped CSV. No string-matching validator
closes it; it would require parsing the citation's assertion and comparing it to
the passage's meaning.

**Coverage remains the binding constraint:** 11 of 25 rows cite `NOT VERIFIED`,
and 3 of 8 quotation claims are still unverifiable.

---

## Verdict

**The limit is not inherent. It is a coverage problem plus a checking problem,
and it is not solved.**

- Evidence-pinning catches the attacks the shipped validator misses â€” but only
  on rows that have a verified passage. **2 of 25.**
- The residual, unfixable-by-matching failure is misquotation.
- Three claims remain unverifiable with the access available. They are reported
  as failures, not reworded into soft positives.

**Superseding the previous audit's framing:** the earlier conclusion was that
the validator *cannot* do evidence verification. That is too strong. It cannot
do it *today*, with 2 covered rows and no DOIs. The mechanism is implementable
and is implemented in the probe. What is missing is coverage.

**Nothing here changes any scientific claim.** Zero value cells changed. H1/H2
remain retired, H3 untested, H4 partially supported, H5 supported.

---

## Concurrent-session collision â€” unresolved

Another session is writing to this repository at the same time. During this
audit it:

- **overwrote `AUDIT_TEST7_2026-10-19.md`** (my version is intact at `46e0bbb`)
- created `TESTWRITE.md` (6 bytes, contents `# Test` â€” a write-size probe)
- created `tools/repro_test7_standalone.py`

I did **not** stage those three. My commit is path-scoped. They remain
untracked/modified in the working tree. **Someone should decide what to do with
them** â€” I did not delete another session's work.

## Test suite

```
build_datasets                       EXIT 0
verify_fractions                     EXIT 0
regression_test                      EXIT 0
test_fisher_exact                    EXIT 0
test_provenance_adversarial          EXIT 0
verify_provenance_edit_is_text_only  EXIT 0
probe_is_limit_fixable               EXIT 0
check_quotes_against_source          EXIT 1   <- correct: 3 claims unverified
```

The last exit 1 is a **real failure**, not a flake. It stays red until the three
remaining claims are verified against retrievable full text.
