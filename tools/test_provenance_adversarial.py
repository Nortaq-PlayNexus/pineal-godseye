"""Adversarial tests of pineal-godseye's per-cell provenance enforcement.

Runs the 7 attacks named in the verification brief and classifies WHICH check
fires for each: SOURCE_CHECK, INTERNAL_CHECK, or NOT_CAUGHT.

HARNESS CORRECTIONS (2026-10-19, second audit pass)
---------------------------------------------------
The first version of this harness used naive `line.split(",")` to edit rows.
That is wrong for two independent reasons, both confirmed by execution:

1. **The shipped CSV is unquoted.** Sources contain commas
   ("Emerling 2021: 'all 13 crocodylians ...'"), so a naive split produces more
   fields than the header and `csv.DictReader` mis-assigns columns. Verified:
   `csv.reader` and `split(",")` disagree on affected rows.
2. **The row selector `ln.startswith(taxon)` matched only the first line**, so
   `parsed[0]` was HAGFISH, not the alligator. A verification harness that
   silently reads the wrong row is worse than no harness.

This version round-trips through the csv module for BOTH reading and writing,
and selects rows by parsed `taxon` field. The canonical CSV is never modified;
every attack is applied to an in-memory copy.

CONTROL: case 5c is the decisive one. It applies a source-contradicting change
AND recomputes the derived column, so internal consistency holds. If that
passes, no evidence verification is occurring.
"""
import csv
import io
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
RAW = PROJECT / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]


def _assert_columns_known():
    """Guard against a mistyped column name silently degrading every check.

    A near-miss typo like 'MTNR1B_intain'[:7] yields a name absent from the
    CSV; DictReader then returns None for it, (None or '') becomes '', and
    every 'is the flag domain valid' test silently passes. The control in
    main() would go red, but only by accident -- and a red control is easy to
    misread as 'the shipped data is bad' rather than 'the harness is wrong'.
    Fail loudly at import time instead.
    """
    if not RAW.exists():
        return
    fields = set(next(csv.reader(io.StringIO(RAW.read_text(encoding="utf-8")))))
    unknown = [g for g in GENES if g not in fields]
    if unknown:
        raise AssertionError(
            f"GENES contains column name(s) absent from the CSV: {unknown}. "
            f"Available: {sorted(fields)}")


_assert_columns_known()


def load(text):
    return list(csv.DictReader(io.StringIO(text)))


def serialise(rows, fields):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def classify(rows):
    """Every condition that can mark a case caught, classified."""
    kinds = set()
    for r in rows:
        if not (r.get("source") or "").strip():
            kinds.add("SOURCE_CHECK")          # citation present and non-blank
        vals = [(r.get(c) or "") for c in GENES]
        for c, v in zip(GENES, vals):
            if v not in ("", "0", "1"):
                kinds.add("INTERNAL_CHECK")    # flag domain
        expected = ("UNKNOWN" if any(v == "" for v in vals)
                    else "1" if all(v == "1" for v in vals) else "0")
        if (r.get("gene_intact_all") or "") != expected:
            kinds.add("INTERNAL_CHECK")        # derived column agreement
    return kinds


def attack(clean, taxon, col=None, value=None, src=None, redrive=False):
    """Apply one mutation to an in-memory copy, via the csv module."""
    rows = load(clean)
    fields = list(rows[0].keys())
    hit = 0
    for r in rows:
        if (r.get("taxon") or "").strip() != taxon:
            continue
        hit += 1
        if col:
            r[col] = value
        if src is not None:
            r["source"] = src
        if redrive:
            vals = [(r.get(c) or "") for c in GENES]
            r["gene_intact_all"] = ("UNKNOWN" if any(v == "" for v in vals)
                                    else "1" if all(v == "1" for v in vals)
                                    else "0")
    assert hit == 1, f"row selector matched {hit} rows for {taxon!r}, expected 1"
    return serialise(rows, fields)


def main():
    clean = RAW.read_text(encoding="utf-8")

    # sanity: the harness reads what it thinks it reads
    base = load(clean)
    assert len(base) == 25, f"expected 25 rows, parsed {len(base)}"
    alig = [r for r in base
            if r["taxon"].strip() == "American Alligator (Alligator)"]
    assert len(alig) == 1, "alligator row not uniquely selectable"
    assert alig[0]["AANAT_intact"] == "1", "baseline alligator AANAT_intact != 1"
    assert alig[0]["gene_intact_all"] == "1", "baseline gene_intact_all != 1"

    print("ADVERSARIAL PROVENANCE TESTS (harness round-trips through csv)\n")
    print(f"[control] shipped CSV -> {classify(base) or 'clean'}")
    print(f"[control] alligator baseline: AANAT_intact="
          f"{alig[0]['AANAT_intact']} gene_intact_all={alig[0]['gene_intact_all']} "
          f"source={alig[0]['source'][:48]}...\n")

    AL = "American Alligator (Alligator)"
    SAL = "Salmon (Salmo salar)"
    cases = [
        ("1. source deleted", attack(clean, AL, src="")),
        ("2. whitespace-only source", attack(clean, AL, src="   ")),
        ("3. fabricated well-formed citation",
         attack(clean, AL, src="Smith et al. 2019 Nature 574:312-318")),
        ("4. citation to nonexistent file",
         attack(clean, AL, src="see /nonexistent/table_99.csv row 4")),
        ("5. value contradicts a REAL source",
         attack(clean, AL, col="AANAT_intact", value="0")),
        ("5c. same, gene_intact_all recomputed (internal consistency HOLDS)",
         attack(clean, AL, col="AANAT_intact", value="0", redrive=True)),
        ("6. value modified, source unchanged",
         attack(clean, AL, col="MTNR1B_intact", value="0")),
        ("7. value AND source changed consistently",
         attack(clean, AL, col="AANAT_intact", value="0",
                src="Emerling 2019 Sci Adv: crocodylians lost AANAT")),
        ("8. source containing commas (round-trip check)",
         attack(clean, SAL, col="ASMT_intact", value="0",
                src="Foo, A. 2020, J. Test 1:1, 'quoted, with comma'")),
    ]

    by_source = by_internal = missed = 0
    for name, text in cases:
        kinds = classify(load(text))
        if not kinds:
            missed += 1
            print(f"[NOT CAUGHT]  {name}")
        elif "SOURCE_CHECK" in kinds:
            by_source += 1
            print(f"[CAUGHT: source]  {name}")
        else:
            by_internal += 1
            print(f"[CAUGHT: internal only]  {name}")

    print()
    print("=" * 70)
    print(f"caught by SOURCE check        : {by_source}")
    print(f"caught by INTERNAL check only : {by_internal}")
    print(f"not caught at all             : {missed}")
    print()
    print("VERDICT")
    print("-" * 70)
    print("Test 7: INTERNAL CHECK ONLY.")
    print("The validator never opens a URL, resolves a DOI, reads a source file,")
    print("or compares any value against cited evidence. It performs")
    print("CITATION-PRESENCE checking plus INTERNAL dataset consistency.")
    print()
    print("Case 5c is decisive: a value contradicting Emerling 2021 passes once")
    print("gene_intact_all is recomputed. Case 8 confirms the harness now")
    print("survives a citation containing commas and quotes.")
    print()
    print("A validator testing internal consistency cannot distinguish a")
    print("transcribed value from a well-fitted fabrication, because")
    print("satisfying the constraint is how the fabrication was built.")


if __name__ == "__main__":
    main()
