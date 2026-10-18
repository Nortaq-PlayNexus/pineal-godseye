"""PHASE 3, corrected: is each attack caught by SOURCE CHECKING or by
INTERNAL CONSISTENCY?

The first run of test_provenance_adversarial.py reported attacks 5, 6 and 7 as
CAUGHT. That verdict was right; the accompanying explanation was wrong. They
were caught by `gene_intact_all disagrees` -- an INTERNAL CONSISTENCY check
between a column and its derived column. That is the exact mechanism that let
the original fabrication through: a constraint satisfied by fitting.

This file re-runs each attack and reports WHICH check fired, so the
distinction the brief asks for is explicit:

  SOURCE CHECK      the validator inspected the citation
  INTERNAL CHECK    the validator compared two values against each other
  NOT CAUGHT        nothing fired

The control is attack 5 with the derived column recomputed, so that internal
consistency holds. If the attack is then NOT caught, the honest conclusion is
that no source checking exists at all.
"""
import csv
import io
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
RAW = PROJECT / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]


def classify(rows):
    """Return the set of violation kinds the validator raises."""
    kinds = set()
    for r in rows:
        if not (r.get("source") or "").strip():
            kinds.add("SOURCE_CHECK")
        vals = [(r.get(c) or "") for c in GENES]
        for c, v in zip(GENES, vals):
            if v not in ("", "0", "1"):
                kinds.add("INTERNAL_CHECK")
        expected = ("UNKNOWN" if any(v == "" for v in vals)
                    else "1" if all(v == "1" for v in vals) else "0")
        if (r.get("gene_intact_all") or "") != expected:
            kinds.add("INTERNAL_CHECK")
    return kinds


def load(t):
    return list(csv.DictReader(io.StringIO(t)))


def edit(clean, taxon, col=None, value=None, src=None):
    hdr = clean.splitlines()[0].split(",")
    out = []
    for ln in clean.splitlines():
        if ln.startswith(taxon):
            p = ln.split(",")
            if col:
                p[hdr.index(col)] = value
            if src is not None:
                p[-1] = src
            ln = ",".join(p)
        out.append(ln)
    return "\n".join(out) + "\n"


def redrive(clean, taxon):
    """Recompute gene_intact_all so internal consistency HOLDS."""
    t = edit(clean, taxon, col="AANAT_intact", value="0")
    rows = load(t)
    hdr = t.splitlines()[0].split(",")
    gi = hdr.index("gene_intact_all")
    lines = t.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith(taxon):
            p = ln.split(",")
            vals = [(p[hdr.index(c)] or "") for c in GENES]
            p[gi] = ("UNKNOWN" if any(v == "" for v in vals)
                     else "1" if all(v == "1" for v in vals) else "0")
            lines[i] = ",".join(p)
    return "\n".join(lines) + "\n"


def main():
    clean = RAW.read_text(encoding="utf-8")
    print("PHASE 3 (corrected) -- which check actually fires?\n")
    print(f"[control] shipped CSV -> {classify(load(clean)) or 'clean'}\n")

    cases = [
        ("1. source deleted", edit(clean, "American Alligator (Alligator)", src="")),
        ("2. whitespace-only source",
         edit(clean, "American Alligator (Alligator)", src="   ")),
        ("3. fabricated well-formed citation",
         edit(clean, "American Alligator (Alligator)",
              src="Smith et al. 2019 Nature 574:312-318")),
        ("4. citation to nonexistent file",
         edit(clean, "American Alligator (Alligator)",
              src="see /nonexistent/table_99.csv row 4")),
        ("5. value contradicts a REAL source",
         edit(clean, "American Alligator (Alligator)",
              col="AANAT_intact", value="0")),
        ("5c. same, but derived column recomputed (internal consistency HOLDS)",
         redrive(clean, "American Alligator (Alligator)")),
        ("6. value modified, source unchanged",
         edit(clean, "American Alligator (Alligator)",
              col="MTNR1B_intact", value="0")),
        ("7. value AND source changed consistently",
         edit(edit(clean, "American Alligator (Alligator)",
                   col="AANAT_intact", value="0"),
              "American Alligator (Alligator)",
              src="Emerling 2019 Sci Adv: crocodylians lost AANAT")),
    ]

    by_source = by_internal = missed = 0
    for name, text in cases:
        kinds = classify(load(text))
        if not kinds:
            missed += 1
            print(f"[NOT CAUGHT]  {name}")
            continue
        if "SOURCE_CHECK" in kinds:
            by_source += 1
            print(f"[CAUGHT: source]  {name}")
        else:
            by_internal += 1
            print(f"[CAUGHT: internal only]  {name}")

    print()
    print("=" * 68)
    print(f"caught by SOURCE check       : {by_source}")
    print(f"caught by INTERNAL check only: {by_internal}")
    print(f"not caught at all            : {missed}")
    print()
    print("CONCLUSION")
    print("-" * 68)
    print("The only attack the SOURCE check catches is a MISSING or BLANK")
    print("citation. Every substantive attack -- fabricated citation, dangling")
    print("citation, value contradicting its source -- passes when internal")
    print("consistency is allowed to hold (case 5c).")
    print()
    print("'Per-cell provenance enforcement' is therefore accurate ONLY as:")
    print("    'no shipped cell lacks a source string'")
    print("It is NOT evidence verification. The validator never resolves a")
    print("citation or compares a value to its source.")
    print()
    print("Note case 5: it appeared 'caught' in the first run, but only because")
    print("the derived column had not been updated. That is an internal")
    print("consistency check -- the same class of constraint that the original")
    print("fabricated rows were built to satisfy.")


if __name__ == "__main__":
    main()
