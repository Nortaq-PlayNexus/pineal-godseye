"""Prove the 2026-10-19 provenance edit changed TEXT ONLY.

The audit replaced five source strings with passages confirmed verbatim
against a retrieved source. If any gene-intactness value moved, the edit
would be a scientific change disguised as a bookkeeping one. This asserts it
was not.

Compares the committed CSV (HEAD) against the working tree, cell by cell,
and fails if any column other than `source` differs.
"""
import csv
import io
import subprocess
import sys

NON_PROVENANCE = ["taxon", "clade", "class_group", "AANAT_intact", "ASMT_intact",
                  "MTNR1A_intact", "MTNR1B_intact", "gene_intact_all"]


def parse(text):
    return list(csv.DictReader(io.StringIO(text)))


def main():
    old = subprocess.run(["git", "show", "HEAD:data/raw/vertebrate_pineal_phylogeny.csv"],
                         capture_output=True, text=True, encoding="utf-8",
                         check=True).stdout
    with open("data/raw/vertebrate_pineal_phylogeny.csv", encoding="utf-8") as fh:
        new = fh.read()

    a, b = parse(old), parse(new)
    if len(a) != len(b):
        print(f"FAIL: row count changed {len(a)} -> {len(b)}")
        sys.exit(1)

    by_taxon_old = {r["taxon"].strip(): r for r in a}
    changed_text, changed_value = [], []

    for r in b:
        t = r["taxon"].strip()
        if t not in by_taxon_old:
            print(f"FAIL: new taxon appeared: {t}")
            sys.exit(1)
        o = by_taxon_old[t]
        for col in NON_PROVENANCE:
            if (o.get(col) or "") != (r.get(col) or ""):
                changed_value.append(f"{t}.{col}: {o.get(col)!r} -> {r.get(col)!r}")
        if (o.get("source") or "") != (r.get("source") or ""):
            changed_text.append(t)

    print("=" * 74)
    print("PROVENANCE-EDIT IMPACT CHECK")
    print("=" * 74)
    print(f"rows                         : {len(b)}")
    print(f"rows with source text changed: {len(changed_text)}")
    for t in changed_text:
        print(f"    {t}")
    print(f"VALUE CELLS changed           : {len(changed_value)}")
    for c in changed_value:
        print(f"    {c}")
    print()
    if changed_value:
        print("FAIL -- this was a data edit, not a provenance edit.")
        sys.exit(1)
    if not changed_text:
        print("FAIL -- nothing changed; the edit did not apply.")
        sys.exit(1)
    print("PASS -- provenance text only. No gene-intactness value moved.")
    print("       No hypothesis status changes on this basis.")


if __name__ == "__main__":
    main()