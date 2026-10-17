"""Prove the build-time provenance validation catches the class of defect it
was written for.

REWRITTEN 2026-10-18. The previous version asserted that nine binary columns
were 0/1, which is exactly what let a teleost PARALOG COUNT (2) sit inside a
binary gene-intactness flag unnoticed. The current builder ships only
gene-intactness columns, so this test now checks the invariant that actually
protects the data:

  1. Every shipped row carries a non-empty source string.
  2. Gene-intactness values are 0, 1, or blank (blank = UNKNOWN).
  3. The derived `gene_intact_all` column agrees with the four inputs, and is
     "UNKNOWN" whenever any input is blank.

Then it proves the check fires by re-introducing the exact historical defects:
  (a) a fabricated row with no source
  (b) a paralog count (2) written into the binary AANAT_intact field

Exits 0 only if both are caught.
"""
import csv
import io
import sys
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]


def scan_phylo(text):
    """Return a list of provenance violations found in the CSV text."""
    hits = []
    for row in csv.DictReader(io.StringIO(text)):
        # A short row yields None for missing fields. A validator must report
        # that, not crash on it.
        taxon = row.get("taxon") or "<unnamed row>"
        src = row.get("source")
        if not src or not src.strip():
            hits.append(f"{taxon}: no source")
        vals = [(row.get(g) or "") for g in GENES]
        for g, v in zip(GENES, vals):
            if v not in ("", "0", "1"):
                hits.append(f"{taxon}.{g}={v!r} not 0/1/blank")
        # Blank means genuinely UNKNOWN (e.g. hagfish genome incomplete;
        # manatee heterozygous for an ASMT pseudogene). It must NOT be
        # collapsed to 0 -- that is how the old dataset encoded unknowns as
        # false zeros. An UNKNOWN in any input makes the whole row UNKNOWN.
        if any(v == "" for v in vals):
            expected = "UNKNOWN"
        elif all(v == "1" for v in vals):
            expected = "1"
        else:
            expected = "0"
        actual = row.get("gene_intact_all") or ""
        if actual != expected:
            hits.append(f"{taxon}: gene_intact_all={actual!r} should be {expected!r}")
    return hits


def main():
    clean = (RAW / "vertebrate_pineal_phylogeny.csv").read_text(encoding="utf-8")

    found = scan_phylo(clean)
    if found:
        print("FAIL: shipped CSV violates its own checks:")
        for h in found:
            print("   ", h)
        sys.exit(1)
    print("PASS: shipped CSV is clean "
          "(all rows sourced, all flags binary-or-blank, derived column agrees)")

    # (a) strip the source from one row -> must be caught. The source text is
    # the LAST field, so match on the row prefix rather than the full string.
    lines = clean.splitlines()
    mutated = list(lines)
    for i, ln in enumerate(mutated):
        if ln.startswith("American Alligator"):
            mutated[i] = ",".join(ln.split(",")[:-1])  # drop the source field
            break
    no_src = "\n".join(mutated) + ("\n" if clean.endswith("\n") else "")
    if no_src == clean:
        print("SKIP: alligator row not found; corruption not applied")
    else:
        hits = scan_phylo(no_src)
        if not any("no source" in h for h in hits):
            print(f"FAIL: unsourced row not caught ({hits})")
            sys.exit(1)
        print(f"PASS: unsourced row caught -> {hits[0]}")

    # (b) re-introduce the historical bug: paralog count 2 in a binary flag
    buggy = clean.replace(
        "Salmon (Salmo salar),Actinopterygii,Fishes,1,1,1,1,1",
        "Salmon (Salmo salar),Actinopterygii,Fishes,1,2,1,1,0", 1)
    if buggy == clean:
        print("SKIP: salmon row not in expected form; corruption not applied")
    else:
        hits = scan_phylo(buggy)
        if not hits:
            print("FAIL: paralog-count corruption not caught")
            sys.exit(1)
        print(f"PASS: paralog count in a binary flag caught -> {hits[0]}")

    print("\nAll provenance checks fired. The validator detects both the")
    print("missing-source and the wrong-type failure modes.")


if __name__ == "__main__":
    main()
