"""Standalone reproduction of adversarial test 7.

Isolated from the harness. Applies the test-7 mutation to an in-memory copy
of the shipped CSV and runs the SAME classify() logic the harness uses, so
the classification is reproducible without depending on the harness file.

Two runs:
  RUN A — mutate AANAT_intact + source, leave gene_intact_all unchanged.
  RUN B — same mutations PLUS recompute gene_intact_all from the four
          component columns (the 5c-style control).

The canonical CSV is never modified.
"""
import csv
import io
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]
TAXON = "American Alligator (Alligator)"


def classify(rows):
    """Same three conditions as test_provenance_adversarial.classify."""
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


def run(label, redrive):
    clean = RAW.read_text(encoding="utf-8")
    rows = list(csv.DictReader(io.StringIO(clean)))
    fields = list(rows[0].keys())
    n = 0
    target = None
    for r in rows:
        if (r.get("taxon") or "").strip() != TAXON:
            continue
        n += 1
        target = r
        # Test 7 mutations
        r["AANAT_intact"] = "0"
        r["source"] = "Emerling 2019 Sci Adv: crocodylians lost AANAT"
        if redrive:
            vals = [(r.get(c) or "") for c in GENES]
            r["gene_intact_all"] = ("UNKNOWN" if any(v == "" for v in vals)
                                    else "1" if all(v == "1" for v in vals)
                                    else "0")
    assert n == 1, f"selector matched {n} rows"

    kinds = classify(rows)
    print(f"--- {label} ---")
    print(f"  taxon              = {target['taxon']}")
    print(f"  AANAT_intact       = {target['AANAT_intact']}")
    print(f"  ASMT_intact        = {target['ASMT_intact']}")
    print(f"  MTNR1A_intact      = {target['MTNR1A_intact']}")
    print(f"  MTNR1B_intact      = {target['MTNR1B_intact']}")
    print(f"  gene_intact_all    = {target['gene_intact_all']}")
    print(f"  source             = {target['source']}")
    print(f"  checks fired       = {sorted(kinds) or ['NONE (NOT CAUGHT)']}")
    print()


def main():
    clean = RAW.read_text(encoding="utf-8")
    base = list(csv.DictReader(io.StringIO(clean)))
    alig = [r for r in base if r["taxon"].strip() == TAXON][0]

    print("STANDALONE TEST 7 REPRODUCTION\n")
    print(f"Baseline alligator row:")
    print(f"  AANAT_intact    = {alig['AANAT_intact']}")
    print(f"  gene_intact_all = {alig['gene_intact_all']}")
    print(f"  source          = {alig['source'][:60]}...")
    print()

    print("RUN A: mutate AANAT_intact=0 + source, gene_intact_all UNCHANGED")
    run("RUN A (gene_intact_all still = 1, inconsistent)", redrive=False)

    print("RUN B: same mutations + recompute gene_intact_all from components")
    run("RUN B (gene_intact_all recomputed = 0, internally consistent)",
        redrive=True)

    print("CONCLUSION")
    print("-" * 60)
    print("RUN A fires INTERNAL_CHECK because gene_intact_all (1) disagrees")
    print("with the component columns (one is now 0).")
    print()
    print("RUN B fires NOTHING. The fabricated value AANAT_intact=0 — which")
    print("contradicts Emerling 2021's statement that all 13 crocodylians")
    print("retain intact melatonin genes — passes every implemented check.")
    print()
    print("The validator checks citation PRESENCE and dataset INTERNAL")
    print("consistency. It does not resolve the citation, retrieve the")
    print("cited source, or compare the value against source evidence.")


if __name__ == "__main__":
    main()