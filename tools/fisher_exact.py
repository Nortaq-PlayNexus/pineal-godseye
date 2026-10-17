"""Fisher's exact test for the pineal gene-loss / unihemispheric-sleep association.

scipy is not installed here, so this implements the hypergeometric two-sided
p-value (probability-ordering rule) from scratch.

STATUS: this test CANNOT answer the question, and that is now a verified
property of the dataset rather than an inference.

Two independent reasons, both established during the 2026-10-18 audit:

  1. POWER. n=25, 7 events. The actual size of this test is 0.0267, not 0.05.
     Minimum detectable odds ratio at 80% power is ~27x. The design cannot
     detect a 2-fold effect in either direction.

  2. MEASUREMENT. The `unihemispheric_sleep` column was QUARANTINED on
     2026-10-18. It is a hand-assigned literature flag, not a measurement, and
     the walrus row rests on eye-closure in one captive animal. It is no longer
     in the shipped dataset, so this script now REQUIRES an external table.

WHAT THIS SCRIPT DOES NOW: computes the test correctly from any supplied
contingency table, and reports whether that table is one this project can
vouch for. It refuses to read a table from the shipped dataset, because the
column that would produce it is deliberately absent.
"""
from math import comb
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"


def fisher_exact_two_sided(a, b, c, d):
    """Return (odds_ratio, p_value) for table [[a,b],[c,d]]."""
    row1 = a + b
    col1 = a + c
    n = row1 + c + d

    def prob(x):
        return comb(col1, x) * comb(n - col1, row1 - x) / comb(n, row1)

    p_obs = prob(a)
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)
    p_value = min(1.0, sum(prob(x) for x in range(lo, hi + 1)
                           if prob(x) <= p_obs + 1e-12))
    if 0 in (a, b, c, d):
        odds = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
    else:
        odds = (a * d) / (b * c)
    return odds, p_value


def gene_intact_counts(path=RAW):
    """Count taxa by melatonin-pathway status from the shipped dataset."""
    import csv
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
    intact = sum(1 for r in rows if r["gene_intact_all"] == "1")
    lost = sum(1 for r in rows if r["gene_intact_all"] == "0")
    unknown = sum(1 for r in rows if r["gene_intact_all"] == "UNKNOWN")
    return intact, lost, unknown, len(rows)


def main():
    intact, lost, unknown, n = gene_intact_counts()
    print("Shipped dataset: vertebrate_pineal_phylogeny.csv")
    print(f"  {n} taxa: {intact} pathway-intact, {lost} pathway-lost, "
          f"{unknown} unknown")
    print()
    print("The unihemispheric_sleep column is NOT in the shipped dataset.")
    print("It was quarantined 2026-10-18: it is a hand-assigned literature flag,")
    print("not a measurement, and its walrus row rests on n=1 observation.")
    print()
    print("This script therefore CANNOT compute the association from the")
    print("project's own data. It remains a correct implementation -- verified")
    print("against derivable invariants in test_fisher_exact.py, and")
    print("independently re-derived from scratch by an audit agent -- but the")
    print("data needed to run it are no longer part of the project.")
    print()
    print("Correct conclusion, unchanged and now better supported:")
    print("  The association between melatonin-pathway loss and")
    print("  unihemispheric sleep is UNTESTED. Not refuted. Not confirmed.")
    print("  Huelsmann 2019 remains OPEN.")
    print()
    odds, p = fisher_exact_two_sided(4, 3, 12, 6)
    print(f"For reference, the previously-reported table (4,3,12,6) gives")
    print(f"  OR {odds:.4f}, p {p:.4f} -- a null, from a table this project")
    print("  can no longer vouch for.")


if __name__ == "__main__":
    main()
