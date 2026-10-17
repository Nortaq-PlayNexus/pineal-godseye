"""Fisher's exact test for the pineal gene-loss / unihemispheric-sleep table.

scipy is not installed on this machine (verified 2026-10-18), so this implements
the exact hypergeometric two-sided p-value from scratch.

IMPLEMENTATION
--------------
Given a 2x2 table [[a,b],[c,d]], the probability of any table with the same
margins under independence is the hypergeometric pmf. The two-sided p-value sums
the pmf over every feasible table whose probability is <= the observed table's
(the probability-ordering rule scipy.stats.fisher_exact uses by default).

VERIFIED: matches exact-Fraction reimplementation to 6 dp on 8 tables, and
passes the derivable invariants in test_fisher_exact.py.

WHY THE TABLE IS COMPUTED, NOT HARDCODED
-----------------------------------------
An earlier version hardcoded `a, b, c, d = 3, 3, 13, 6` and labelled it "from
the shipped dataset (verified 2026-10-18)". It was not: that is the
walrus-EXCLUDED table from validate_usws_flag.py. The shipped CSV yields
4, 3, 12, 6. Both the odds ratio and the p-value were wrong as reported
(0.4615/0.6300 vs the correct 0.6667/0.6729), and the power analysis was
computed on the wrong event count (E=6 vs E=7).

Hardcoding a contingency table next to a claim it is supposed to support is how
a number stops describing the artefact it is attributed to. So the table is
derived from the CSV every run.
"""
from math import comb
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]


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
    p_value = 0.0
    for x in range(lo, hi + 1):
        if prob(x) <= p_obs + 1e-12:
            p_value += prob(x)
    p_value = min(p_value, 1.0)

    if 0 in (a, b, c, d):
        odds = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
    else:
        odds = (a * d) / (b * c)
    return odds, p_value


def table_from_dataset(path=RAW):
    """Compute the 2x2 table directly from the shipped phylogeny CSV.

    Rows    = unihemispheric sleep present / absent
    Columns = melatonin pathway intact / lost
    """
    import csv

    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
    a = b = c = d = 0
    for r in rows:
        intact = all(r[g] == "1" for g in GENES)
        usws = r["unihemispheric_sleep"] == "1"
        if usws and intact:
            a += 1
        elif usws:
            b += 1
        elif intact:
            c += 1
        else:
            d += 1
    return a, b, c, d, len(rows)


def main():
    a, b, c, d, n = table_from_dataset()
    odds, p = fisher_exact_two_sided(a, b, c, d)

    print(f"Table COMPUTED FROM {RAW.name} ({n} taxa)")
    print("=" * 60)
    print("                          intact   lost")
    print(f"  unihemispheric sleep       {a:4d}   {b:4d}")
    print(f"  bihemispheric only         {c:4d}   {d:4d}")
    print()
    print(f"  n = {n}   events = {a + b}")
    print(f"  odds ratio = {odds:.4f}")
    print(f"  Fisher exact two-sided p = {p:.4f}")
    print()
    if p > 0.05:
        print("  NULL: no significant association in this sample.")
    else:
        print("  Association detected -- follow up required.")
    print()
    print("  CAVEATS: n is small, phylogenetic independence is ignored, and")
    print("  the unihemispheric_sleep column is a hand-assigned literature")
    print("  flag with at least one row (walrus) resting on n=1 observation.")
    print("  See tools/power_analysis.py -- this design cannot test the")
    print("  hypothesis in the direction Huelsmann 2019 proposes.")


if __name__ == "__main__":
    main()
