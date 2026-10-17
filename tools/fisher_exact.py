"""Fisher's exact test for the 2x2 table, implemented in stdlib only.

scipy is not installed on this machine (verified 2026-10-18:
ModuleNotFoundError: No module named 'scipy'). Rather than skip the statistic
or report a number from memory, this implements the exact hypergeometric
two-sided p-value.

Method: given a 2x2 table [[a,b],[c,d]] with fixed margins, the probability of
observing a table with the same margins under the null (independence) is the
hypergeometric pmf. The two-sided p-value sums the pmf over every table whose
probability is <= the observed table's probability (the standard "probability
ordering" definition, which is what scipy.stats.fisher_exact implements with
its default alternative='two-sided').

Row 1 = unihemispheric sleep present; columns = [melatonin pathway intact,
melatonin pathway lost].
"""
from math import comb


def fisher_exact_two_sided(a, b, c, d):
    """Return (odds_ratio, p_value) for table [[a,b],[c,d]].

    a = sleep=1 & intact, b = sleep=1 & lost
    c = sleep=0 & intact, d = sleep=0 & lost
    """
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    n = row1 + row2

    def prob(x):
        """P(table with row1=x) under hypergeometric(col1, row1, n)."""
        return comb(col1, x) * comb(n - col1, row1 - x) / comb(n, row1)

    p_obs = prob(a)
    # Sum pmf over all feasible x with prob(x) <= p_obs (with small tolerance
    # for float accumulation error).
    p_value = 0.0
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)
    for x in range(lo, hi + 1):
        px = prob(x)
        if px <= p_obs + 1e-12:
            p_value += px
    p_value = min(p_value, 1.0)

    # Odds ratio with Haldane-Anscombe correction if any cell is zero.
    if 0 in (a, b, c, d):
        odds = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
    else:
        odds = (a * d) / (b * c)

    return odds, p_value


def main():
    # From the shipped dataset (verified 2026-10-18).
    # USWS present: intact=3, lost=3. USWS absent: intact=13, lost=6.
    a, b, c, d = 3, 3, 13, 6
    odds, p = fisher_exact_two_sided(a, b, c, d)
    print("2x2 table (rows: USWS present / absent; cols: intact / lost)")
    print(f"  present: intact={a}  lost={b}")
    print(f"  absent:  intact={c}  lost={d}")
    print(f"n = {a + b + c + d}")
    print(f"odds ratio = {odds:.4f}")
    print(f"Fisher exact two-sided p = {p:.4f}")
    print()
    if p > 0.05:
        print("p > 0.05: NO significant association between melatonin-pathway")
        print("state and unihemispheric sleep in this 25-taxon sample.")
        print("This is a NULL result and is reported as such.")
    else:
        print("p <= 0.05: association detected (would warrant follow-up).")
    print()
    print("Caveats: n=25, single 2x2, taxa not independent (phylogenetic")
    print("structure ignored), and the unihemispheric_sleep flag is a")
    print("hand-assigned literature annotation with two known weak rows.")


if __name__ == "__main__":
    main()
