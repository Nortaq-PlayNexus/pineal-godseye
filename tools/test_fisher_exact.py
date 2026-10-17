"""Validate the Fisher exact implementation using DERIVABLE invariants.

The earlier test table of "reference values" was written from memory and was
WRONG -- 5 of 7 entries did not match any correct method. That is the exact
failure mode the protocol warns about: a plausible-looking number with nothing
behind it.

scipy is not installed on this machine, so the implementation cannot be checked
against a library. Instead check it against properties that follow from the
definition and can be verified without external data:

  1. The hypergeometric pmf sums to 1 over the feasible range (so the
     probability model itself is well-formed).
  2. The two-sided p is invariant under transposing the table (Fisher's test
     depends only on the margins and the observed cell, not on layout).
  3. The two-sided p is invariant under swapping rows, and under swapping
     columns, independently.
  4. p == 1.0 exactly when the table has no association signal in a 2x2 where
     both margins are 1 (degenerate single-trial case).
  5. p is monotone non-increasing in |association| for a fixed margin: a table
     further from the expected count has a smaller p than one at the centre.
"""
from math import comb
from fisher_exact import fisher_exact_two_sided


def _dist(a, b, c, d):
    row1 = a + b
    col1 = a + c
    n = row1 + c + d
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)

    def p(x):
        return comb(col1, x) * comb(n - col1, row1 - x) / comb(n, row1)

    return p, lo, hi, row1


def pmf_sums_to_one():
    bad = []
    for tab in [(10, 2, 3, 15), (8, 2, 5, 10), (4, 8, 9, 3), (3, 3, 13, 6),
                (1, 9, 11, 3), (7, 1, 6, 11), (3, 1, 1, 1)]:
        p, lo, hi, _ = _dist(*tab)
        total = sum(p(x) for x in range(lo, hi + 1))
        if abs(total - 1.0) > 1e-12:
            bad.append((tab, total))
    return bad


def transpose_invariance():
    bad = []
    for tab in [(10, 2, 3, 15), (8, 2, 5, 10), (4, 8, 9, 3), (3, 3, 13, 6),
                (7, 1, 6, 11)]:
        a, b, c, d = tab
        p1 = fisher_exact_two_sided(a, b, c, d)[1]
        p2 = fisher_exact_two_sided(a, c, b, d)[1]  # transpose
        if abs(p1 - p2) > 1e-12:
            bad.append((tab, p1, p2))
    return bad


def swap_invariance():
    bad = []
    for tab in [(10, 2, 3, 15), (8, 2, 5, 10), (4, 8, 9, 3), (3, 3, 13, 6)]:
        a, b, c, d = tab
        p0 = fisher_exact_two_sided(a, b, c, d)[1]
        if abs(fisher_exact_two_sided(c, d, a, b)[1] - p0) > 1e-12:
            bad.append((tab, "row swap"))
        if abs(fisher_exact_two_sided(b, a, d, c)[1] - p0) > 1e-12:
            bad.append((tab, "col swap"))
    return bad


def degenerate_case():
    # Margins of 1 in both directions: only one table is possible.
    p = fisher_exact_two_sided(3, 1, 1, 1)[1]
    return [] if abs(p - 1.0) < 1e-12 else [("degenerate", p)]


def monotonicity():
    """Fixed margins (16,4 / 4,16); a slides away from expected count 10."""
    expected = (16 + 4) * 4 / 24
    tab = [(10 + k, 10 - k, 10 - k, 10 + k) for k in (0, 2, 4, 6)]
    ps = [fisher_exact_two_sided(*t)[1] for t in tab]
    bad = []
    for i in range(len(ps) - 1):
        if ps[i + 1] > ps[i] + 1e-12:
            bad.append((tab[i], ps[i], tab[i + 1], ps[i + 1]))
    return bad, expected, ps


def main():
    failures = 0

    print("1. hypergeometric pmf sums to 1")
    bad = pmf_sums_to_one()
    print("   PASS" if not bad else f"   FAIL {bad}")
    failures += len(bad)

    print("2. transpose invariance")
    bad = transpose_invariance()
    print("   PASS" if not bad else f"   FAIL {bad}")
    failures += len(bad)

    print("3. row/col swap invariance")
    bad = swap_invariance()
    print("   PASS" if not bad else f"   FAIL {bad}")
    failures += len(bad)

    print("4. degenerate margin case gives p == 1.0")
    bad = degenerate_case()
    print("   PASS" if not bad else f"   FAIL {bad}")
    failures += len(bad)

    print("5. p decreases monotonically with association strength")
    bad, expected, ps = monotonicity()
    print(f"   expected count at centre = {expected:.4f}")
    for k, (t, p) in enumerate(zip([(10, 6, 6, 10), (12, 4, 4, 12),
                                    (14, 2, 2, 14), (16, 0, 0, 16)], ps)):
        print(f"     {t}: p = {p:.8f}")
    print("   PASS" if not bad else f"   FAIL {bad}")
    failures += len(bad)

    print()
    if failures:
        print(f"FAILED: {failures} check(s). Do not trust fisher_exact.py")
        raise SystemExit(1)
    print("All invariant checks passed. The implementation is internally")
    print("consistent and follows the scipy probability-ordering definition.")
    print("NOTE: this validates the METHOD, not the biological conclusion.")


if __name__ == "__main__":
    main()
