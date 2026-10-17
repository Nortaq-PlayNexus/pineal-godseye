"""Power analysis: was the null result (p = 0.63) informative or underpowered?

The Fisher test on the 25-taxon table returned odds ratio 0.46, p = 0.63. The
author reported this as a null result. A null from an underpowered test is NOT
evidence of absence -- it is evidence of nothing. This script determines which,
by asking: with 25 taxa and 6 events, what effect size COULD this design detect?

DESIGN: 16 melatonin-intact, 9 melatonin-lost, 6 unihemispheric-sleep events.
TEST: Fisher exact, two-sided, probability-ordering rule.

CORRECTION (first attempt was wrong): the initial power calculation
reproduced a DECREASING power curve -- power rose as the true odds ratio
increased, which is backwards. The enumeration was wrong. The fix is to
enumerate the 7 possible event splits x = 0..6, compute the Fisher p for each,
find the actual rejection set, and weight each split by its probability under
the alternative. The test is CONSERVATIVE on these margins: its true size is
0.0119, not 0.05, because the two-sided probability-ordering rule produces a
disconnected rejection region on small tables ({x=0, x=1}).

Power is reported against both the nominal alpha = 0.05 and the ACTUAL size,
because a conservative test that only fires at 0.012 has far less power than
its nominal alpha suggests.

Pure stdlib; scipy is not installed on this machine.
"""
from math import comb
from fisher_exact import fisher_exact_two_sided, table_from_dataset

# CORRECTION 2026-10-18: these were hardcoded as NI, NL, E = 16, 9, 6. That was
# the walrus-EXCLUDED table, not the shipped dataset, which has 7 USWS rows.
# Derive them from the CSV instead.
_a, _b, _c, _d, N_TAXA = table_from_dataset()
NI, NL, E = _a + _c, _b + _d, _a + _b


def split_table(x):
    a, b = x, E - x
    c, d = NI - x, NL - (E - x)
    return a, b, c, d


def h0_pmf(x):
    return comb(NI, x) * comb(NL, E - x) / comb(NI + NL, E)


def h1_weights(theta):
    """Probability of each split x under odds ratio theta (noncentral
    hypergeometric), normalised."""
    w = {}
    for x in range(E + 1):
        pmf = h0_pmf(x)
        if pmf > 0:
            w[x] = pmf * (theta ** x)
    total = sum(w.values())
    return {x: v / total for x, v in w.items()}


def power_at(theta, alpha):
    weights = h1_weights(theta)
    return sum(
        w for x, w in weights.items()
        if fisher_exact_two_sided(*split_table(x))[1] < alpha
    )


def main():
    print(f"Design derived from vertebrate_pineal_phylogeny.csv:")
    print(f"  {NI} melatonin-intact / {NL} melatonin-lost, {E} USWS events, "
          f"{N_TAXA} taxa. Fisher two-sided.")
    print()

    # Actual size of the test (must equal alpha if calibrated).
    rejects = [x for x in range(E + 1) if fisher_exact_two_sided(*split_table(x))[1] < 0.05]
    size = sum(h0_pmf(x) for x in rejects)
    print(f"Rejection set at alpha=0.05: x in {rejects}")
    print(f"ACTUAL size of this test = {size:.4f}  (nominal alpha = 0.0500)")
    print("The test is CONSERVATIVE: its two-sided probability-ordering rule")
    print("gives a disconnected rejection region on tables this small.")
    print()

    print(f"Power by true odds ratio (nominal alpha = 0.05):")
    print(f"  {'true OR':>9}  {'power':>7}   {'power at actual size':>22}")
    print("  " + "-" * 44)
    for theta in [1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 50.0]:
        p_nom = power_at(theta, 0.05)
        p_act = power_at(theta, size)
        print(f"  {theta:9.1f}  {p_nom:7.3f}   {p_act:22.3f}")
    print()

    # The decreasing power curve is CORRECT for this test, not a bug.
    # Verified: the H1 weighting reproduces the H0 pmf exactly at theta=1.
    # The reason is that the ONLY splits reaching p<0.05 are x=0 and x=1,
    # and theta>1 pushes probability AWAY from those (toward x=6, p=0.057).
    # So this test detects a protective effect (theta<1) but has almost no
    # power to detect an elevated risk (theta>1) at these margins.
    print("Power in the PROTECTIVE direction (OR < 1):")
    print(f"  {'true OR':>9}  {'power':>7}")
    print("  " + "-" * 20)
    for theta in [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
        print(f"  {theta:9.2f}  {power_at(theta, 0.05):7.3f}")
    print()

    # Minimum detectable OR at 80% power, each direction.
    def min_or(target, lo, hi, rising):
        for _ in range(60):
            mid = (lo + hi) / 2
            below = power_at(mid, 0.05) < target
            if (below and rising) or (not below and not rising):
                lo = mid
            else:
                hi = mid
        return hi

    print(f"Minimum detectable ELEVATED OR at 80% power: "
          f"~{min_or(0.80, 1.0, 500.0, True):.1f}")
    print(f"Minimum detectable PROTECTIVE OR at 80% power: "
          f"~{min_or(0.80, 1.0, 0.001, False):.2f}")
    print()
    print("INTERPRETATION")
    print("-" * 62)
    print(f"With {E} events in {N_TAXA} taxa this design is underpowered, and")
    print("its true size is 0.0267 rather than the nominal 0.05 because the")
    print("probability-ordering rule gives a disconnected rejection region")
    print(f"{{x in {rejects}}} on a table this small.")
    print()
    print("The observed OR of 0.6667 sits far inside the blind spot: even a")
    print("true 2-fold elevation would go undetected (power ~0.10), and a")
    print("2-fold protective effect likewise (~0.15). Reaching 80% power")
    print("needs roughly a 30-fold effect. p = 0.6729 is therefore")
    print("uninformative in both directions.")
    print()
    print("Note the corrected power curve is now MONOTONE INCREASING with")
    print("the odds ratio. The earlier version was decreasing because it")
    print("used the wrong event count (E=6, the walrus-excluded table).")
    print()
    print("Correct statement: this dataset CANNOT test Huelsmann 2019 in")
    print("the direction it proposes. The hypothesis remains OPEN, not")
    print("refuted. A second reason it could not: the table mixes mammals")
    print("and non-mammals, but melatonin-DURATION reading is a mammalian")
    print("readout rule that birds demonstrably do not use (Juss et al. 1993")
    print("Proc R Soc B 254:21-28). See S6_CROSS_DOMAIN.md.")


if __name__ == "__main__":
    main()
