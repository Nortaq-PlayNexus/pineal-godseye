"""Test Huelsmann et al. 2019's hypothesis against the project phylogeny dataset.

Huelsmann (Sci Adv 5:eaaw6671) proposed that loss of melatonin synthesis and
receptor genes "may have been a precondition for the evolution of
unihemispheric sleep" in cetaceans.

This reads the shipped vertebrate_pineal_phylogeny.csv (the same file whose
errors were corrected 2026-10-17), defines a melatonin-pathway-intact flag,
and tests association with unihemispheric_sleep and thermal regime.

Bounded: 25 taxa, no allocation beyond a dict. Deliberately NOT a causal
claim -- a contingency table cannot establish causation. Reported as
association only, with the small-N caveat made explicit.
"""
import csv
import sys
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]


def main():
    rows = list(csv.DictReader(open(RAW, newline="", encoding="utf-8")))
    if len(rows) != 25:
        sys.exit(f"expected 25 taxa, got {len(rows)}")

    for r in rows:
        r["melatonin_intact"] = all(r[g] == "1" for g in GENES)

    intact = [r for r in rows if r["melatonin_intact"]]
    lost = [r for r in rows if not r["melatonin_intact"]]

    print(f"Taxa with intact melatonin pathway: {len(intact)}/25")
    print(f"Taxa with lost/broken pathway:     {len(lost)}/25")
    print()

    # Contingency: pathway state vs unihemispheric_sleep
    table = {"intact": {"1": 0, "0": 0}, "lost": {"1": 0, "0": 0}}
    for r in rows:
        key = "intact" if r["melatonin_intact"] else "lost"
        table[key][r["unihemispheric_sleep"]] += 1
    print("Pathway state  vs  unihemispheric_sleep (1 = present):")
    print(f"  intact, sleep=1: {table['intact']['1']}")
    print(f"  intact, sleep=0: {table['intact']['0']}")
    print(f"  lost,   sleep=1: {table['lost']['1']}")
    print(f"  lost,   sleep=0: {table['lost']['0']}")
    print()

    # Which taxa lost the pathway
    print("Taxa with lost/broken melatonin pathway:")
    for r in lost:
        broken = [g for g in GENES if r[g] != "1"]
        print(f"  {r['taxon']:<32} broken: {', '.join(broken)}")
    print()

    # Thermal regime distribution
    from collections import Counter
    regimes = Counter(r["thermal_metabolic_regime"] for r in rows)
    print("Thermal/metabolic regimes across all 25 taxa:")
    for reg, n in regimes.most_common():
        print(f"  {reg:<34} {n}")

    # Fisher exact test for the 2x2, if scipy is available
    try:
        from scipy.stats import fisher_exact
        odds, p = fisher_exact([[table["intact"]["1"], table["intact"]["0"]],
                                [table["lost"]["1"], table["lost"]["0"]]])
        print(f"\nFisher exact (pathway state vs unihemispheric sleep):")
        print(f"  odds ratio = {odds:.3f}, p = {p:.4f}")
        print("  NOTE: n=25, single 2x2, no correction for multiple testing.")
        print("  Association only -- cannot establish causation.")
    except ImportError:
        print("\nscipy not available; skipping Fisher exact.")


if __name__ == "__main__":
    main()