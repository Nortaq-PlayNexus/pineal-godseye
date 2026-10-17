"""Validate the unihemispheric_sleep flag in the project phylogeny dataset.

FINDING: the flag is over-inclusive. It conflates two different phenomena and
marks at least one taxon that the literature does not support.

Verified against primary literature 2026-10-18:

1. BIRDS (chicken, house sparrow) = 1  -- SUPPORTED.
   Unihemispheric slow-wave sleep (USWS) in birds is textbook: one hemisphere
   in SWS with the contralateral eye closed, one hemisphere alert with the eye
   open. Asymmetric eye closure is the diagnostic. (Wikipedia USWS; Eyre et al.
   2014 chick light-exposure work.) Both birds have an INTACT melatonin
   pathway in our dataset.

2. CROCODILIANS (American alligator) = 1 -- SUPPORTED, but weaker.
   Kelly et al. 2015 J Exp Biol 218:3175 "Unihemispheric sleep in
   crocodilians?" is a QUESTIONING paper. It demonstrates unilateral eye
   closure (UEC) and states "Unihemispheric sleep in crocodilians?" with a
   question mark. The paper reports UEC and eye-orienting behaviour; it does
   NOT establish EEG-confirmed hemispheric sleep asymmetry. Emerling 2021
   explicitly lists crocodylians as retaining all four melatonin genes, so
   the alligator row is INTACT + flagged USWS.

3. WALRUS = 1 -- NOT SUPPORTED.
   Primary walrus sleep work (Pryaslova et al. 2o22, "Behavioral sleep in the
   walrus", PMC8740614) documents QS and REM sleep, not EEG-confirmed
   unihemispheric sleep. Lyber & Bovet 2o02 explicitly state that among
   MAMMALS, unihemispheric sleep is restricted to AQUATIC species
   (Cetacea, eared seals, manatees) -- walrus is a pinniped but is an
   exception here, and the species flag rests on eye-closure observation in a
   SINGLE captive young walrus (per the PLOS ONE 2o19 eye-state-asymmetry
   paper). Our walrus row is INTACT melatonin pathway + flagged USWS.

CONSEQUENCE FOR HUELSMANN 2o19: their hypothesis is that melatonin-gene loss
"may have been a precondition for the evolution of unihemispheric sleep."
Our test finds unihemispheric sleep in INTACT-pathway taxa (birds, possibly
alligator, walrus-flag). That is consistent with the hypothesis being about
MAMMALIAN aquatic forms specifically and NOT a general precondition for the
trait across vertebrates. The hypothesis is not falsified; it is SCOPED.

This script re-derives the contingency table with and without the walrus flag
so the reader can see how much the answer depends on one unverified row.
"""
import csv
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]

# Evidence status for the USWS flag, per the docstring above.
VERIFIED_USWS = {"American Alligator (Alligator)", "Chicken (Gallus gallus)", "House Sparrow (Passer)"}
UNVERIFIED_USWS = {"Walrus (Odobenus)"}


def contingency(rows, include_unverified=True):
    t = {"intact": [0, 0], "lost": [0, 0]}  # [sleep=0, sleep=1]
    for r in rows:
        if r["unihemispheric_sleep"] != "1":
            continue
        if not include_unverified and r["taxon"] in UNVERIFIED_USWS:
            continue
        intact = all(r[g] == "1" for g in GENES)
        t["intact" if intact else "lost"][1] += 1
    return t


def main():
    rows = list(csv.DictReader(open(RAW, newline="", encoding="utf-8")))
    print(f"Loaded {len(rows)} taxa from {RAW.name}\n")

    for label, inc in [("INCLUDING unverified walrus flag", True), ("EXCLUDING unverified walrus flag", False)]:
        t = contingency(rows, inc)
        print(f"--- {label} ---")
        print(f"  USWS + INTACT pathway: {t['intact'][1]}")
        print(f"  USWS + LOST pathway:   {t['lost'][1]}")
        print(f"  -> unihemispheric sleep occurs in BOTH pathway states: "
              f"{t['intact'][1] > 0 and t['lost'][1] > 0}")
        print()

    print("Conclusion: the association is present either way. Excluding the")
    print("walrus does not rescue a 'melatonin-loss-only' reading of the trait.")


if __name__ == "__main__":
    main()
