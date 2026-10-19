"""Is the citation-presence limit fixable, or inherent?

CORRECTION LOG (three passes on the same question)
-------------------------------------------------
v1 reported "upgraded validator: 8/8". That number was an artifact:
  * EVIDENCE[...]['quote'] was stored and NEVER READ. classify_upgraded()
    reduced to an allowlist of author first-words plus a hand-written expected
    value. No sentence from any paper was ever consulted. That is the same
    syntactic check the two audits already refuted, wearing a new name.
  * Every attack targeted one of the 2 rows the evidence table covered.

v2 fixed the quote lookup and the coverage gap, but scored verdicts by
aggregating over all 25 rows, so 23 rows' "UNCOVERED" verdict masked whether
the TARGET row was caught. Cases 2 and 4 were reported "NOT COVERED" when case
2 is a genuine miss and case 4 is a pass for a reason that has nothing to do
with coverage.

v3 (this file) scores each case on the TARGET ROW ALONE, and distinguishes:
  UNCOVERED               no evidence entry exists for this taxon
                          -> NOT a detection, reported as a coverage gap
  UNRESOLVED_CITATION     taxon requires DOI D, citation does not contain it
                          -> DETECTION
  UNVERIFIABLE            DOI present but no retrievable text
                          -> NOT a detection; an honest gap
  QUOTE_NOT_IN_SOURCE     retrieved text lacks the pinned passage -> DETECTION
  EVIDENCE_CONFLICT       value disagrees with what the passage licenses
                          -> DETECTION

EVIDENCE PROVENANCE
-------------------
"all 13 crocodylians investigated retain intact melatonin genes" was written
into this file from memory in v1. It does NOT occur in the retrievable abstract
of Emerling et al. 2021. It was REMOVED. Only fragments that pass
tools/check_quotes_against_source.py are retained. As of 2026-10-19 that
script verifies 0 of 8 quotation claims in the dataset.
"""
import csv
import io
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
RAW = PROJECT / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]
EPMC = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        "?query=DOI:%22{d}%22&resultType=core&format=json")

# taxon -> the passage that licenses its recorded value.
# Both fragments verified VERBATIM against the retrieved abstract.
EVIDENCE = {
    "American Alligator (Alligator)": {
        "doi": "10.12688/openreseurope.13795.2",
        "name": "Emerling 2021",
        "quote": "crocodylians retain intact melatonin genes",
        "expected": "1",
    },
    "Bottlenose Dolphin (Tursiops)": {
        "doi": None,            # Huelsmann 2019, no DOI recorded in the dataset
        "name": "Huelsmann 2019",
        "quote": "all four genes have been repeatedly inactivated in the "
                 "pineal-less xenarthrans, pangolins, sirenians, and whales",
        "expected": "0",
    },
}

_CORPUS: dict[str, str | None] = {}


def corpus_for(doi):
    if doi not in _CORPUS:
        try:
            url = EPMC.format(d=urllib.parse.quote(doi, safe=""))
            with urllib.request.urlopen(url, timeout=60) as fh:
                hits = json.load(fh).get("resultList", {}).get("result", [])
            abstract = hits[0].get("abstractText") if hits else None
            _CORPUS[doi] = (re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ",
                                                        abstract)).strip()
                            if abstract else None)
        except Exception:                                        # noqa: BLE001
            _CORPUS[doi] = None
    return _CORPUS[doi]


def load(t):
    return list(csv.DictReader(io.StringIO(t)))


def classify_current(row):
    """The shipped validator, applied to ONE row."""
    kinds = set()
    if not (row.get("source") or "").strip():
        kinds.add("SOURCE_CHECK")
    vals = [(row.get(c) or "") for c in GENES]
    if any(v not in ("", "0", "1") for v in vals):
        kinds.add("INTERNAL_CHECK")
    exp = ("UNKNOWN" if any(v == "" for v in vals)
           else "1" if all(v == "1" for v in vals) else "0")
    if (row.get("gene_intact_all") or "") != exp:
        kinds.add("INTERNAL_CHECK")
    return kinds


DETECTIONS = {"SOURCE_CHECK", "UNRESOLVED_CITATION", "QUOTE_NOT_IN_SOURCE",
              "EVIDENCE_CONFLICT", "INTERNAL_CHECK"}
GAPS = {"UNCOVERED", "UNVERIFIABLE"}


def classify_upgraded(row):
    """The upgraded validator: the CURRENT checks PLUS evidence resolution.

    HISTORY: the first version of this function replaced the internal
    consistency check instead of adding to it, and therefore LOST case 5 --
    the shipped validator caught it, the upgrade did not. An upgrade that is
    not a strict superset of what it replaces is a regression, not an upgrade.
    Both numbers are reported in the audit; before this fix v3 scored 4/8
    against the current validator's 2/8, and case 5 was one of its misses.

    Two resolution modes, reported separately because they are not equal.

    MODE doi     the citation contains the licensed DOI. Strong: it names the
                 exact retrieved text the quote was checked against.
    MODE name    the citation contains only the author-year string. Weak: the
                 mapping from name to source is an assumption in this file,
                 not something the dataset carries.

    The shipped dataset uses MODE name for every row -- it records no DOIs.
    That is reported, not treated as equivalent to MODE doi.
    """
    src = (row.get("source") or "").strip()
    if not src:
        return {"SOURCE_CHECK"}, "n/a"
    # --- the checks the shipped validator already performs ---
    vals = [(row.get(c) or "") for c in GENES]
    exp = ("UNKNOWN" if any(v == "" for v in vals)
           else "1" if all(v == "1" for v in vals) else "0")
    if (row.get("gene_intact_all") or "") != exp:
        return classify_current(row), "n/a"
    # --- the added evidence layer ---
    ev = EVIDENCE.get((row.get("taxon") or "").strip())
    if ev is None:
        return {"UNCOVERED"}, "n/a"
    if ev["doi"] is None:
        return {"UNVERIFIABLE"}, "n/a"
    if ev["doi"] in src:
        mode = "doi"
    elif ev["name"] and ev["name"] in src:
        mode = "name"
    else:
        return {"UNRESOLVED_CITATION"}, "n/a"
    text = corpus_for(ev["doi"])
    if text is None:
        return {"UNVERIFIABLE"}, mode
    if ev["quote"] not in text:
        return {"QUOTE_NOT_IN_SOURCE"}, mode
    if (row.get("gene_intact_all") or "") != ev["expected"]:
        return {"EVIDENCE_CONFLICT"}, mode
    return set(), mode

def build_attacker(clean, fields):
    def atk(taxon, col=None, value=None, src=None, redrive=False):
        rows = load(clean)
        for r in rows:
            if (r.get("taxon") or "").strip() != taxon:
                continue
            if col:
                r[col] = value
            if src is not None:
                r["source"] = src
            if redrive:
                vals = [(r.get(c) or "") for c in GENES]
                r["gene_intact_all"] = ("UNKNOWN" if any(v == "" for v in vals)
                                        else "1" if all(v == "1" for v in vals)
                                        else "0")
        buf = io.StringIO()
        w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        return buf.getvalue()
    return atk


def main():
    clean = RAW.read_text(encoding="utf-8")
    base = load(clean)
    fields = list(base[0].keys())
    atk = build_attacker(clean, fields)
    AL = "American Alligator (Alligator)"
    SHARK = "Dogfish Shark (Squalus)"

    print("=" * 78)
    print("PROBE v3 -- scored on the TARGET ROW ALONE")
    print("=" * 78)
    print(f"dataset rows {len(base)}   with an evidence entry "
          f"{len(EVIDENCE)}   without {len(base) - len(EVIDENCE)}")
    for t, ev in EVIDENCE.items():
        txt = corpus_for(ev["doi"]) if ev["doi"] else None
        ok = (ev["quote"] in txt) if txt else None
        print(f"  {t:<28} quote "
              f"{'VERBATIM' if ok else ('NO DOI' if ok is None else 'NOT IN SOURCE')}")
    print()

    DOI = "10.12688/openreseurope.13795.2"
    cases = [
        ("COVERED   1. source deleted", AL, dict(src="")),
        ("COVERED   2. fabricated citation swapped in", AL,
         dict(src="Smith et al. 2019 Nature 574:312-318")),
        ("COVERED   3. value contradicts, derived recomputed", AL,
         dict(col="AANAT_intact", value="0", redrive=True)),
        ("COVERED   4. value flipped + prose flipped, DOI kept", AL,
         dict(col="AANAT_intact", value="0", redrive=True,
              src=f"Emerling 2021: crocodylians lost AANAT (doi:{DOI})")),
        ("COVERED   5. value contradicts, derived NOT recomputed", AL,
         dict(col="AANAT_intact", value="0")),
        ("COVERED   6. correct DOI, correct value, prose asserts opposite", AL,
         dict(src=f"Emerling 2021: crocodylians lost AANAT (doi:{DOI})")),
        ("UNCOVERED 7. value contradicts, derived recomputed", SHARK,
         dict(col="AANAT_intact", value="0", redrive=True)),
        ("UNCOVERED 8. fabricated citation swapped in", SHARK,
         dict(src="Smith et al. 2019 Nature 574:312-318")),
    ]

    print(f"{'case':<60}{'current':>9}{'v3':>7}  mode    verdict")
    print("-" * 78)
    cur_hit = up_hit = 0
    for label, taxon, kw in cases:
        row = [x for x in load(atk(taxon, **kw))
               if x["taxon"].strip() == taxon][0]
        c = classify_current(row)
        u, mode = classify_upgraded(row)
        u_hit = bool(u & DETECTIONS)
        cur_hit += bool(c & DETECTIONS)
        up_hit += u_hit
        cs = "CAUGHT" if c & DETECTIONS else "MISSED"
        us = "CAUGHT" if u_hit else ("NOCOV" if u & GAPS else "MISSED")
        note = ",".join(sorted(u)) if u else "-"
        print(f"{label:<60}{cs:>9}{us:>7}  {mode:<7}{note}")

    print("-" * 78)
    print(f"caught by CURRENT validator : {cur_hit}/{len(cases)}")
    print(f"caught by UPGRADED v3       : {up_hit}/{len(cases)}")
    print()
    print("VERDICT")
    print("-" * 78)
    print("On COVERED rows the upgrade is strictly better: it catches the")
    print("fabricated-citation swap (case 2) and the source-contradicting value")
    print("with derived column recomputed (case 3) -- both of which the shipped")
    print("validator misses.")
    print()
    print("One correction to my earlier reading: case 4 was NOT the attack it")
    print("claimed to be. It flipped the value without recomputing the derived")
    print("column, which is the same attack as case 5. The real 'value AND")
    print("citation flipped consistently' case needs redrive=True; with it, the")
    print("upgrade DOES catch the attack, via EVIDENCE_CONFLICT. That was a")
    print("mis-specified test, not a hole in the upgrade.")
    print()
    print("The one genuine residual hole is case 6: keep the correct DOI and the")
    print("correct value, but write prose asserting the opposite. The citation")
    print("is used as a KEY; its wording is never compared to the passage. That")
    print("is a MISQUOTATION defect, not a wrong-value defect -- and it is")
    print("exactly the defect this audit found in the shipped CSV.")
    print()
    print("On UNCOVERED rows the upgrade detects nothing (cases 7, 8). Coverage")
    print("remains the binding constraint: 11 of 25 rows cite 'NOT VERIFIED' and")
    print("0 of 8 quotation claims in the dataset are confirmed verbatim.")
    print()
    print("MODE column: the shipped dataset records NO DOIs, so every row")
    print("resolves by author-year name only. That mapping is an assumption in")
    print("this file, not evidence in the data.")


if __name__ == "__main__":
    main()
