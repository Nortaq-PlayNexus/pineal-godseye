"""Check every quoted fragment in the dataset against the ACTUAL source text.

WHY THIS EXISTS
---------------
Two earlier audits (AUDIT_2026-10-18, AUDIT_TEST7_2026-10-19) established that
the shipped provenance validator cannot detect a fabricated citation. That left
an obvious question: do the citations the dataset ACTUALLY ships survive
contact with their sources?

v1 of this script hardcoded a hand-picked list of fragments. That is a defect:
a hand-picked list is exactly what missed the rows, and it can be curated to
whatever the author already believes. v2 EXTRACTS every single-quoted fragment
from the source column automatically, so a new row cannot escape the check.

VERDICT VOCABULARY
------------------
  VERBATIM   the fragment occurs exactly in retrieved source text
  NOT FOUND  it does not occur -- it is not a quotation of that text
  NO SOURCE  no retrievable text for that citation (e.g. no DOI recorded)

'NOT FOUND' does NOT prove the underlying VALUE is wrong. It proves the string
is not a quotation of the text that could be retrieved. A passage in the paper
body would not be seen, because Europe PMC serves no full text for these
records (no PMCID, isOpenAccess=N). That caveat is printed with the results.
"""
import csv
import io
import json
import re
import sys
import urllib.parse
import urllib.request

EPMC = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        "?query=DOI:%22{d}%22&resultType=core&format=json")

# Which citation each source string refers to, and whether we can retrieve it.
# A citation with no DOI here gets NO SOURCE -- honestly, not optimistically.
RETRIEVAL = {
    "Emerling 2021": "10.12688/openreseurope.13795.2",
}

_CACHE: dict[str, str | None] = {}


def fetch(doi):
    if doi in _CACHE:
        return _CACHE[doi]
    url = EPMC.format(d=urllib.parse.quote(doi, safe=""))
    try:
        with urllib.request.urlopen(url, timeout=60) as fh:
            data = json.load(fh)
        hits = data.get("resultList", {}).get("result", [])
        abstract = hits[0].get("abstractText") if hits else None
    except Exception as exc:                                  # noqa: BLE001
        print(f"  ! retrieval failed for {doi}: {exc}")
        abstract = None
    text = (re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", abstract)).strip()
            if abstract else None)
    _CACHE[doi] = text
    return text


def normalise(s):
    return re.sub(r"\s+", " ", s).strip()


def fragments(source):
    """Every single-quoted fragment -- the dataset's claims of VERBATIM."""
    return re.findall(r"'([^']+)'", source or "")


def verdict(frag, text):
    if text is None:
        return "NO SOURCE"
    f = normalise(frag)
    if f in text:
        return "VERBATIM"
    # an elided quotation ('a ... b') can only match piecewise
    parts = [normalise(p) for p in f.split("...") if normalise(p)]
    if parts and all(p in text for p in parts):
        return "VERBATIM (elided)"
    return "NOT FOUND"


def main():
    raw = open("data/raw/vertebrate_pineal_phylogeny.csv", encoding="utf-8").read()
    rows = list(csv.DictReader(io.StringIO(raw)))

    texts = {name: fetch(doi) for name, doi in RETRIEVAL.items()}
    print("=" * 78)
    print("QUOTATION AUDIT -- every single-quoted fragment in the source column")
    print("=" * 78)
    for name, t in texts.items():
        print(f"  {name:<18} "
              f"{'retrieved ' + str(len(t)) + ' chars' if t else 'NO TEXT RETRIEVED'}")
    print()

    total = ok = nf = nos = 0
    findings = []
    for r in rows:
        src = r.get("source") or ""
        frags = fragments(src)
        if not frags:
            continue
        match = next((n for n in RETRIEVAL if n in src), None)
        text = texts.get(match) if match else None
        for frag in frags:
            v = verdict(frag, text)
            total += 1
            ok += v.startswith("VERBATIM")
            nf += v == "NOT FOUND"
            nos += v == "NO SOURCE"
            findings.append((r["taxon"], v, normalise(frag)))

    print(f"{'verdict':<18}{'taxon':<28}fragment")
    print("-" * 78)
    for taxon, v, frag in findings:
        shown = frag if len(frag) <= 38 else frag[:35] + "..."
        print(f"{v:<18}{taxon:<28}{shown}")

    print("-" * 78)
    print(f"quotation claims found : {total}")
    print(f"  VERBATIM in source   : {ok}")
    print(f"  NOT FOUND            : {nf}")
    print(f"  NO SOURCE (no DOI)   : {nos}")
    print()

    # rows asserting a value with no retrievable citation at all
    unbacked = [r["taxon"] for r in rows
                if (r.get("source") or "").strip() == "NOT VERIFIED"]
    print(f"rows asserting a value with NO source at all : {len(unbacked)}/{len(rows)}")
    print()

    print("CAVEAT -- read before acting on NOT FOUND")
    print("-" * 78)
    print("A NOT FOUND verdict means the string is not a quotation of the")
    print("RETRIEVABLE text. Europe PMC serves only the abstract for these")
    print("records (no PMCID, isOpenAccess=N), so a genuine sentence from the")
    print("paper body would also read NOT FOUND.")
    print()
    print("What IS established regardless of that caveat:")
    print("  1. The dataset presents these strings inside quotation marks, so")
    print("     it asserts they are verbatim. Nothing in the build verified")
    print("     that assertion, and the audit could not confirm it.")
    print("  2. The same audit reproduces one of them: writing a fragment")
    print("     from memory produced a string the source does not contain.")
    print("     The defect is in the process that produces citations, not in")
    print("     one stale row.")
    sys.exit(1 if (nf or nos) else 0)


if __name__ == "__main__":
    main()