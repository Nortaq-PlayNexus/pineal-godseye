"""Prove the build-time validation catches the exact regression it was written for.

Re-introduces the original bug (AANAT_intact=2 for zebrafish, a paralog count
written into a binary flag) into a copy of the shipped CSV and confirms the
binary-integrity check fires on it. Exits 0 only if the check catches it.
"""
import csv
import io
import sys
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
BINARY_COLS = [
    "parietal_eye", "pineal_organ", "direct_photosensitive",
    "intrinsic_oscillator", "AANAT_intact", "ASMT_intact",
    "MTNR1A_intact", "MTNR1B_intact", "unihemispheric_sleep",
]


def scan(text):
    hits = []
    for row in csv.DictReader(io.StringIO(text)):
        for c in BINARY_COLS:
            if row[c] not in ("0", "1"):
                hits.append(f"{row['taxon']}.{c}={row[c]!r}")
    return hits


def main():
    clean = (RAW / "vertebrate_pineal_phylogeny.csv").read_text(encoding="utf-8")
    assert scan(clean) == [], "shipped CSV unexpectedly fails its own check"

    # Re-introduce the original bug: teleost paralog count (2) in a binary flag.
    # The CSV is unquoted, so match the bare-field form exactly.
    corrupted = clean.replace(
        "Zebrafish (Danio rerio),Actinopterygii,Fishes,0,1,1,1,4,1,1,1,1,0",
        "Zebrafish (Danio rerio),Actinopterygii,Fishes,0,1,1,1,4,2,1,1,1,0",
        1,
    )
    assert corrupted != clean, "replacement failed to apply"

    caught = scan(corrupted)
    if not caught:
        print("FAIL: validation did not catch the re-introduced bug")
        sys.exit(1)
    if not any(h.endswith("AANAT_intact='2'") for h in caught):
        print(f"FAIL: caught {caught} but not the zebrafish AANAT_intact=2 case")
        sys.exit(1)

    print("PASS: clean CSV scans clean; corrupted CSV caught as:", caught)
    sys.exit(0)


if __name__ == "__main__":
    main()