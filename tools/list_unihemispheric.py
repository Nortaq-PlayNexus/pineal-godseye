"""List the taxa flagged unihemispheric_sleep=1 and their melatonin state.

The key question for Huelsmann's hypothesis: is unihemispheric sleep exclusive
to gene-loss lineages? If intact-pathway taxa also show it, the "melatonin loss
is a precondition for unihemispheric sleep" framing is weaker than the paper
implies, because sleep appears without loss.
"""
import csv
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "vertebrate_pineal_phylogeny.csv"
GENES = ["AANAT_intact", "ASMT_intact", "MTNR1A_intact", "MTNR1B_intact"]

rows = list(csv.DictReader(open(RAW, newline="", encoding="utf-8")))

print("Taxa flagged unihemispheric_sleep=1:")
for r in rows:
    if r["unihemispheric_sleep"] == "1":
        intact = all(r[g] == "1" for g in GENES)
        state = "INTACT" if intact else "LOST"
        broken = [g.replace("_intact", "") for g in GENES if r[g] != "1"]
        b = f" (broken: {', '.join(broken)})" if broken else ""
        print(f"  [{state}] {r['taxon']:<30} {r['clade']}{b}")

print()
print("Reconciliation: the CSV column is a hand-assigned literature flag, not")
print("a measured trait. Several rows are questionable and are noted below.")
print("The association test is therefore only as good as those flags.")
