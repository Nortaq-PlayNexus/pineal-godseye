import csv
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"

print("=== scRNAseq_pineal_celltypes.csv ===")
rows = list(csv.DictReader(open(RAW / "scRNAseq_pineal_celltypes.csv",
                                newline="", encoding="utf-8")))
print(f"{'cell_type':<28}{'frac_total':<13}{'within':<9}{'sourced':<9}")
for r in rows:
    yn = "yes" if r["source"].strip() else "NO"
    print(f"{r['cell_type']:<28}{r['fraction_of_total_cells']:<13}"
          f"{r['fraction_within_lineage']:<9}{yn:<9}")
print()
print("Mays' reported shares: pinealocytes 90%, astrocytes 7%, microglia 1%,")
print("VLMC 2%, endothelial 0.1%  -> sum 1.001")

print()
print("=== vertebrate_pineal_phylogeny.csv ===")
rows2 = list(csv.DictReader(open(RAW / "vertebrate_pineal_phylogeny.csv",
                                 newline="", encoding="utf-8")))
intact = sum(1 for r in rows2 if r["gene_intact_all"] == "1")
lost = sum(1 for r in rows2 if r["gene_intact_all"] == "0")
unk = sum(1 for r in rows2 if r["gene_intact_all"] == "UNKNOWN")
nv = sum(1 for r in rows2 if "NOT VERIFIED" in r["source"])
print(f"rows: {len(rows2)}   intact: {intact}   lost: {lost}   unknown: {unk}")
print(f"rows explicitly marked NOT VERIFIED: {nv}")
print()
print("Previously shipped 5 binary morphology + 1 opsin + 2 free-text columns")
print("are GONE. Only gene-intactness remains, with per-row provenance.")
