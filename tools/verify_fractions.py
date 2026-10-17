"""Verify the scRNAseq lineage fractions against Mays et al. 2018.

Mays values below were read from the PLOS ONE full text (PMC6197868):
  - "Pinealocytes accounted for 90% of the profiled cells"
  - "Astrocytes accounted for 7% of profiled cells"
  - "alpha-, beta-, and gamma-Astrocytes accounted for 85%, 7%, and 8%"
  - "Microglia accounted for 1% of profiled cells"
  - "alpha- and beta-Microglia accounted for 64% and 36%"
  - "VLMCs accounted for 2% of profiled cells"
  - "Endothelial cells accounted for 0.1% of profiled cells"
"""
import csv
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "scRNAseq_pineal_celltypes.csv"

# Lineage -> (mays_total_fraction, dataset_expected)
MAYS = {
    "Pinealocyte": 0.90,
    "Astrocyte": 0.07,
    "Microglia": 0.01,
    "Vascular/Leptomeningeal": 0.02,
    "Endothelial": 0.001,
}


def main():
    rows = list(csv.DictReader(open(RAW, newline="", encoding="utf-8")))

    total = sum(float(r["fraction_of_total_cells"]) for r in rows)
    print(f"TOTAL across all {len(rows)} rows = {total:.6f}")
    print()

    print("LINEAGE SHARE OF TOTAL CELLS")
    print(f"{'lineage':<28}{'dataset':>10}{'Mays':>10}   verdict")
    print("-" * 58)
    failures = 0
    for lin, m in MAYS.items():
        d = sum(float(r["fraction_of_total_cells"]) for r in rows
                if r["major_lineage"] == lin)
        ok = abs(d - m) < 0.005
        if not ok:
            failures += 1
        print(f"{lin:<28}{d:>10.3f}{m:>10.3f}   {'OK' if ok else 'MISMATCH'}")
    print()

    print("WITHIN-LINEAGE SPLIT: astrocytes")
    print("  Mays 2018: alpha 85% / beta 7% / gamma 8%")
    astro = [r for r in rows if r["major_lineage"] == "Astrocyte"]
    s = sum(float(r["fraction_of_total_cells"]) for r in astro)
    for r in astro:
        v = float(r["fraction_of_total_cells"]) / s * 100
        print(f"    {r['cell_type']:<24}{v:>7.1f}%")
    print()

    nonpine = sum(float(r["fraction_of_total_cells"]) for r in rows
                  if r["major_lineage"] != "Pinealocyte")
    print(f"Sum of the 7 NON-pinealocyte rows = {nonpine:.6f}")
    print(f"Sum of the 2 pinealocyte rows      = {1 - nonpine:.6f}   "
          f"(Mays: 0.900)")
    print()

    if failures:
        print(f"RESULT: {failures} lineage fraction(s) DISAGREE with Mays 2018.")
        print()
        print("The pinealocyte block matches the source exactly. The remaining")
        print("seven rows sum to exactly the residual complement (0.100) and")
        print("all nine sum to exactly 1.000 -- i.e. they were distributed to")
        print("satisfy the sum-to-one validator, not transcribed from the")
        print("paper. The build-time validator PASSES this construction;")
        print("it does not detect it.")
    else:
        print("RESULT: all lineage fractions match Mays 2018.")


if __name__ == "__main__":
    main()
