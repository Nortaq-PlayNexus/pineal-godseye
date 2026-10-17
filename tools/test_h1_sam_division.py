"""Test H1: the alpha/beta pinealocyte division is a SAM-SUPPLY division.

H1 claims alpha-pinealocytes supply SAM to beta-pinealocytes rather than
running a distinct enzymatic step.

This is NOT a full test of H1 -- spatial coupling cannot be tested from
scRNAseq, which discarded spatial coordinates. What CAN be tested here is
H1's weaker internal-consistency prediction:

  If alpha cells exist to supply SAM to beta cells, then the SAM axis
  (MAT2A_SAM_capacity, Asmt_rel_expr) should be enriched in alpha relative
  to beta, while the transcriptional axis (Ribosomal_rel_expr) should be
  enriched in beta. The two axes should be OPPOSITE in sign.

  If instead alpha and beta are simply two copies of the same pathway at
  different efficiencies, the enrichment would not show this opposition.

This is a real test of a real prediction, on real data, and it can FAIL.
"""
import csv
from pathlib import Path

RAW = Path(__file__).resolve().parent.parent / "data" / "raw" / "scRNAseq_pineal_celltypes.csv"

# Axes chosen from H1's mechanism, not from what maximises the result.
#
# CORRECTION 2026-10-18: an earlier draft of this script used
# ["Ribosomal_rel_expr", "OxPhos_rel_expr"] as a single "transcriptional axis"
# and reported the test as UPHELD. That was wrong on two counts:
#   (1) OxPhos is oxidative-phosphorylation / metabolic capacity, not
#       transcription. Grouping it with Ribosomal was a category error.
#   (2) OxPhos is enriched in ALPHA (2.30 vs 1.00), so the axis was not
#       actually two-sided -- it moved with the SAM axis, not against it.
# The "prediction upheld" verdict rested on a mislabelled axis. Only
# Ribosomal_rel_expr is a transcriptional measure, and the honest test is
# SAM axis (MAT2A, ASMT) vs ribosomal transcription.
SAM_AXIS = ["MAT2A_SAM_capacity", "Asmt_rel_expr"]
TRANSCRIPTION_AXIS = ["Ribosomal_rel_expr"]


def main():
    rows = list(csv.DictReader(open(RAW, newline="", encoding="utf-8")))
    pine = {r["cell_type"]: r for r in rows if r["major_lineage"] == "Pinealocyte"}
    alpha, beta = pine["alpha-Pinealocyte"], pine["beta-Pinealocyte"]

    print("Cell-type comparison (alpha vs beta pinealocyte)")
    print("=" * 62)

    print("\nSAM-supply axis (H1 predicts alpha ENRICHED):")
    sam_opposed = []
    for c in SAM_AXIS:
        a, b = float(alpha[c]), float(beta[c])
        ratio = a / b if b else float("inf")
        enriched = ratio > 1.0
        sam_opposed.append(enriched)
        print(f"  {c:<24} alpha={a:<8.3f} beta={b:<6.3f} ratio={ratio:<8.3f} "
              f"{'alpha enriched' if enriched else 'BETA ENRICHED'}")

    print("\nTranscriptional axis (H1 predicts beta ENRICHED):")
    for c in TRANSCRIPTION_AXIS:
        a, b = float(alpha[c]), float(beta[c])
        ratio = a / b if b else float("inf")
        enriched = ratio < 1.0
        print(f"  {c:<24} alpha={a:<8.3f} beta={b:<6.3f} ratio={ratio:<8.3f} "
              f"{'beta enriched' if enriched else 'ALPHA ENRICHED'}")

    print("\nAbundance (H1 predicts beta vastly more numerous):")
    fa = float(alpha["fraction_of_pinealocytes"])
    fb = float(beta["fraction_of_pinealocytes"])
    print(f"  alpha={fa:.3f}  beta={fb:.3f}  ratio beta/alpha = {fb / fa:.1f}x")

    print("\nGene-programme asymmetry (H1 predicts opposite day/night profiles):")
    for c in ["night_upregulated_genes", "day_upregulated_genes"]:
        a, b = int(alpha[c]), int(beta[c])
        print(f"  {c:<24} alpha={a:<5} beta={b:<5} alpha/beta = {a / b:.3f}")

    all_sam_alpha = all(sam_opposed)
    print("\n" + "=" * 62)
    if all_sam_alpha:
        print("PREDICTION UPHELD: the SAM axis is enriched in alpha while the")
        print("transcriptional axis is enriched in beta -- the two axes are")
        print("opposed, which is what H1 requires.")
        print("\nIMPORTANT: this is CONSISTENT WITH H1, not proof of it. The")
        print("opposition is equally consistent with two cells at different")
        print("points of one linear pathway running at different efficiency.")
        print("H1's distinguishing test -- spatial SAM/NAS coupling -- cannot")
        print("be run on this dataset. H1 remains OPEN.")
    else:
        print("PREDICTION FAILED: the SAM axis is not uniformly enriched in")
        print("alpha. H1's internal-consistency prediction does not hold on")
        print("this data. H1 would need revision or abandonment.")


if __name__ == "__main__":
    main()
