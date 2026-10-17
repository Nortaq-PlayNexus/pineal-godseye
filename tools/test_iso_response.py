"""Do alpha and beta pinealocytes respond differently to adrenergic stimulation?

REWRITTEN after the first version was discarded. That version joined vehicle
and isoproterenol on cell barcodes and found 5 cells in common, then reported
percent-change figures derived from them. Those numbers were noise: 1 "alpha"
cell cannot support any comparison. The bug was in the abort condition, which
triggered on zero shared barcodes rather than on too few cells to mean anything.

THE DESIGN IS UNPAIRED. These are different animals -- the two samples share
essentially no barcodes (5 of 1,576 / 1,996), so cells cannot be matched
between conditions. The correct comparison is between GROUPS of cells, with a
test that accounts for the fact that thousands of comparisons are being made.

WHAT IS ACTUALLY TESTABLE HERE, AND WHAT IS NOT
Isoproterenol-induced changes in this dataset are between-animal. So a
difference in alpha-vs-beta response CAN be estimated, but it is confounded
with animal-to-animal variation and cannot be attributed to the drug with the
same confidence as a within-animal design. That limitation is inherent to the
public data, not to this script.

Mays et al. 2018 report that drug-induced changes concentrate in pinealocytes
(99%) but do not report whether alpha and beta respond DIFFERENTLY. That is
the gap this fills.
"""
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
NC = ROOT / "data" / "geo" / "normcounts"
FILES = {
    "vehicle": "GSM3188351_Vehicle_ExtendedRef_filtered_normcounts.txt.gz",
    "iso": "GSM3188352_Isoproterenol_ExtendedRef_filtered_normcounts.txt.gz",
}
ALPHA_FRAC = 0.05
MIN_EXPR = 0.1          # CP10K detection floor
MIN_CELLS_PER_GROUP = 30  # below this, no test is reported at all
FDR = 0.05


def load(fname):
    path = NC / fname
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = [h.strip('"').split("_", 1)[1] for h in header[2:]]
        syms, rows = [], []
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            syms.append(p[1].strip('"'))
            rows.append(p[2:])
    mat = np.array(rows, dtype=np.float32).T
    return barcodes, syms, mat


def welch(a, b):
    """Welch's t. Returns (t, dof). Normal approximation for large n."""
    na, nb = len(a), len(b)
    ma, mb = a.mean(), b.mean()
    va, vb = a.var(ddof=1), b.var(ddof=1)
    se = np.sqrt(va / na + vb / nb)
    if se == 0:
        return 0.0, 1.0
    t = (ma - mb) / se
    dof = (va / na + vb / nb) ** 2 / (
        (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    return float(t), float(dof)


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        val = min(prev, p[i] * n / (rank + 1))
        adj[i] = val
        prev = val
    return adj


def norm_p_from_t(t, dof):
    """Two-sided p from t using a normal approximation (dof > 50 here)."""
    from math import erfc, sqrt
    return erfc(abs(t) / sqrt(2.0))


def split_groups(syms, mat):
    idx = {s: i for i, s in enumerate(syms)}
    asmt = mat[:, idx["Asmt"]]
    pineal = np.where(asmt > MIN_EXPR)[0]
    order = pineal[np.argsort(asmt[pineal])[::-1]]
    k = max(1, int(round(ALPHA_FRAC * len(order))))
    alpha = order[:k]
    beta = order[k:]
    return idx, pineal, alpha, beta


def main():
    print("Isoproterenol response: alpha vs beta pinealocytes (UNPAIRED)\n")
    load = {}
    for k, f in FILES.items():
        load[k] = load_fn(f)
        print(f"  {k:<9} {load[k][2].shape[0]:>5} cells x "
              f"{load[k][2].shape[1]} genes")

    bv, bi = load["vehicle"][0], load["iso"][0]
    shared = set(bv) & set(bi)
    print(f"\n  shared barcodes between conditions: {len(shared)}")
    print("  -> conditions are DIFFERENT ANIMALS. Unpaired comparison only.")

    syms = load["vehicle"][1]
    assert syms == load["iso"][1], "gene lists differ between samples"
    n_genes = len(syms)
    gi = {s: i for i, s in enumerate(syms)}

    groups = {}
    for cond in ("vehicle", "iso"):
        _, _, mat = load[cond]
        _, pineal, alpha, beta = split_groups(syms, mat)
        groups[cond] = dict(pineal=pineal, alpha=alpha, beta=beta, mat=mat)
        print(f"\n  {cond}: pinealocytes {len(pineal)}  alpha {len(alpha)}  "
              f"beta {len(beta)}")

    minn = min(len(groups[c]["alpha"]) for c in groups)
    if minn < MIN_CELLS_PER_GROUP:
        print(f"\nABORT: smallest alpha group has {minn} cells, below the "
              f"{MIN_CELLS_PER_GROUP}-cell floor. No test reported.")
        return

    print(f"\nAlpha cell counts pass the {MIN_CELLS_PER_GROUP}-cell floor.")

    # Per-gene vehicle->iso change, then is alpha's change different from beta's?
    print(f"\nGENE-LEVEL: does alpha's drug response differ from beta's?")
    print(f"  Welch t on per-cell expression, vehicle vs iso, within group,")
    print(f"  then testing whether the alpha and beta t-statistics differ is not")
    print(f"  the right test -- instead we test alpha-vs-beta directly in each")
    print(f"  condition and look for genes where the SIGN of the drug effect")
    print(f"  differs between subtypes.\n")

    # Correct framing: for each gene compute drug effect in alpha and in beta,
    # then test whether the two effects differ (difference-in-differences).
    eff = {}
    for g in ("alpha", "beta"):
        a_rows = groups["vehicle"][g]
        b_rows = groups["iso"][g]
        V = groups["vehicle"]["mat"][a_rows]
        I = groups["iso"]["mat"][b_rows]
        eff[g] = I.mean(axis=0) - V.mean(axis=0)

    # Vectorised per-gene p for each group.
    def per_gene_p(mat_iso, mat_veh):
        a = mat_iso
        b = mat_veh
        na, nb = len(a), len(b)
        se = np.sqrt(a.var(axis=0, ddof=1) / na + b.var(axis=0, ddof=1) / nb)
        se[se == 0] = np.inf
        t = (a.mean(axis=0) - b.mean(axis=0)) / se
        from math import erfc, sqrt
        return np.array([erfc(abs(x) / sqrt(2.0)) for x in t]), t

    pa, ta = per_gene_p(groups["iso"]["mat"][groups["iso"]["alpha"]],
                        groups["vehicle"]["mat"][groups["vehicle"]["alpha"]])
    pb, tb = per_gene_p(groups["iso"]["mat"][groups["iso"]["beta"]],
                        groups["vehicle"]["mat"][groups["vehicle"]["beta"]])

    # Test the difference of effects with a conservative z on the effect
    # difference, using the pooled SE of both groups' means.
    n_a1, n_a0 = len(groups["iso"]["alpha"]), len(groups["vehicle"]["alpha"])
    n_b1, n_b0 = len(groups["iso"]["beta"]), len(groups["vehicle"]["beta"])
    A1 = groups["iso"]["mat"][groups["iso"]["alpha"]]
    A0 = groups["vehicle"]["mat"][groups["vehicle"]["alpha"]]
    B1 = groups["iso"]["mat"][groups["iso"]["beta"]]
    B0 = groups["vehicle"]["mat"][groups["vehicle"]["beta"]]
    se_diff = np.sqrt(A1.var(axis=0, ddof=1) / n_a1 + A0.var(axis=0, ddof=1) / n_a0
                      + B1.var(axis=0, ddof=1) / n_b1 + B0.var(axis=0, ddof=1) / n_b0)
    se_diff[se_diff == 0] = np.inf
    z_diff = (A1.mean(axis=0) - A0.mean(axis=0)
              - (B1.mean(axis=0) - B0.mean(axis=0))) / se_diff
    from math import erfc, sqrt
    pdiff = np.array([erfc(abs(x) / sqrt(2.0)) for x in z_diff])
    qdiff = bh(pdiff)

    n_sig = int((qdiff < FDR).sum())
    print(f"  genes tested: {n_genes}")
    print(f"  genes with a significantly DIFFERENT drug effect in alpha vs "
          f"beta (BH q<{FDR}): {n_sig}")
    print(f"  expected by chance at q<{FDR}: ~{int(n_genes*FDR)}")
    print()
    if n_sig > 0:
        idx = np.where(qdiff < FDR)[0]
        idx = idx[np.argsort(qdiff[idx])]
        print(f"  top subtyping-responsive genes:")
        for j in idx[:20]:
            ea = eff["alpha"][j]
            eb = eff["beta"][j]
            print(f"    {syms[j]:<12} alpha {ea:+8.2f}  beta {eb:+8.2f}  "
                  f"diff {ea-eb:+8.2f}  q={qdiff[j]:.2e}")
    else:
        print("  -> No gene shows a subtype-specific drug response at BH q<0.05.")
        print("     Alpha and beta respond to isoproterenol INDISTINGUISHABLY")
        print("     at the transcriptional level in this dataset.")

    # Sanity check: are the groups actually responding at all?
    print(f"\nSANITY CHECK -- is there a drug effect at all?")
    for g, (I_, V_) in (("alpha", (A1, A0)), ("beta", (B1, B0))):
        q = bh(per_gene_p(I_, V_)[0])
        print(f"  {g}: {int((q < FDR).sum())} genes responsive to drug "
              f"(BH q<{FDR}), of {n_genes}")
    print("  (Mays report Aanat, Crem, Drd4, Pde10a among the responders)")

    print(f"\n  Mays-named responders, raw mean CP10K:")
    print(f"    {'gene':<9}{'a_veh':>9}{'a_iso':>9}{'b_veh':>9}{'b_iso':>9}")
    for gname in ("Aanat", "Crem", "Drd4", "Pde10a", "Asmt", "Tph1"):
        if gname not in gi:
            continue
        j = gi[gname]
        print(f"    {gname:<9}{A0[:,j].mean():>9.2f}{A1[:,j].mean():>9.2f}"
              f"{B0[:,j].mean():>9.2f}{B1[:,j].mean():>9.2f}")

    out = ROOT / "results" / "iso_response.json"
    out.write_text(json.dumps({
        "design": "unpaired, different animals",
        "shared_barcodes": len(shared),
        "alpha_cells": {c: int(len(groups[c]["alpha"])) for c in groups},
        "beta_cells": {c: int(len(groups[c]["beta"])) for c in groups},
        "n_genes_tested": n_genes,
        "n_subtype_specific_q05": n_sig,
        "expected_by_chance": int(n_genes * FDR),
    }, indent=2))
    print(f"\n  wrote {out}")


def load_fn(fname):
    return load(fname)


if __name__ == "__main__":
    main()
