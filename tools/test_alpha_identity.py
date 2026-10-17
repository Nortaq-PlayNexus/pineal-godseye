"""Are the Asmt-high cells actually pinealocytes, and does the result replicate?

TWO QUESTIONS, BOTH ANSWERABLE FROM THE RAW MATRICES.

Q1. IDENTITY. The alpha population was identified by ranking on Asmt. If those
cells are a contaminating non-pinealocyte type that merely expresses Asmt, the
whole alpha/beta finding is meaningless. Test: are Asmt-high cells enriched for
INDEPENDENT pinealocyte markers that were NOT used to define them (Tph1, Crx,
Gngt1, Gngt2, Nr1d1, Rax, Otx2)? Enrichment on markers that played no part in
the selection is not circular.

Q2. REPLICATION. The first run used Day1 only. Does the complexity signature
reproduce in Night1, an independent animal and circadian state?

Also reports whether the "transcriptionally quiet" description of alpha cells
holds, given they carry MORE total UMI than any other group. Quiet in
ribosomal terms is not the same as quiet overall.
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
GEO = ROOT / "data" / "geo"
MIN_UMI, MIN_GENES = 500, 250

# Rnor6.0 symbols. "Mt-nd1"-style for mitochondria (NOT "mt-").
MITO = "Mt-"
RIBO = ("Rpl", "Rps")

# Pinealocyte markers that played NO part in defining the alpha slice.
#
# CORRECTION: an earlier version used Tph1, Gngt1 and Gngt2 as identity controls
# and scored 2/7 and 3/7 enriched. That was a badly designed test, not evidence
# against identity. Mays 2018 states alpha and beta share an "overlapping set of
# marker genes" and specifically lists Tph1, Gngt1, Gngt2, Crx, Cngb1, Cnga1
# and Pde6c as expressed in BOTH subtypes. A marker shared by both subtypes
# cannot discriminate between them, so near-1.0 ratios are the EXPECTED result
# and carry no information about contamination.
#
# The right control is the REST of the dataset: if Asmt-high cells are
# pinealocytes, they must be far more pinealocyte-like than non-pinealocytes.
# That contrast is not circular -- the slice was chosen on Asmt alone.
INDEPENDENT_MARKERS = ["Tph1", "Crx", "Gngt1", "Gngt2", "Nr1d1", "Rax",
                       "Cngb1", "Pde6c", "Sag"]


def load(sample):
    d = GEO / sample
    syms = []
    with gzip.open(d / f"{sample}_genes.tsv.gz", "rt") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                syms.append(p[1])
    with gzip.open(d / f"{sample}_matrix.mtx.gz", "rt") as f:
        for ln in f:
            q = ln.split()
            if len(q) == 3 and all(x.isdigit() for x in q):
                ng, nc, nz = (int(x) for x in q)
                break
        rr = np.empty(nz, np.int32)
        cc = np.empty(nz, np.int32)
        vv = np.empty(nz, np.int32)
        k = 0
        for ln in f:
            if not ln or ln.startswith("%"):
                continue
            a, b, v = ln.split()
            rr[k] = int(a) - 1
            cc[k] = int(b) - 1
            vv[k] = int(v)
            k += 1
    umi = np.zeros(nc, np.int64)
    det = np.zeros(nc, np.int64)
    np.add.at(umi, cc, vv)
    np.add.at(det, cc, 1)
    keep = (umi >= MIN_UMI) & (det >= MIN_GENES)
    remap = np.full(nc, -1, np.int64)
    remap[keep] = np.arange(int(keep.sum()))
    ek = keep[cc]
    r, c, v = rr[ek], remap[cc[ek]], vv[ek]
    return np.array(syms[:ng]), r, c, v, umi[keep], int(keep.sum())


def analyse(sample, frac=0.05):
    print("=" * 74)
    print(f"SAMPLE {sample}")
    print("=" * 74)
    sym, r, c, v, umi, n = load(sample)
    print(f"  cells after filter: {n:,}")

    idx = {s: i for i, s in enumerate(sym)}
    asmt_i = idx["Asmt"]
    mito_sel = np.array([s.startswith(MITO) for s in sym])
    ribo_sel = np.array([s.startswith(RIBO) for s in sym])
    print(f"  Mt- genes: {int(mito_sel.sum())}  Rpl/Rps genes: {int(ribo_sel.sum())}")

    def col(mask):
        out = np.zeros(n, np.int64)
        np.add.at(out, c, np.where(mask[r], v, 0))
        return out

    asmt = col(sym[r] == asmt_i) if False else None  # placeholder, replaced below
    asmt = np.zeros(n, np.int64)
    np.add.at(asmt, c, np.where(r == asmt_i, v, 0))
    mito = col(mito_sel)
    ribo = col(ribo_sel)

    pin = asmt > 0
    pidx = np.where(pin)[0]
    order = pidx[np.argsort(asmt[pidx])[::-1]]
    k = max(1, int(round(frac * len(order))))
    hi, rest = order[:k], order[k:]

    res = {"sample": sample, "n_cells": n, "n_pineal": int(pin.sum()),
           "n_alpha_slice": k}

    print(f"\n  pinealocytes (Asmt+): {int(pin.sum()):,}   alpha slice (top "
          f"{int(frac*100)}%): {k}")

    # --- Q1 identity ---------------------------------------------------------
    print(f"\nQ1. IDENTITY -- enrichment for markers NOT used to define the slice")
    print(f"    {'marker':<10}{'Asmt-high':>12}{'all pineal':>13}{'ratio':>9}")
    print("    " + "-" * 44)
    ident = {}
    print(f"    {'marker':<10}{'Asmt-high':>12}{'pinealocyte':>14}"
          f"{'non-pineal':>13}{'vs non-pineal':>15}")
    print("    " + "-" * 64)
    nonpin = np.where(~pin)[0]
    for m in INDEPENDENT_MARKERS:
        if m not in idx:
            print(f"    {m:<10}{'(absent from annotation)':>40}")
            continue
        mc = np.zeros(n, np.int64)
        np.add.at(mc, c, np.where(r == idx[m], v, 0))
        # counts per 10k UMI, so library depth cannot drive the comparison
        a = 1e4 * mc[hi].mean() / umi[hi].mean()
        b = 1e4 * mc[pidx].mean() / umi[pidx].mean()
        o = 1e4 * mc[nonpin].mean() / umi[nonpin].mean() if len(nonpin) else float("nan")
        vs_non = a / o if o and o > 0 else float("nan")
        ident[m] = float(vs_non)
        flag = "" if (np.isnan(vs_non) or vs_non >= 1.5) else "  <-- not elevated"
        print(f"    {m:<10}{a:>12.3f}{b:>14.3f}{o:>13.3f}{vs_non:>15.2f}{flag}")
    res["identity_vs_nonpinealocyte"] = ident
    strong = [m for m, v in ident.items() if v >= 1.5]
    print(f"\n    {len(strong)}/{len(ident)} markers elevated >=1.5x versus "
          f"non-pinealocytes")
    if len(strong) >= 4:
        print("    -> Asmt-high cells are far more pinealocyte-like than the")
        print("       rest of the tissue. Identity is consistent with pinealocyte.")
    else:
        print("    -> IDENTITY NOT ESTABLISHED. Treat the alpha result as")
        print("       provisional pending a better marker panel.")

    # --- Q2 complexity signature --------------------------------------------
    print(f"\nQ2. COMPLEXITY SIGNATURE (replication of the Day1 result)")
    print(f"    {'group':<14}{'UMI':>10}{'ribo%':>9}{'mito%':>9}")
    print("    " + "-" * 42)
    for name, idxs in [("Asmt-high", hi), ("all pinealocyte", pidx)]:
        print(f"    {name:<14}{umi[idxs].mean():>10.0f}"
              f"{100*ribo[idxs].mean()/umi[idxs].mean():>8.2f}%"
              f"{100*mito[idxs].mean()/umi[idxs].mean():>8.2f}%")
    rf_hi = 100 * ribo[hi].mean() / umi[hi].mean()
    rf_all = 100 * ribo[pidx].mean() / umi[pidx].mean()
    umi_hi = float(umi[hi].mean())
    umi_all = float(umi[pidx].mean())
    res.update(rf_hi=float(rf_hi), rf_all=float(rf_all),
               umi_hi=umi_hi, umi_all=umi_all,
               depth_ratio=umi_hi / umi_all)
    print(f"    depth ratio (Asmt-high / all pinealocyte) = {umi_hi/umi_all:.2f}")
    print(f"    ribosomal fraction ratio                   = {rf_hi/rf_all:.2f}")
    if umi_hi > umi_all and rf_hi < rf_all:
        print("    -> REPLICATED: deeper sequencing AND lower ribosomal fraction.")
        print("       Opposite of the low-complexity artefact signature.")
    else:
        print("    -> DID NOT REPLICATE on this sample.")

    # --- the 'transcriptionally quiet' question ------------------------------
    nonribo_hi = 100 - 100 * ribo[hi].mean() / umi[hi].mean()
    nonribo_all = 100 - 100 * ribo[pidx].mean() / umi[pidx].mean()
    print(f"\n    'Transcriptionally quiet' check:")
    print(f"      non-ribosomal transcript share: Asmt-high {nonribo_hi:.1f}% vs "
          f"all pinealocyte {nonribo_all:.1f}%")
    print(f"      absolute ribo counts/cell    : Asmt-high "
          f"{ribo[hi].mean():.0f} vs all {ribo[pidx].mean():.0f}")
    if ribo[hi].mean() > ribo[pidx].mean():
        print("      -> NOT quiet in absolute terms: Asmt-high cells carry MORE")
        print("         ribosomal counts, not fewer. 'Quiet' is true only in")
        print("         RELATIVE terms, and they are the deepest-sequenced cells.")
    else:
        print("      -> genuinely fewer ribosomal counts; 'quiet' holds absolutely.")
    res["ribo_abs_hi"] = float(ribo[hi].mean())
    res["ribo_abs_all"] = float(ribo[pidx].mean())
    print()
    return res


def main():
    all_res = []
    for s in ("GSM3188343", "GSM3188347"):
        all_res.append(analyse(s))
        print()
    print("=" * 74)
    print("CROSS-SAMPLE SUMMARY")
    print("=" * 74)
    print(f"{'sample':<14}{'cells':>8}{'depth r':>10}{'ribo r':>9}"
          f"{'markers vs non-pineal':>22}")
    for r in all_res:
        vals = r["identity_vs_nonpinealocyte"]
        enr = sum(1 for v in vals.values() if v >= 1.5)
        print(f"{r['sample']:<14}{r['n_cells']:>8}{r['depth_ratio']:>10.2f}"
              f"{r['rf_hi']/r['rf_all']:>9.2f}{f'{enr}/{len(vals)}':>22}")
    out = ROOT / "results" / "identity_and_replication.json"
    out.write_text(json.dumps(all_res, indent=2))
    print(f"\n  wrote {out}")


if __name__ == "__main__":
    main()
