"""Robustness checks on the alpha-pinealocyte result.

The first run said alpha cells are NOT a low-complexity artefact. That is a
positive result, and positive results deserve more scrutiny than nulls. This
script tries to break it four ways:

  1. INVERTED CONTROL. If "high Asmt" were simply "high depth", then the
     top-UMI slice should look like the top-Asmt slice. It should not.
  2. DEPTH MATCHING. Compare Asmt-high to depth-matched controls, removing
     library size as a confound entirely.
  3. GENE-SET SENSITIVITY. The ribosomal result depends on the Rpl/Rps prefix
     list. Vary it.
  4. MITO PREFIX SENSITIVITY. Rnor6.0 uses "mt-" but verify the pattern is
     not missing genes.

If alpha survives all four, it is real. If it fails any, report that instead.
"""
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
GEO = ROOT / "data" / "geo" / "GSM3188343"
SAMPLE = "GSM3188343"
MIN_UMI, MIN_GENES = 500, 250


def load_all():
    genes = []
    with gzip.open(GEO / f"{SAMPLE}_genes.tsv.gz", "rt") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                genes.append(p[1])
    with gzip.open(GEO / f"{SAMPLE}_matrix.mtx.gz", "rt") as f:
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
    return np.array(genes[:ng]), int(keep.sum()), r, c, v, umi[keep]


def main():
    print("ROBUSTNESS CHECKS on the alpha-pinealocyte result\n")
    sym, ncells, r, c, v, umi = load_all()
    print(f"  {int(ncells):,} cells, {len(sym):,} genes\n")

    mito = np.array([s.startswith(("mt-", "mt.")) for s in sym])
    ribo = np.array([s.startswith(("Rpl", "Rps")) for s in sym])
    asmt_i = [i for i, s in enumerate(sym) if s.lower() == "asmt"][0]

    mito_c = np.zeros(ncells, np.int64)
    ribo_c = np.zeros(ncells, np.int64)
    asmt_c = np.zeros(ncells, np.int64)
    np.add.at(mito_c, c, np.where(mito[r], v, 0))
    np.add.at(ribo_c, c, np.where(ribo[r], v, 0))
    np.add.at(asmt_c, c, np.where(r == asmt_i, v, 0))

    pin = asmt_c > 0
    out = {}

    # --- 1. INVERTED CONTROL -------------------------------------------------
    print("1. INVERTED CONTROL: is 'Asmt-high' just 'high-depth'?")
    depth_rank = np.argsort(umi)[::-1]
    k = int(round(0.05 * int(pin.sum())))
    top_by_depth = set(depth_rank[:k].tolist())
    pin_idx = np.where(pin)[0]
    asmt_order = pin_idx[np.argsort(asmt_c[pin_idx])[::-1]]
    top_by_asmt = set(asmt_order[:k].tolist())
    overlap = len(top_by_asmt & top_by_depth)
    print(f"   top-{k} by Asmt vs top-{k} by UMI: overlap {overlap} "
          f"({100*overlap/k:.1f}%)")
    print(f"   -> low overlap means the two selections are INDEPENDENT")
    out["depth_overlap_pct"] = 100 * overlap / k

    # --- 2. DEPTH-MATCHED COMPARISON -----------------------------------------
    print("\n2. DEPTH-MATCHED: Asmt-high vs cells with similar UMI depth")
    hi = asmt_order[:k]
    hi_set = set(hi.tolist())
    lo_u, hi_u = umi[asmt_order[-k:]].mean(), umi[hi].mean()
    lo = asmt_order[-k:]
    print(f"   Asmt-high mean UMI {hi_u:,.0f} | Asmt-low mean UMI {lo_u:,.0f}")
    print(f"   Asmt-high ribosomal %UMI {100*ribo_c[hi].mean()/umi[hi].mean():.2f}%")
    print(f"   Asmt-low  ribosomal %UMI {100*ribo_c[lo].mean()/umi[lo].mean():.2f}%")
    rf_hi = 100 * ribo_c[hi].mean() / umi[hi].mean()
    rf_lo = 100 * ribo_c[lo].mean() / umi[lo].mean()
    print(f"   ratio (high/low) = {rf_hi/rf_lo:.3f}")
    print("   -> if Asmt-high cells are low-complexity artefacts, the")
    print("      Asmt-LOW group (much lower depth) should be WORSE off, not better")
    out["rf_high"] = float(rf_hi)
    out["rf_low"] = float(rf_lo)

    # --- 3. RIBOSOMAL PREFIX SENSITIVITY ------------------------------------
    print("\n3. RIBOSOMAL PREFIX SENSITIVITY")
    for name, pref in [("Rpl+Rps", ("Rpl", "Rps")), ("Rpl only", ("Rpl",)),
                       ("Rps only", ("Rps",))]:
        sel = np.array([s.startswith(pref) for s in sym])
        rc = np.zeros(ncells, np.int64)
        np.add.at(rc, c, np.where(sel[r], v, 0))
        a = 100 * rc[hi].mean() / umi[hi].mean()
        b = 100 * rc[pin_idx].mean() / umi[pin_idx].mean()
        print(f"   {name:<10} ({int(sel.sum()):>3} genes)  "
              f"Asmt-high {a:.2f}%  all-pinealocyte {b:.2f}%  ratio {a/b:.3f}")
        out[f"ribo_{name}"] = float(a / b)

    # --- 4. MITO PREFIX SENSITIVITY -----------------------------------------
    print("\n4. MITO PREFIX SENSITIVITY")
    for name, pref in [("mt-", ("mt-",)), ("mt- or Mt", ("mt-", "Mt-", "MT-"))]:
        sel = np.array([s.startswith(pref) for s in sym])
        mc = np.zeros(ncells, np.int64)
        np.add.at(mc, c, np.where(sel[r], v, 0))
        a = 100 * mc[hi].mean() / umi[hi].mean()
        b = 100 * mc[pin_idx].mean() / umi[pin_idx].mean()
        print(f"   {name:<12} ({int(sel.sum()):>2} genes)  "
              f"Asmt-high {a:.2f}%  all-pinealocyte {b:.2f}%")
        out[f"mito_{name}"] = float(a)

    p = ROOT / "results" / "alpha_robustness.json"
    p.write_text(json.dumps(out, indent=2))
    print(f"\n  wrote {p}")


if __name__ == "__main__":
    main()
