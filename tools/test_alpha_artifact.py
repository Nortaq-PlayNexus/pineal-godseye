"""Is the alpha-pinealocyte cluster a real cell type or a droplet artefact?

THE QUESTION
Mays et al. 2018 define alpha-pinealocytes by ~3.4-fold higher Asmt, ~2.3-fold
higher OxPhos, and ~8.2-fold LOWER ribosomal expression than beta. Low ribosomal
plus high mitochondrial content is the canonical low-complexity / dying-cell
signature in droplet scRNA-seq, and a well-known driver of spurious subclusters.
Mays reported n=275 alpha cells and no per-cluster UMI distributions.

If alpha cells are an artefact of low library complexity, they should show:
  - systematically LOWER total UMI counts
  - systematically HIGHER mitochondrial fraction
  - a ribosomal deficit that is explained by overall complexity rather than
    being a specific biological feature

If alpha is real, the Asmt-high cells should have NORMAL complexity, and their
metabolic profile should differ without the whole cell being low-quality.

METHOD
The MTX file has no cell-type labels, so alpha cells are identified the way a
sceptic would: rank pinealocytes by Asmt expression and take the top slice. If
the alpha population is an artefact, those cells ARE the artefact. This makes
the test conservative -- it does not depend on trusting Mays' cluster labels.

RATIO-BASED, NOT CELL-RANK-BASED, COMPARISON
Mays' 8.2x ribosomal figure is a ratio of MEANS of normalised expression.
A low-complexity artefact inflates that ratio on its own: ribosomal genes drop
out first because they are numerous and low-expression, so as a cell's library
complexity falls, ribosomal detection falls with it. So a high ribosomal ratio
is NOT evidence of biology on its own. The discriminating quantity is the
ABSOLUTE ribosomal count per cell, compared against the UMI depth of the same
cell. A real alpha-pinealocyte runs its own ribosomes and should retain
ribosomal counts proportional to its library. An artefact should show ribosomal
counts falling WITH library depth, not holding a constant fraction.
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
GEO = ROOT / "data" / "geo" / "GSM3188343"
SAMPLE = "GSM3188343"

MITO_PREFIXES = ("mt-", "mt.", "Mt-")
RIBO_PREFIXES = ("Rpl", "Rps")

MIN_UMI = 500
MIN_GENES = 250


def load():
    genes = []
    with gzip.open(GEO / f"{SAMPLE}_genes.tsv.gz", "rt") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                genes.append(p[1])

    with gzip.open(GEO / f"{SAMPLE}_matrix.mtx.gz", "rt") as f:
        for ln in f:
            parts = ln.split()
            if len(parts) == 3 and all(p.isdigit() for p in parts):
                n_genes, n_cells, nnz = (int(p) for p in parts)
                break
        else:
            raise ValueError("no dimension line")
        rr = np.empty(nnz, np.int32)
        cc = np.empty(nnz, np.int32)
        vv = np.empty(nnz, np.int32)
        k = 0
        for ln in f:
            if not ln or ln.startswith("%"):
                continue
            a, b, v = ln.split()
            rr[k] = int(a) - 1
            cc[k] = int(b) - 1
            vv[k] = int(v)
            k += 1
    assert k == nnz

    sym = np.array(genes[:n_genes], dtype=object)
    mito = np.array([any(s.startswith(p) for p in MITO_PREFIXES) for s in sym])
    ribo = np.array([s.startswith(RIBO_PREFIXES) for s in sym])
    print(f"  mitochondrial genes matched : {int(mito.sum())}")
    print(f"  ribosomal genes matched     : {int(ribo.sum())}")

    umi = np.zeros(n_cells, np.int64)
    detected = np.zeros(n_cells, np.int64)
    np.add.at(umi, cc, vv)
    np.add.at(detected, cc, 1)
    keep = (umi >= MIN_UMI) & (detected >= MIN_GENES)

    # `keep` is indexed by CELL; the entry arrays are indexed by NONZERO.
    # Expand it to entry space before masking, or the shapes will not match.
    entry_keep = keep[cc]
    cc_k = cc[entry_keep]
    vv_k = vv[entry_keep]
    rr_k = rr[entry_keep]
    umi_k = umi[keep]

    # Remap original cell IDs to positions within the kept subset, so the
    # per-cell accumulators below are indexed consistently.
    remap = np.full(n_cells, -1, dtype=np.int64)
    remap[keep] = np.arange(int(keep.sum()))
    cc_k = remap[cc_k]

    mito_counts = np.zeros(len(umi_k), np.int64)
    ribo_counts = np.zeros(len(umi_k), np.int64)
    np.add.at(mito_counts, cc_k, np.where(mito[rr_k], vv_k, 0))
    np.add.at(ribo_counts, cc_k, np.where(ribo[rr_k], vv_k, 0))

    asmt_idx = [i for i, s in enumerate(sym[:n_genes]) if s.lower() == "asmt"]
    asmt_counts = np.zeros(len(umi_k), np.int64)
    if asmt_idx:
        np.add.at(asmt_counts, cc_k, np.where(rr_k == asmt_idx[0], vv_k, 0))

    return dict(umi=umi_k, detected=detected[keep], mito=mito_counts,
                ribo=ribo_counts, asmt=asmt_counts, n_kept=len(umi_k))


def pct(a, q):
    return float(np.percentile(a, q))


def main():
    print("GSE115723 / GSM3188343 (Day1) -- alpha-pinealocyte artefact test\n")
    d = load()
    n = d["n_kept"]
    print(f"  cells passing filter (UMI>={MIN_UMI}, genes>={MIN_GENES}): {n:,}")
    print()

    pineal = d["asmt"] > 0
    print(f"  cells with detectable Asmt: {int(pineal.sum()):,} "
          f"({100*pineal.sum()/n:.1f}% of kept cells)")
    print()

    # Rank Asmt-positive pinealocytes by Asmt counts and take the top slice.
    # Mays reports alpha as 5% of pinealocytes. Slice at several depths so the
    # conclusion does not depend on one cutoff.
    pmask = pineal
    order = np.argsort(d["asmt"][pmask])[::-1]
    asmt_sorted = d["asmt"][pmask][order]
    umi_s = d["umi"][pmask][order]
    ribo_s = d["ribo"][pmask][order]
    mito_s = d["mito"][pmask][order]
    n_pin = len(asmt_sorted)

    print("ASMT-HIGH vs REMAINDER (within Asmt-positive pinealocytes)")
    print("=" * 74)
    print(f"{'slice':<10}{'n':>7}{'UMI/cell':>11}{'ribo/cell':>12}"
          f"{'ribo %UMI':>12}{'mito %UMI':>12}")
    print("-" * 74)

    results = {}
    for frac in (0.02, 0.05, 0.10, 0.20):
        k = max(1, int(round(frac * n_pin)))
        top_u, top_r, top_m = umi_s[:k].mean(), ribo_s[:k].mean(), mito_s[:k].mean()
        rest_u = umi_s[k:].mean()
        rest_r = ribo_s[k:].mean()
        rest_m = mito_s[k:].mean()
        rf_top = 100 * top_r / top_u
        rf_rest = 100 * rest_r / rest_u
        mf_top = 100 * top_m / top_u
        mf_rest = 100 * rest_m / rest_u
        results[frac] = dict(
            k=k, umi_top=float(top_u), umi_rest=float(rest_u),
            ribo_top=float(top_r), ribo_rest=float(rest_r),
            ribo_pct_top=float(rf_top), ribo_pct_rest=float(rf_rest),
            mito_pct_top=float(mf_top), mito_pct_rest=float(mf_rest))
        print(f"top {int(frac*100):<6}%{k:>7}{top_u:>11.0f}{top_r:>12.1f}"
              f"{rf_top:>11.2f}%{mf_top:>11.2f}%")
        print(f"{'rest':<10}{n_pin-k:>7}{rest_u:>11.0f}{rest_r:>12.1f}"
              f"{rf_rest:>11.2f}%{mf_rest:>11.2f}%")
        print()

    r5 = results[0.05]
    print("INTERPRETATION (top 5% slice, matching Mays' alpha fraction)")
    print("-" * 74)
    umi_ratio = r5["umi_top"] / r5["umi_rest"]
    print(f"  UMI/cell          Asmt-high {r5['umi_top']:.0f} vs rest {r5['umi_rest']:.0f}"
          f"   ratio {umi_ratio:.2f}")
    print(f"  ribosomal %UMI    Asmt-high {r5['ribo_pct_top']:.2f}% vs "
          f"rest {r5['ribo_pct_rest']:.2f}%   ratio "
          f"{r5['ribo_pct_top']/r5['ribo_pct_rest']:.2f}")
    print(f"  mitochondrial %   Asmt-high {r5['mito_pct_top']:.2f}% vs "
          f"rest {r5['mito_pct_rest']:.2f}%")
    print()

    verdict = []
    if umi_ratio < 0.85:
        verdict.append("Asmt-high cells have LOWER library depth -> artefact signature")
    else:
        verdict.append("Asmt-high cells have comparable library depth -> NOT an "
                       "artefact of low complexity")
    if r5["mito_pct_top"] > r5["mito_pct_rest"] * 1.25:
        verdict.append("Asmt-high cells have substantially HIGHER mitochondrial "
                       "fraction -> artefact signature")
    else:
        verdict.append("mitochondrial fraction is comparable -> not a dying-cell "
                       "artefact")
    print("  " + "\n  ".join(verdict))

    # Also report the raw correlation between depth and ribosomal fraction,
    # which is the mechanistic reason a ratio alone is uninformative.
    depth = d["umi"][pmask].astype(float)
    ribo_frac = d["ribo"][pmask] / np.maximum(depth, 1)
    if depth.std() > 0 and ribo_frac.std() > 0:
        r = float(np.corrcoef(depth, ribo_frac)[0, 1])
        print(f"\n  corr(UMI depth, ribosomal fraction) across pinealocytes "
              f"= {r:+.3f}")
        print("  A strongly positive correlation means ribosomal fraction is")
        print("  driven by library complexity, so ANY group with lower depth")
        print("  will show an inflated ribosomal RATIO without that being biology.")

    out = ROOT / "results" / "alpha_artifact_test.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(results, indent=2))
    print(f"\n  wrote {out}")


if __name__ == "__main__":
    main()
