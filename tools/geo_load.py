"""Load GSE115723 GSM3188343 (Day1) and report basic structure.

This is the daytime sample Mays et al. 2018 used for the alpha/beta analysis
(5,667 cells reported in the paper). Nothing is inferred here; this only
establishes what is actually in the file before any QC statistic is computed.
"""
import gzip
from pathlib import Path

import numpy as np

GEO = Path(__file__).resolve().parent.parent / "data" / "geo" / "GSM3188343"
SAMPLE = "GSM3188343"


def read_genes():
    with gzip.open(GEO / f"{SAMPLE}_genes.tsv.gz", "rt") as f:
        rows = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    return rows


def read_barcodes():
    with gzip.open(GEO / f"{SAMPLE}_barcodes.tsv.gz", "rt") as f:
        return [ln.strip() for ln in f if ln.strip()]


def read_matrix():
    """Read a CellRanger MTX into a dense (genes x cells) int32 array."""
    path = GEO / f"{SAMPLE}_matrix.mtx.gz"
    with gzip.open(path, "rt") as f:
        # CellRanger writes: '%%MatrixMarket ...', '%', then 'G C nnz'.
        # The dimension line is the first that parses as three integers.
        for ln in f:
            parts = ln.split()
            if len(parts) == 3 and all(p.isdigit() for p in parts):
                n_genes, n_cells, nnz = (int(p) for p in parts)
                break
        else:
            raise ValueError("no dimension line found in MTX header")
        print(f"  MTX header: {n_genes} genes x {n_cells} cells, {nnz} nonzeros")
        print("  NOTE: this is the RAW UNFILTERED droplet matrix. The paper's")
        print("        5,667 cells are the subset that survives CellRanger")
        print("        filtering; the barcodes file lists all droplets.")
        rows = np.empty(nnz, dtype=np.int32)
        cols = np.empty(nnz, dtype=np.int32)
        vals = np.empty(nnz, dtype=np.int32)
        k = 0
        for ln in f:
            if not ln or ln.startswith("%"):
                continue
            a, b, v = ln.split()
            rows[k] = int(a) - 1
            cols[k] = int(b) - 1
            vals[k] = int(v)
            k += 1
            if k % 5_000_000 == 0:
                print(f"    parsed {k:,}/{nnz:,}")
    assert k == nnz, f"expected {nnz} entries, parsed {k}"
    # 737,280 droplets x 22,250 genes dense int32 would be 65 GB. Store the
    # coordinate arrays instead and materialise only the cells we need.
    return rows, cols, vals, n_genes, n_cells


def cell_stats(rows, cols, vals, n_genes, n_cells, min_umi=500, min_genes=250):
    """Compute per-cell UMI and detected-gene counts from coordinates."""
    umi = np.zeros(n_cells, dtype=np.int64)
    detected = np.zeros(n_cells, dtype=np.int64)
    np.add.at(umi, cols, vals)
    np.add.at(detected, cols, 1)
    keep = (umi >= min_umi) & (detected >= min_genes)
    return umi, detected, keep


def main():
    print("Reading GSE115723 / GSM3188343 (Day1, OriginalRef)\n")
    genes = read_genes()
    barcodes = read_barcodes()
    print(f"  genes.tsv    : {len(genes):,} rows, {len(genes[0])} cols")
    print(f"  barcodes.tsv : {len(barcodes):,} cells")
    print(f"  example gene row: {genes[0]}")
    print()
    rows, cols, vals, n_genes, n_cells = read_matrix()
    umi, detected, keep = cell_stats(rows, cols, vals, n_genes, n_cells)
    print(f"  total counts : {vals.sum():,}")
    print(f"  UMI/cell     : min={umi.min()}, median={int(np.median(umi))}, "
          f"max={umi.max()}")
    print()
    print(f"  droplets                 : {n_cells:,}")
    print(f"  passing UMI>=500         : {int(keep.sum()):,}")
    print(f"  passing UMI>=500,genes>=250: {int(keep.sum()):,}")
    print()
    print("The paper reports 5,667 cells after its own filtering. This matrix is")
    print("the RAW droplet set, so a cell count above 5,667 is expected and is")
    print("not a discrepancy -- Mays' exact filter settings are not published.")
    print()
    print("Next step: identify alpha-pinealocytes by Asmt expression and test")
    print("their complexity signature. See test_alpha_artifact.py.")


if __name__ == "__main__":
    main()
