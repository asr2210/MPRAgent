"""
Experiment 009 — ccre_active_only

50 000 sequences drawn uniformly at random from cCRE classes that have
explicit "active regulatory" chromatin signatures: PLS (promoter-like),
pELS (proximal enhancer-like), CA-H3K4me3 (open + H3K4me3). Pool size
~376 K.

Excluded: dELS (distal enhancer, weaker signal), CA / CA-CTCF / CA-TF
(only DNase / one accessory mark), TF-only.

Hypothesis: per-element information density is higher in actively
regulatory classes than in dELS-heavy uniform cCRE. Test whether the
"shrink the pool to highly active elements" intuition (which failed
for DHS in exp 004 because it correlated with cell-type-invariance)
works when the pool is *functionally* defined rather than signal-
thresholded.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SEED = 0
N = 50_000
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"
ACTIVE = {"PLS", "pELS", "CA-H3K4me3"}


def main():
    rng = np.random.default_rng(SEED)
    meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    keep = meta["cls"].isin(ACTIVE).to_numpy()
    pool = np.where(keep)[0]
    print(f"Active-cCRE pool: {len(pool):,} (PLS+pELS+CA-H3K4me3)")
    sel = rng.choice(pool, size=N, replace=False)
    chosen = seqs[sel]
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")


if __name__ == "__main__":
    main()
