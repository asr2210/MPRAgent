"""
Experiment 019 — pure orthogonal DHS, 50k

Isolate the orth-DHS effect: 50 k drawn uniformly from the 950 k
DHS-non-cCRE pool. No cCRE half.

Baselines:
- 002 (full DHS uniform 50k): 0.5627
- 015 (25k orth + 25k cCRE): 0.5736
- 016 (35k orth + 15k cCRE): 0.5689
- 006 (50k cCRE class-bal): 0.5637

If 019 > 002: orth-DHS pool is intrinsically better than full DHS
(not just because it deduplicates cCRE — it's better content).
If 019 < 002: orth-DHS only helps in combination with cCRE.
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
WINDOW = 200
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def non_ccre_dhs_mask(dhs_meta, ccre_meta, window):
    mask = np.ones(len(dhs_meta), dtype=bool)
    ccre_by_chrom = {c: np.sort(g["mid"].to_numpy()) for c, g in ccre_meta.groupby("chrom")}
    dhs_chrom = dhs_meta["chrom"].to_numpy()
    dhs_summit = dhs_meta["summit"].to_numpy()
    for i, (chrom, summit) in enumerate(zip(dhs_chrom, dhs_summit)):
        mids = ccre_by_chrom.get(chrom)
        if mids is None: continue
        j = np.searchsorted(mids, summit)
        if j < len(mids) and mids[j] - summit < window:
            mask[i] = False; continue
        if j > 0 and summit - mids[j-1] < window:
            mask[i] = False
    return mask


def main():
    rng = np.random.default_rng(SEED)
    dhs_meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    mask = non_ccre_dhs_mask(dhs_meta, ccre_meta, WINDOW)
    pool = np.where(mask)[0]
    print(f"Orth DHS pool: {len(pool):,}")
    idx = rng.choice(pool, size=N, replace=False)
    rows = np.asarray(dhs_seqs[idx])
    rng.shuffle(rows, axis=0)
    with open(OUT, "wb") as f:
        for row in rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,}")


if __name__ == "__main__":
    main()
