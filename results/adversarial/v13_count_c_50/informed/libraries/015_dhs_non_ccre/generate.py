"""
Experiment 015 — DHS non-overlapping cCRE + cCRE class-balanced

The 008/011/012/014 plateau at 0.567–0.569 suggests effective unique
diversity in those mixes may be lower than 50 k because DHS and cCRE
overlap heavily (cCREs are derived from DHS + chromatin marks).

This experiment forces the DHS half to be *orthogonal* to cCRE by
removing any DHS whose summit is within 200 bp of any cCRE midpoint.
The DHS-not-in-cCRE pool represents open chromatin without the
classical regulatory marks SCREEN uses — includes potential
insulators, undiscovered regulatory regions, or low-confidence sites.

Design:
  25 k DHS-only (filtered to non-cCRE-overlapping) uniform
  25 k cCRE class-balanced
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
N_DHS = 25_000
N_CCRE = 25_000
PER_CLS = N_CCRE // 8
WINDOW = 200  # consider any DHS within this distance of a cCRE midpoint as "overlapping"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def non_ccre_dhs_mask(dhs_meta: pd.DataFrame, ccre_meta: pd.DataFrame, window: int) -> np.ndarray:
    """Boolean mask of DHS rows where summit is >= `window` from any cCRE midpoint on same chrom."""
    mask = np.ones(len(dhs_meta), dtype=bool)
    # Build dict: chrom -> sorted np.array of cCRE midpoints
    ccre_by_chrom = {c: np.sort(g["mid"].to_numpy()) for c, g in ccre_meta.groupby("chrom")}
    dhs_chrom = dhs_meta["chrom"].to_numpy()
    dhs_summit = dhs_meta["summit"].to_numpy()
    for i, (chrom, summit) in enumerate(zip(dhs_chrom, dhs_summit)):
        mids = ccre_by_chrom.get(chrom)
        if mids is None:
            continue  # no cCRE on this chrom; keep
        # find insertion point
        j = np.searchsorted(mids, summit)
        # check left and right neighbors
        if j < len(mids) and mids[j] - summit < window:
            mask[i] = False
            continue
        if j > 0 and summit - mids[j-1] < window:
            mask[i] = False
        if i % 500_000 == 0:
            print(f"  progress {i:,}/{len(dhs_meta):,}")
    return mask


def main():
    rng = np.random.default_rng(SEED)
    dhs_meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")

    print(f"DHS={len(dhs_meta):,}, cCRE={len(ccre_meta):,}")
    print(f"Computing DHS non-overlap mask (window={WINDOW} bp)...")
    mask = non_ccre_dhs_mask(dhs_meta, ccre_meta, WINDOW)
    pool = np.where(mask)[0]
    print(f"DHS-non-cCRE pool: {len(pool):,} ({100*mask.mean():.1f}% of DHS)")

    dhs_idx = rng.choice(pool, size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    parts = []
    for c in sorted(ccre_meta["cls"].unique()):
        cpool = np.where(ccre_meta["cls"].to_numpy() == c)[0]
        sel = rng.choice(cpool, size=PER_CLS, replace=False) if len(cpool) >= PER_CLS else cpool.copy()
        parts.append(sel)
    ccre_idx = np.concatenate(parts)
    if len(ccre_idx) < N_CCRE:
        used = set(ccre_idx.tolist())
        extra_pool = np.array([i for i in np.where(ccre_meta["cls"].to_numpy() == "dELS")[0] if i not in used])
        ccre_idx = np.concatenate(
            [ccre_idx, rng.choice(extra_pool, size=N_CCRE - len(ccre_idx), replace=False)]
        )
    ccre_rows = np.asarray(ccre_seqs[ccre_idx])

    all_rows = np.concatenate([dhs_rows, ccre_rows], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,}")


if __name__ == "__main__":
    main()
