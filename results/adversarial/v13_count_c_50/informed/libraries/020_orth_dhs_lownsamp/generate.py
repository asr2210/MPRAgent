"""
Experiment 020 — orth-DHS cell-type-specific subset

Push the orthogonality direction. 015 confirmed orth-DHS adds value
because it brings non-K562 cell-type-specific accessibility. The
DHS metadata has numsamples (per-element biosample count) — peaks
seen in few samples ARE cell-type-specific by definition.

This experiment further filters orth-DHS to numsamples ≤ 5:
"non-cCRE AND cell-type-specific" — the maximally cell-type-
informative orthogonal subset.

Design:
  25 k from orth-DHS-numsamples≤5 pool + 25 k cCRE class-bal.

Predicts: HepG2/SK-N-SH further lift (cell-type-specific signal is
even stronger); K562 may dip slightly.

Note: exp 001's "w∝1/sqrt(numsamples)" with FULL DHS pool was a
slight loss (-0.003). This is different because we're sampling from
the orth pool first, then filtering. The orth pool is already
biased toward non-K562 — adding a hard numsamples cap should amplify
that bias rather than dilute it.
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
WINDOW = 200
MAX_NSAMP = 5
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
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")

    mask = non_ccre_dhs_mask(dhs_meta, ccre_meta, WINDOW)
    nsamp = dhs_meta["numsamples"].to_numpy()
    full_mask = mask & (nsamp <= MAX_NSAMP)
    pool = np.where(full_mask)[0]
    print(f"orth-DHS pool: {mask.sum():,}; orth & nsamp<={MAX_NSAMP}: {len(pool):,}")
    if len(pool) < N_DHS:
        print("  pool too small — relaxing MAX_NSAMP")
        # fall back: use orth-DHS uniform
        pool = np.where(mask)[0]
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
