"""
Experiment 018 — orth-DHS light + cCRE heavy (15/35)

015 (25/25) is the current best at eval_01=0.5736.
016 (35/15 orth-heavy) was 0.5689 (worse).
This tests the OTHER direction: 15k orth-DHS + 35k cCRE class-bal.

If 018 > 015: more cCRE side, less orth-DHS = better. Optimal ratio
is cCRE-heavy.
If 018 < 015: 25/25 is the optimum; orth-DHS contributes real
incremental value at that mix.
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
N_DHS = 15_000
N_CCRE = 35_000
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
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")

    mask = non_ccre_dhs_mask(dhs_meta, ccre_meta, WINDOW)
    pool = np.where(mask)[0]
    print(f"Orth DHS pool: {len(pool):,}")
    dhs_idx = rng.choice(pool, size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    classes = sorted(ccre_meta["cls"].unique())
    per_cls = N_CCRE // len(classes)
    parts = []
    for c in classes:
        cpool = np.where(ccre_meta["cls"].to_numpy() == c)[0]
        sel = rng.choice(cpool, size=per_cls, replace=False) if len(cpool) >= per_cls else cpool.copy()
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
    print(f"Wrote {N:,} ({N_DHS} orth + {N_CCRE} cCRE)")


if __name__ == "__main__":
    main()
