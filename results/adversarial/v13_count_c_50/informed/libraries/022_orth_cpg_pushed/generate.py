"""
Experiment 022 — push CpG-rich slice further

021 added 6 k PLS+pELS to 015's mix, ties eval_01 (0.5746) + eval_08
+0.044. This explores whether more CpG content can keep eval_01
while further lifting eval_08.

Design:
  20 k orth-DHS (numsamples ≤ 5)
  20 k cCRE class-balanced
  10 k PLS+pELS supplement (5 k each)

Test point: 20 % CpG slice — double 021's 12 %.

If eval_01 still holds: more CpG is free; push toward 015→021 family.
If eval_01 falls: the slice that maximizes eval_01 is small (~6 k);
021's design was the sweet spot.
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
N_DHS = 20_000
N_CCRE = 20_000
N_CPG = 10_000
PER_CLS = N_CCRE // 8  # 2500
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
    pool = np.where(mask & (nsamp <= MAX_NSAMP))[0]
    print(f"orth-DHS pool: {len(pool):,}")
    dhs_idx = rng.choice(pool, size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    cls_array = ccre_meta["cls"].to_numpy()
    parts = []
    cb_used = set()
    for c in sorted(np.unique(cls_array)):
        cpool = np.where(cls_array == c)[0]
        sel = rng.choice(cpool, size=PER_CLS, replace=False) if len(cpool) >= PER_CLS else cpool.copy()
        parts.append(sel)
        cb_used.update(sel.tolist())
    ccre_idx = np.concatenate(parts)
    if len(ccre_idx) < N_CCRE:
        extra_pool = np.array([i for i in np.where(cls_array == "dELS")[0] if i not in cb_used])
        ccre_idx = np.concatenate(
            [ccre_idx, rng.choice(extra_pool, size=N_CCRE - len(ccre_idx), replace=False)]
        )
        cb_used.update(ccre_idx.tolist())
    ccre_rows = np.asarray(ccre_seqs[ccre_idx])

    cpg_pool = np.array([i for i in np.where(np.isin(cls_array, ["PLS", "pELS"]))[0] if i not in cb_used])
    print(f"CpG-rich pool (PLS+pELS, not yet used): {len(cpg_pool):,}")
    cpg_idx = rng.choice(cpg_pool, size=N_CPG, replace=False)
    cpg_rows = np.asarray(ccre_seqs[cpg_idx])

    all_rows = np.concatenate([dhs_rows, ccre_rows, cpg_rows], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} ({N_DHS} orth + {N_CCRE} cCRE + {N_CPG} PLS/pELS)")


if __name__ == "__main__":
    main()
