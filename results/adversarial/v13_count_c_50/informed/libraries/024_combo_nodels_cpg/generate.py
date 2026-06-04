"""
Experiment 024 — stack all four winning priors

Combines:
- 015: orthogonal DHS-non-cCRE
- 020: numsamples ≤ 5 on orth pool
- 023: cCRE class-balanced WITHOUT dELS
- 021: 6 k PLS+pELS supplement

Design:
  22 k orth-DHS (numsamples ≤ 5)
  22 k cCRE balanced across 7 non-dELS classes (3142 each)
  6 k PLS+pELS supplement (3 k each)

But wait: 023's no-dELS cCRE already gives 3571 PLS, 3571 pELS. Adding
6 k more on top means ~6571 PLS and 6571 pELS total. That's a lot of
PLS+pELS. Let me check this doesn't crash diversity.

Alternative: drop the extra PLS+pELS slice — 023's design already
upweights PLS+pELS by removing dELS. Maybe 023 IS the right CpG
balance and 021's slice is redundant given 023.

Tests:
  024 vs 023: does adding CpG slice on top of no-dELS still help?
  024 vs 021: does no-dELS cCRE + CpG slice beat with-dELS cCRE + CpG slice?
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
N_DHS = 22_000
N_CCRE = 22_000
N_CPG = 6_000
WINDOW = 200
MAX_NSAMP = 5
EXCLUDE = {"dELS"}
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
    dhs_idx = rng.choice(pool, size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    cls_array = ccre_meta["cls"].to_numpy()
    classes = sorted([c for c in np.unique(cls_array) if c not in EXCLUDE])
    per_cls = N_CCRE // len(classes)  # 22000/7 = 3142
    parts = []
    cb_used = set()
    for c in classes:
        cpool = np.where(cls_array == c)[0]
        sel = rng.choice(cpool, size=per_cls, replace=False) if len(cpool) >= per_cls else cpool.copy()
        parts.append(sel)
        cb_used.update(sel.tolist())
        print(f"  cCRE {c}: pool={len(cpool):,}, sampled={len(sel):,}")
    ccre_idx = np.concatenate(parts)
    if len(ccre_idx) < N_CCRE:
        extra_pool = np.array([i for i in np.where(cls_array == "pELS")[0] if i not in cb_used])
        ccre_idx = np.concatenate(
            [ccre_idx, rng.choice(extra_pool, size=N_CCRE - len(ccre_idx), replace=False)]
        )
        cb_used.update(ccre_idx.tolist())
    ccre_rows = np.asarray(ccre_seqs[ccre_idx])

    cpg_pool = np.array([i for i in np.where(np.isin(cls_array, ["PLS", "pELS"]))[0] if i not in cb_used])
    print(f"CpG-rich extra pool: {len(cpg_pool):,}")
    cpg_idx = rng.choice(cpg_pool, size=N_CPG, replace=False)
    cpg_rows = np.asarray(ccre_seqs[cpg_idx])

    all_rows = np.concatenate([dhs_rows, ccre_rows, cpg_rows], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} ({N_DHS} orth + {N_CCRE} cCRE-no-dELS + {N_CPG} CpG)")


if __name__ == "__main__":
    main()
