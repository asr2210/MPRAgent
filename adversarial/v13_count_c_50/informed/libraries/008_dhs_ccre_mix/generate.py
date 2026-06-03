"""
Experiment 008 — dhs_ccre_mix

25 000 uniform random DHS + 25 000 class-balanced cCRE (3125/class × 8).

Combines the two strongest sources I've found:
- DHS provides cell-type / accessibility-program coverage. It carries
  what wins on eval_07 and eval_13.
- cCRE class-balanced provides functional-class coverage (PLS, CTCF,
  promoters, distal enhancers). It carries what wins on eval_04/08/09.

If the two are independently informative, the mix should beat either
alone. If they overlap, it'll tie one of them.
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
N_DHS = 25_000
N_CCRE = 25_000
PER_CLS = N_CCRE // 8  # 3125


def main():
    rng = np.random.default_rng(SEED)
    # DHS portion
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    dhs_idx = rng.choice(dhs_seqs.shape[0], size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])
    # cCRE portion (class-balanced)
    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    parts = []
    for c in sorted(ccre_meta["cls"].unique()):
        pool = np.where(ccre_meta["cls"].to_numpy() == c)[0]
        if len(pool) >= PER_CLS:
            sel = rng.choice(pool, size=PER_CLS, replace=False)
        else:
            sel = pool.copy()
        parts.append(sel)
    ccre_idx = np.concatenate(parts)
    if len(ccre_idx) < N_CCRE:
        extra_pool = np.setdiff1d(
            np.where(ccre_meta["cls"].to_numpy() == "dELS")[0], ccre_idx
        )
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
    print(f"Wrote {N:,}: {N_DHS} DHS + {N_CCRE} cCRE class-balanced")


if __name__ == "__main__":
    main()
