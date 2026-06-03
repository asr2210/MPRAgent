"""
Experiment 005 — dhs_component_stratified

3 125 DHS elements per NMF component (16 components × 3 125 = 50 000),
uniform random within each. Forces equal representation across all 16
chromatin accessibility programs, overriding the raw pool's imbalance
(Primitive/embryonic ≈ 17.4 % of the pool, Stromal A only 1.6 %).

Direct test of the post-exp-004 "diversity dominates per-element
quality" hypothesis. Predicts a small gain over exp 002 (uniform DHS).
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
PER_COMP = N // 16  # 3125
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def main():
    rng = np.random.default_rng(SEED)
    meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    comps = meta["component"].to_numpy()
    unique = sorted(np.unique(comps))
    print(f"{len(unique)} components; {PER_COMP} per component")
    assert len(unique) == 16, f"expected 16, got {len(unique)}"

    all_idx = []
    for c in unique:
        pool = np.where(comps == c)[0]
        sel = rng.choice(pool, size=PER_COMP, replace=False)
        all_idx.append(sel)
        print(f"  {c}: pool={len(pool):,}, sampled={PER_COMP}")
    idx = np.concatenate(all_idx)
    rng.shuffle(idx)
    chosen = seqs[idx]
    assert chosen.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")


if __name__ == "__main__":
    main()
