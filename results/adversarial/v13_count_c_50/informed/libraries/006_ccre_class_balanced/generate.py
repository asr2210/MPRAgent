"""
Experiment 006 — ccre_class_balanced

Sample 50 000 sequences from the SCREEN cCRE registry, equal counts
across the 8 chromatin-mark classes (6 250 per class, except where
the pool is smaller, in which case we use everything available and
re-balance the remainder uniformly).

Why cCRE vs DHS:
- DHS is "any open chromatin"; cCRE is "open chromatin + filtered by
  ChIP-seq marks indicating actual regulatory function."
- cCRE class breakdown carries explicit regulatory grammar info (PLS =
  promoter, dELS = distal enhancer, CTCF, TF, etc.), where DHS just has
  empirical co-accessibility topics.

Hypothesis: classes are more directly aligned with the cell-type-
generalizable regulatory features than DHS topics, so this should
either tie or beat uniform DHS at 0.5627.
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


def main():
    rng = np.random.default_rng(SEED)
    meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    classes = sorted(meta["cls"].unique())
    target_per_class = N // len(classes)  # 6250
    print(f"Classes: {classes}; target {target_per_class}/class")

    all_idx = []
    leftover = 0
    for c in classes:
        pool = np.where(meta["cls"].to_numpy() == c)[0]
        if len(pool) < target_per_class:
            sel = pool.copy()
            leftover += target_per_class - len(pool)
        else:
            sel = rng.choice(pool, size=target_per_class, replace=False)
        all_idx.append(sel)
        print(f"  {c}: pool={len(pool):,}, sampled={len(sel):,}")

    # Top up shortfall with extra uniform draws from the largest classes
    if leftover > 0:
        print(f"Topping up {leftover} from dELS (largest pool)")
        already = set(np.concatenate(all_idx).tolist())
        pool = np.where(meta["cls"].to_numpy() == "dELS")[0]
        pool = np.array([i for i in pool if i not in already])
        extra = rng.choice(pool, size=leftover, replace=False)
        all_idx.append(extra)

    idx = np.concatenate(all_idx)
    rng.shuffle(idx)
    assert len(idx) == N, f"got {len(idx)}"
    chosen = seqs[idx]
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")


if __name__ == "__main__":
    main()
