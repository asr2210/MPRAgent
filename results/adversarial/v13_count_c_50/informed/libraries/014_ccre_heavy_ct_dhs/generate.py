"""
Experiment 014 — ccre_heavy + celltype-targeted DHS

Mosaic strategy combining proven priors:
  - 35 k cCRE class-balanced  (012 design — lifts eval_04/09)
  - 15 k DHS component-targeted (011 design scaled — lifts eval_07/13)

DHS component allocation (15 k total), scaled from 011's 25 k:
  Neural                3000  (was 5000)
  Cancer / epithelial   3000  (was 5000)
  Digestive             3000  (was 5000)
  Renal / cancer        1800  (was 3000)
  Organ devel./renal    1200  (was 2000)
  Tissue invariant      1500  (was 2500)
  other 10 components   150 each (= 1500 total)
  total ≈ 15 000

Tests whether the two diversifying priors compose additively or
whether they trade off (mutually exclusive lifts).
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

DHS_TARGET = {
    "Neural": 3000,
    "Cancer / epithelial": 3000,
    "Digestive": 3000,
    "Renal / cancer": 1800,
    "Organ devel. / renal": 1200,
    "Tissue invariant": 1500,
}
N_OTHER = 1500
N_DHS = 15_000
N_CCRE = 35_000


def main():
    rng = np.random.default_rng(SEED)
    dhs_meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    comps = dhs_meta["component"].to_numpy()
    used_comps = set(DHS_TARGET.keys())
    other_comps = [c for c in np.unique(comps) if c not in used_comps]
    per_other = N_OTHER // len(other_comps)
    targets = dict(DHS_TARGET)
    for c in other_comps:
        targets[c] = per_other
    total = sum(targets.values())
    if total < N_DHS:
        targets["Tissue invariant"] += N_DHS - total

    parts = []
    for c, n in targets.items():
        pool = np.where(comps == c)[0]
        n = min(n, len(pool))
        sel = rng.choice(pool, size=n, replace=False) if n > 0 else np.array([], dtype=int)
        parts.append(sel)
        print(f"  DHS {c}: pool={len(pool):,}, sampled={len(sel):,}")
    dhs_idx = np.concatenate(parts)[:N_DHS]
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])
    print(f"DHS total: {len(dhs_idx)}")

    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    classes = sorted(ccre_meta["cls"].unique())
    per_cls = N_CCRE // len(classes)  # 4375
    parts = []
    for c in classes:
        pool = np.where(ccre_meta["cls"].to_numpy() == c)[0]
        sel = rng.choice(pool, size=per_cls, replace=False) if len(pool) >= per_cls else pool.copy()
        parts.append(sel)
    ccre_idx = np.concatenate(parts)
    if len(ccre_idx) < N_CCRE:
        used = set(ccre_idx.tolist())
        extra_pool = np.array([i for i in np.where(ccre_meta["cls"].to_numpy() == "dELS")[0]
                               if i not in used])
        ccre_idx = np.concatenate(
            [ccre_idx, rng.choice(extra_pool, size=N_CCRE - len(ccre_idx), replace=False)]
        )
    ccre_rows = np.asarray(ccre_seqs[ccre_idx])
    print(f"cCRE total: {len(ccre_idx)}")

    all_rows = np.concatenate([dhs_rows, ccre_rows], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200), f"shape {all_rows.shape}"
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,}")


if __name__ == "__main__":
    main()
