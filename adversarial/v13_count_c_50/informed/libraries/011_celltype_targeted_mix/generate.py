"""
Experiment 011 — celltype_targeted_mix

Across all my libraries so far, K562 scores 0.05–0.10 higher than HepG2
and SK-N-SH on eval_01. K562 is over-represented in the DHS pool
(erythroid + many lymphoid biosamples). HepG2 (liver carcinoma) and
SK-N-SH (neuroblastoma) cell types map to underweighted NMF components.

Design:
  25 000 cCRE class-balanced (unchanged from 008)
  25 000 DHS, oversampling components likely to be HepG2/SK-N-SH-active:
    - Neural          (5 000)   ← SK-N-SH-relevant
    - Cancer/epithelial(5 000)  ← HepG2-relevant
    - Digestive       (5 000)   ← HepG2-relevant (liver / gut)
    - Renal / cancer  (3 000)   ← HepG2-adjacent
    - Organ devel./renal(2 000)
    - Tissue invariant(2 500)   ← cross-cell
    - other 10 components, uniform split (2 500 / 10 = 250 each)

Predicts: HepG2 + SK-N-SH eval_01 scores lift; K562 holds (because cCRE
+ tissue-invariant carry the K562-relevant content); overall eval_01 up.
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
    "Neural": 5000,
    "Cancer / epithelial": 5000,
    "Digestive": 5000,
    "Renal / cancer": 3000,
    "Organ devel. / renal": 2000,
    "Tissue invariant": 2500,
}
N_OTHER = 2500          # split among the remaining 10 components

N_CCRE = 25_000
PER_CLS = N_CCRE // 8


def main():
    rng = np.random.default_rng(SEED)
    # DHS
    dhs_meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    comps = dhs_meta["component"].to_numpy()
    used_comps = set(DHS_TARGET.keys())
    other_comps = [c for c in np.unique(comps) if c not in used_comps]
    per_other = N_OTHER // len(other_comps)

    targets = dict(DHS_TARGET)
    for c in other_comps:
        targets[c] = per_other

    # Top up to 25 k
    total = sum(targets.values())
    diff = 25_000 - total
    if diff > 0:
        targets["Tissue invariant"] = targets.get("Tissue invariant", 0) + diff

    parts = []
    for c, n in targets.items():
        pool = np.where(comps == c)[0]
        if n <= 0:
            continue
        if len(pool) >= n:
            sel = rng.choice(pool, size=n, replace=False)
        else:
            sel = pool.copy()
        parts.append(sel)
        print(f"  DHS {c}: pool={len(pool):,}, sampled={len(sel):,}")
    dhs_idx = np.concatenate(parts)
    assert dhs_idx.shape[0] == 25_000, f"dhs idx wrong: {dhs_idx.shape[0]}"
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    # cCRE class-balanced
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
