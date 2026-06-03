"""
Experiment 017 — orth-DHS component-targeted + cCRE class-balanced

Combine two proven priors:
- 015's orthogonality (DHS not within 200bp of any cCRE midpoint)
- 011's HepG2/SK-N-SH-relevant component targeting

Design: 25 k orth-DHS with component targeting + 25 k cCRE class-bal.
DHS components weighted as in 011 (Neural, Cancer/epi, Digestive, etc.)
but drawn from the orthogonal pool (~950k).

Tests whether the two priors compose additively. If 017 > 015 (0.5736),
orthogonality + cell-type targeting are independent levers.
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
N_OTHER = 2500


def non_ccre_dhs_mask(dhs_meta: pd.DataFrame, ccre_meta: pd.DataFrame, window: int) -> np.ndarray:
    mask = np.ones(len(dhs_meta), dtype=bool)
    ccre_by_chrom = {c: np.sort(g["mid"].to_numpy()) for c, g in ccre_meta.groupby("chrom")}
    dhs_chrom = dhs_meta["chrom"].to_numpy()
    dhs_summit = dhs_meta["summit"].to_numpy()
    for i, (chrom, summit) in enumerate(zip(dhs_chrom, dhs_summit)):
        mids = ccre_by_chrom.get(chrom)
        if mids is None:
            continue
        j = np.searchsorted(mids, summit)
        if j < len(mids) and mids[j] - summit < window:
            mask[i] = False
            continue
        if j > 0 and summit - mids[j-1] < window:
            mask[i] = False
    return mask


def main():
    rng = np.random.default_rng(SEED)
    dhs_meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")

    print("Computing orthogonal DHS mask...")
    mask = non_ccre_dhs_mask(dhs_meta, ccre_meta, WINDOW)
    pool_idx = np.where(mask)[0]
    print(f"Orthogonal DHS pool: {len(pool_idx):,}")

    comps_full = dhs_meta["component"].to_numpy()
    comps_orth = comps_full[pool_idx]  # components of orthogonal pool
    used_comps = set(DHS_TARGET.keys())
    other_comps = [c for c in np.unique(comps_full) if c not in used_comps]
    per_other = N_OTHER // len(other_comps)
    targets = dict(DHS_TARGET)
    for c in other_comps:
        targets[c] = per_other
    total = sum(targets.values())
    if total < N_DHS:
        targets["Tissue invariant"] += N_DHS - total

    parts = []
    for c, n in targets.items():
        # restrict pool to (orthogonal AND component == c)
        in_comp = pool_idx[comps_orth == c]
        n = min(n, len(in_comp))
        sel = rng.choice(in_comp, size=n, replace=False) if n > 0 else np.array([], dtype=int)
        parts.append(sel)
        print(f"  DHS {c}: orth-pool={len(in_comp):,}, sampled={len(sel):,}")
    dhs_idx = np.concatenate(parts)[:N_DHS]
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])
    print(f"DHS total: {len(dhs_idx)}")

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
    assert all_rows.shape == (N, 200), f"shape {all_rows.shape}"
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,}")


if __name__ == "__main__":
    main()
