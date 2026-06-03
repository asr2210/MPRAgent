"""
Experiment 029 — 025 design + small CpG slice

Steal 3k from cCRE budget, add 3k CpG-rich PLS+pELS to lift eval_08
without dragging eval_01.

Design:
  25 k orth-DHS (numsamples ≤ 5)
  22 k cCRE balanced across 6 classes (no dELS, no CA): ~3666/class
   3 k top-CpG PLS+pELS

021 (with 8-class cCRE + CpG) tied 020 at 0.5745. With 6-class pruning
that gave +0.0017, the CpG slice may now lift eval_08 ~+0.005 without
sacrificing eval_01.
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
N_CCRE = 22_000
N_CPG = 3_000
WINDOW = 200
MAX_NSAMP = 5
EXCLUDE = {"dELS", "CA"}
BYTES_ACGT = {b"A"[0]: 0, b"C"[0]: 1, b"G"[0]: 2, b"T"[0]: 3}
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


def count_cpg(rows):
    # rows is (N, 200) uint8 of ASCII bytes — count C(67) followed by G(71)
    C, G = 67, 71
    c = rows[:, :-1] == C
    g = rows[:, 1:] == G
    return (c & g).sum(axis=1)


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
    classes = sorted([c for c in np.unique(cls_array) if c not in EXCLUDE])
    per_cls = N_CCRE // len(classes)  # 3666
    parts = []
    used_ccre = set()
    for c in classes:
        cpool = np.where(cls_array == c)[0]
        sel = rng.choice(cpool, size=per_cls, replace=False)
        parts.append(sel)
        used_ccre.update(sel.tolist())
    ccre_idx = np.concatenate(parts)
    # top up to N_CCRE if shy
    if len(ccre_idx) < N_CCRE:
        extra_pool = np.array(
            [i for i in np.where(cls_array == "pELS")[0] if i not in used_ccre]
        )
        topup = rng.choice(extra_pool, size=N_CCRE - len(ccre_idx), replace=False)
        ccre_idx = np.concatenate([ccre_idx, topup])
        used_ccre.update(topup.tolist())
    ccre_rows = np.asarray(ccre_seqs[ccre_idx])

    # CpG slice: top-CpG PLS+pELS not already used
    cpg_pool = np.array(
        [i for i in np.where(np.isin(cls_array, ["PLS", "pELS"]))[0] if i not in used_ccre]
    )
    print(f"CpG candidate pool: {len(cpg_pool):,}")
    cand_rows = np.asarray(ccre_seqs[cpg_pool])
    cpg_counts = count_cpg(cand_rows)
    # take top 3*N_CPG by CpG count, then random N_CPG to add diversity
    top_n = 3 * N_CPG
    top_idx_local = np.argpartition(-cpg_counts, top_n)[:top_n]
    cpg_idx = cpg_pool[rng.choice(top_idx_local, size=N_CPG, replace=False)]
    cpg_rows = np.asarray(ccre_seqs[cpg_idx])
    print(f"CpG slice mean CpG/200bp: {count_cpg(cpg_rows).mean():.1f} "
          f"(orth-DHS: {count_cpg(dhs_rows).mean():.1f}, cCRE: {count_cpg(ccre_rows).mean():.1f})")

    all_rows = np.concatenate([dhs_rows, ccre_rows, cpg_rows], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,}")


if __name__ == "__main__":
    main()
