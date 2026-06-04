"""
Experiment 010 — dhs_ccre_mix_rc_aug

Take the exp 008 mix design (DHS + cCRE class-balanced) but cut element
count in half (25 000 unique) and augment with reverse-complement of
each, for a final library of 50 000 sequences.

Tests whether the regulatory-grammar model benefits from explicit
forward+reverse-strand exposure (TFs bind on both strands; the model
shouldn't have to learn strand symmetry implicitly).

If RC-augmented beats 008 (0.5671), strand symmetry was a missing prior.
If it loses, element diversity is doing more work than augmentation.
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
N_UNIQUE = 25_000
N_DHS = 12_500
N_CCRE = 12_500
PER_CLS = N_CCRE // 8  # 1562, but we'll top up
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"

# A=65 C=67 G=71 T=84
COMP = np.zeros(256, dtype=np.uint8)
COMP[65] = 84
COMP[67] = 71
COMP[71] = 67
COMP[84] = 65


def reverse_complement(arr: np.ndarray) -> np.ndarray:
    return COMP[arr][:, ::-1]


def main():
    rng = np.random.default_rng(SEED)
    # DHS
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    dhs_idx = rng.choice(dhs_seqs.shape[0], size=N_DHS, replace=False)
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

    unique = np.concatenate([dhs_rows, ccre_rows], axis=0)
    assert unique.shape == (N_UNIQUE, 200)
    rc = reverse_complement(unique)
    all_rows = np.concatenate([unique, rc], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} ({N_UNIQUE} unique × 2 strands)")


if __name__ == "__main__":
    main()
