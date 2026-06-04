"""
Experiment 013 — mutation_augmented

008/011/012 are all ≈0.568 on eval_01. Element-level sampling has
plateaued. Time to test if local sequence-space densification adds
information that raw unique sampling doesn't.

Design:
  25 k unique elements (12.5 k DHS uniform + 12.5 k cCRE class-bal)
  25 k mutated versions: each = original + 5 random SNPs

Total: 50 k sequences in pairs (original, +5-SNP) across loci.

Rationale:
- TFBS motifs are small (6-12 bp). 5 random SNPs distributed across a
  200 bp window will typically *not* disrupt a motif (motif occupies
  ~6% of sequence), but will create paired examples where most of the
  regulatory grammar is preserved.
- Tests whether the model benefits from seeing "near-by" sequences
  that share most regulatory features but differ in non-coding flank.
- Analogous to data augmentation in vision (small crops/shifts).

Sanity: RC augmentation (exp 010) was a wash. Mutation augmentation
is qualitatively different — RC preserves the *exact* motif content
on the other strand; mutation breaks 5 random positions which could
hit motifs or flank.

If 013 > 008: local sequence diversity adds learning signal.
If 013 < 008: unique elements dominate.
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
PER_CLS = N_CCRE // 8
N_MUT = 5  # SNPs per mutated copy
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"

BASES = np.array([65, 67, 71, 84], dtype=np.uint8)  # A C G T


def mutate(rows: np.ndarray, n_snps: int, rng: np.random.Generator) -> np.ndarray:
    out = rows.copy()
    n, L = out.shape
    # pick n_snps positions per row (without replacement within a row)
    for i in range(n):
        positions = rng.choice(L, size=n_snps, replace=False)
        for p in positions:
            current = out[i, p]
            # pick a different base
            choices = BASES[BASES != current]
            out[i, p] = rng.choice(choices)
    return out


def main():
    rng = np.random.default_rng(SEED)
    dhs_seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    dhs_idx = rng.choice(dhs_seqs.shape[0], size=N_DHS, replace=False)
    dhs_rows = np.asarray(dhs_seqs[dhs_idx])

    ccre_meta = pd.read_parquet(DATA / "ccre_meta.parquet")
    ccre_seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    parts = []
    for c in sorted(ccre_meta["cls"].unique()):
        pool = np.where(ccre_meta["cls"].to_numpy() == c)[0]
        sel = rng.choice(pool, size=PER_CLS, replace=False) if len(pool) >= PER_CLS else pool.copy()
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
    print(f"Unique: {N_UNIQUE}")
    mut = mutate(unique, n_snps=N_MUT, rng=rng)
    print(f"Mutated: {N_UNIQUE} ({N_MUT} SNPs each)")
    all_rows = np.concatenate([unique, mut], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} ({N_UNIQUE} unique + {N_UNIQUE} mutated)")


if __name__ == "__main__":
    main()
