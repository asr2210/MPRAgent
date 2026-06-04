#!/usr/bin/env python3
"""
009_mpra_real — Real measured MPRA sequences from Sharpr (Ernst et al. 2016).

Sharpr-MPRA: 914,348 sequences of 145bp tested in K562/HepG2 reporter assays.
12-dim activity output. We pick the top-25k and bottom-25k by mean activity
(maximum contrast), then pad 145bp → 200bp with random ACGT flanks.

Test: do REAL measured MPRA enhancer/silencer pairs produce signal on v14?
If this also returns noise, v14 evaluator does not respond to genuine
activity-bearing sequences.
"""
import h5py
import numpy as np
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
H5 = DATA / "sharpr_train.hdf5"
OUT = HERE / "sequences_0.txt"

L = 200
SHARPR_L = 145
PAD_TOTAL = L - SHARPR_L  # 55
PAD_LEFT = PAD_TOTAL // 2  # 27
PAD_RIGHT = PAD_TOTAL - PAD_LEFT  # 28
N_HALF = 25_000
SEED = 0
ACGT = "ACGT"


def main():
    rng = random.Random(SEED)
    np_rng = np.random.default_rng(SEED)
    with h5py.File(H5, "r") as f:
        X = f["X/sequence"]  # (N, 145, 4) one-hot
        Y = f["Y/output"]    # (N, 12)
        print("Loading Y...")
        y = Y[:]
        # mean activity per sequence
        mean_act = y.mean(axis=1)
        print(f"mean activity: min={mean_act.min():.3f} max={mean_act.max():.3f}")

        # Top and bottom indices by mean activity
        order = np.argsort(mean_act)
        bot_idx = order[:N_HALF]
        top_idx = order[-N_HALF:]
        all_idx = np.concatenate([bot_idx, top_idx])
        # Sort for HDF5 fancy indexing efficiency
        all_idx_sorted = np.sort(all_idx)
        print(f"Loading {len(all_idx_sorted)} sequences from HDF5...")
        X_sel = X[all_idx_sorted]
        print(f"X_sel shape: {X_sel.shape}")

        # Decode one-hot to ACGT. Argmax along last axis.
        idx = X_sel.argmax(axis=-1)  # (N, 145)
        # Check for any rows with all-zero (no base)
        has_base = X_sel.sum(axis=-1) > 0  # (N, 145)
        # Convert
        seqs = []
        for i in range(idx.shape[0]):
            row = idx[i]
            mask = has_base[i]
            chars = [ACGT[row[j]] if mask[j] else rng.choice(ACGT) for j in range(SHARPR_L)]
            core = "".join(chars)
            # Pad with random ACGT
            left = "".join(rng.choices(ACGT, k=PAD_LEFT))
            right = "".join(rng.choices(ACGT, k=PAD_RIGHT))
            seqs.append(left + core + right)

        print(f"Decoded {len(seqs)} sequences (each {len(seqs[0])}bp)")
        rng.shuffle(seqs)
        with open(OUT, "w") as f_out:
            for s in seqs:
                f_out.write(s + "\n")
    print(f"Wrote {len(seqs)} -> {OUT}")


if __name__ == "__main__":
    main()
