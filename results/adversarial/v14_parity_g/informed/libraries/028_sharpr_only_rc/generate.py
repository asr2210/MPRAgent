#!/usr/bin/env python3
"""
028_sharpr_only_rc — Ablation: only Sharpr-MPRA poles + RC.
25k Sharpr (12.5k top, 12.5k bottom) parents + 25k RC = 50k.
"""
import random
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SHARPR = ROOT / "data" / "sharpr_train.hdf5"
OUT = HERE / "sequences_0.txt"

L = 200
N_HALF = 12_500
SHARPR_L = 145
PAD_LEFT = 27
PAD_RIGHT = 28
SEED = 0
ACGT = "ACGT"
RC = str.maketrans("ACGT", "TGCA")


def main():
    rng = random.Random(SEED)
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)
        sel = np.concatenate([order[:N_HALF], order[-N_HALF:]])
        sel_sorted = np.sort(sel)
        X_sel = f["X/sequence"][sel_sorted]
    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0
    parents = []
    for i in range(idx.shape[0]):
        chars = [ACGT[idx[i, j]] if has_base[i, j] else rng.choice(ACGT)
                 for j in range(SHARPR_L)]
        left = "".join(rng.choices(ACGT, k=PAD_LEFT))
        right = "".join(rng.choices(ACGT, k=PAD_RIGHT))
        parents.append(left + "".join(chars) + right)
    out = []
    for s in parents:
        out.append(s)
        out.append(s.translate(RC)[::-1])
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} Sharpr+RC -> {OUT}")


if __name__ == "__main__":
    main()
