#!/usr/bin/env python3
"""
014_revcomp_augmented — Sharpr top+bottom 12,500 each + their reverse
complements. Tests whether strand-augmented training helps the model
learn strand-invariant motif representations.

Rationale: convolutional MPRA models often benefit from reverse-complement
augmentation. Including both strands of each parent sequence doubles motif
exposure without changing biology.
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
N_PARENT_HALF = 12_500
SHARPR_L = 145
PAD_LEFT = 27
PAD_RIGHT = 28
SEED = 0
ACGT_STR = "ACGT"
RC = str.maketrans("ACGT", "TGCA")


def main():
    rng = random.Random(SEED)
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)
        bot_idx = order[:N_PARENT_HALF]
        top_idx = order[-N_PARENT_HALF:]
        all_idx = np.concatenate([bot_idx, top_idx])
        sorted_idx = np.sort(all_idx)
        X_sel = f["X/sequence"][sorted_idx]

    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0
    parents = []
    for i in range(idx.shape[0]):
        chars = [ACGT_STR[idx[i, j]] if has_base[i, j] else rng.choice(ACGT_STR)
                 for j in range(SHARPR_L)]
        left = "".join(rng.choices(ACGT_STR, k=PAD_LEFT))
        right = "".join(rng.choices(ACGT_STR, k=PAD_RIGHT))
        parents.append(left + "".join(chars) + right)

    out = []
    for p in parents:
        out.append(p)
        out.append(p.translate(RC)[::-1])

    assert len(out) == 50_000, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} (parent+RC) sequences -> {OUT}")


if __name__ == "__main__":
    main()
