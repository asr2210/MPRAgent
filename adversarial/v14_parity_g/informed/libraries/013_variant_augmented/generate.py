#!/usr/bin/env python3
"""
013_variant_augmented — Take Sharpr-MPRA top 5,000 high-activity sequences
and create 10 single-substitution variants each (50k total).

Rationale: variant-effect prediction is a primary deliverable of any
MPRA model. Training on dense local neighborhoods of real high-activity
sequences gives the model many comparable inputs with controlled local
differences — exactly the signal a saturation-mutagenesis benchmark needs.
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
N_PARENT = 5_000
N_VAR_PER = 10
SHARPR_L = 145
PAD_LEFT = 27
PAD_RIGHT = 28
SEED = 0
ACGT_STR = "ACGT"


def main():
    rng = random.Random(SEED)
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)[::-1]
        top_idx = order[:N_PARENT]
        sorted_idx = np.sort(top_idx)
        X_sel = f["X/sequence"][sorted_idx]

    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0
    out = []
    for i in range(idx.shape[0]):
        chars = [ACGT_STR[idx[i, j]] if has_base[i, j] else rng.choice(ACGT_STR)
                 for j in range(SHARPR_L)]
        left = "".join(rng.choices(ACGT_STR, k=PAD_LEFT))
        right = "".join(rng.choices(ACGT_STR, k=PAD_RIGHT))
        parent = left + "".join(chars) + right
        out.append(parent)
        # Generate N_VAR_PER - 1 variants
        parent_list = list(parent)
        for _ in range(N_VAR_PER - 1):
            v = list(parent_list)
            # 1-3 random substitutions
            n_sub = rng.randint(1, 3)
            positions = rng.sample(range(L), n_sub)
            for p in positions:
                alt_choices = [c for c in ACGT_STR if c != v[p]]
                v[p] = rng.choice(alt_choices)
            out.append("".join(v))

    assert len(out) == N_PARENT * N_VAR_PER, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} variant-augmented sequences -> {OUT}")


if __name__ == "__main__":
    main()
