"""
Experiment 003 — dhs80_synth20

80 % uniform-random DHS + 20 % i.i.d. uniform random synthetic sequences.

Rationale: in exp 002, eval_08 collapsed to 0.17 while every other eval
sat at 0.52+. eval_08 looks synthetic-OOD. The cheapest possible test of
"does any non-genomic content unlock eval_08" is to swap 10k of the 50k
DHS slots for fully random sequences. If eval_08 jumps meaningfully and
eval_01 holds, the model gains a lot from a small extrapolation set.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SEED = 0
N = 50_000
N_SYNTH = 10_000          # 20 %
N_DHS = N - N_SYNTH       # 40 000
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def main():
    rng = np.random.default_rng(SEED)
    seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    idx = rng.choice(seqs.shape[0], size=N_DHS, replace=False)
    dhs = seqs[idx]

    # i.i.d. uniform synthetic, ASCII bytes
    bases = np.array([65, 67, 71, 84], dtype=np.uint8)  # A C G T
    synth = rng.choice(bases, size=(N_SYNTH, 200))

    all_rows = np.concatenate([dhs, synth], axis=0)
    rng.shuffle(all_rows, axis=0)
    assert all_rows.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in all_rows:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} sequences ({N_DHS:,} DHS + {N_SYNTH:,} synth)")


if __name__ == "__main__":
    main()
