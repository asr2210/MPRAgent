#!/usr/bin/env python3
"""
Experiment 004 — Dirichlet(0.3) extreme-composition synthetic.

50,000 sequences, each with base frequencies ~ Dirichlet(0.3,0.3,0.3,0.3)
(concentrates near simplex corners → many sequences dominated by 1-2 bases),
then 200bp i.i.d. from that composition.

Hypothesis: per-sequence compositional variance is the dominant learnable
signal at 50k. Dirichlet(0.3) gives more extreme per-sequence compositions
than Dirichlet(1.0) (baseline = 0.0768 on eval_01), so it should push
eval_01 higher if the variance hypothesis holds.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for i in range(N_TOTAL):
        p = rng.dirichlet([ALPHA, ALPHA, ALPHA, ALPHA])
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        seqs.append("".join(ALPHABET[idx]))

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"

    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences (Dirichlet alpha={ALPHA}) to {OUT}")


if __name__ == "__main__":
    main()
