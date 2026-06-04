#!/usr/bin/env python3
"""
Experiment 020 — Exact-composition Dirichlet(0.3) sequences.

For each Dirichlet(0.3) target composition p, build a sequence with EXACTLY
round(200·p_i) of base i (adjusted to sum to 200), shuffled to random positions.
Removes multinomial sampling noise from sequence generation.
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


def exact_counts(p, total):
    """Largest-remainder method: round to integers summing to total."""
    raw = p * total
    floors = np.floor(raw).astype(int)
    remainder = total - floors.sum()
    if remainder > 0:
        frac = raw - floors
        top_idx = np.argsort(-frac)[:remainder]
        floors[top_idx] += 1
    return floors


def main():
    rng = np.random.default_rng(SEED)
    seqs = []
    for _ in range(N_TOTAL):
        p = rng.dirichlet([ALPHA] * 4)
        counts = exact_counts(p, SEQ_LEN)
        assert counts.sum() == SEQ_LEN
        # Build pool of bases
        chars = np.concatenate([np.repeat(ALPHABET[i], counts[i]) for i in range(4)])
        rng.shuffle(chars)
        seqs.append("".join(chars))
    seqs = np.array(seqs)
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(seqs))}")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
