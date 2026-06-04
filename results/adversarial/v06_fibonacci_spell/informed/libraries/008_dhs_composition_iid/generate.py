"""
Experiment 008: real-DHS composition-matched iid sequences.

Take 50,000 SynthSeqs (topic-weighted). For each, compute its per-sequence
base composition (pA, pC, pG, pT). Replace the sequence with 200 iid
bases drawn from that composition.

This decomposes the gap between real DHS (exp 001 = 0.1319) and
dirichlet (exp 002 = 0.1395):
- If 008 ≈ 0.1319 → real DHS composition distribution is the bottleneck
  (not enough diversity); structural patterns within real seqs don't matter
- If 008 ≈ 0.1395 → real DHS composition is fine, but structure HURTS;
  removing structure recovers the dirichlet score
- If 008 < 0.1319 → composition WAS doing some work in real DHS;
  iid version is even worse
"""
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
BASES = np.array(list("ACGT"))

COMPONENT_FREQ = {
    1: 0.0439, 2: 0.0156, 3: 0.1744, 4: 0.1127,
    5: 0.0780, 6: 0.0401, 7: 0.0738, 8: 0.1285,
    9: 0.0331, 10: 0.0443, 11: 0.0268, 12: 0.0604,
    13: 0.0403, 14: 0.0236, 15: 0.0520, 16: 0.0525,
}
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def main():
    rng = np.random.default_rng(SEED)

    syn = pd.read_csv(DATA / "train_synthseqs.csv.gz", sep="\t",
                      compression="gzip")
    weights = syn["component"].map(COMPONENT_FREQ).to_numpy()
    weights = weights / weights.sum()
    idx = rng.choice(len(syn), size=N, replace=False, p=weights)

    # Vectorized count of A/C/G/T per chosen sequence
    chosen = syn["raw_sequence"].iloc[idx].to_numpy()
    counts = np.zeros((N, 4), dtype=np.int32)
    for i, s in enumerate(chosen):
        a = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        # A=65, C=67, G=71, T=84
        counts[i, 0] = (a == 65).sum()
        counts[i, 1] = (a == 67).sum()
        counts[i, 2] = (a == 71).sum()
        counts[i, 3] = (a == 84).sum()
    probs = counts / counts.sum(axis=1, keepdims=True)

    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} composition-matched iid seqs to {OUT}")


if __name__ == "__main__":
    main()
