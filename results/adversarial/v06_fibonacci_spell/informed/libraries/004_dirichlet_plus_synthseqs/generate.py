"""
Experiment 004: 50% Dirichlet synthetic + 50% real DHS SynthSeqs (topic-weighted).

Tests the central question: does real biological sequence add value when
combined with the best synthetic diversity strategy?
- Dirichlet alone (exp 002) → 0.1395
- SynthSeqs alone (exp 001) → 0.1319

If mix > 0.1395: real biology contributes information beyond synthetic diversity.
If mix < 0.1395: synthetic diversity dominates; real sequences are dead weight.
"""
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N_TOTAL = 50_000
N_DIRI = N_TOTAL // 2
N_DHS = N_TOTAL - N_DIRI
LEN = 200
BASES = np.array(list("ACGT"))

COMPONENT_FREQ = {
    1: 0.0439, 2: 0.0156, 3: 0.1744, 4: 0.1127,
    5: 0.0780, 6: 0.0401, 7: 0.0738, 8: 0.1285,
    9: 0.0331, 10: 0.0443, 11: 0.0268, 12: 0.0604,
    13: 0.0403, 14: 0.0236, 15: 0.0520, 16: 0.0525,
}


def main():
    rng = np.random.default_rng(SEED)

    # Dirichlet half
    probs = rng.dirichlet((0.5, 0.5, 0.5, 0.5), size=N_DIRI)
    diri_seqs = []
    for i in range(N_DIRI):
        b = rng.choice(4, size=LEN, p=probs[i])
        diri_seqs.append("".join(BASES[b]))

    # DHS half: topic-weighted sampling from synthseqs
    syn = pd.read_csv(DATA / "train_synthseqs.csv.gz", sep="\t",
                      compression="gzip")
    weights = syn["component"].map(COMPONENT_FREQ).to_numpy()
    weights = weights / weights.sum()
    idx = rng.choice(len(syn), size=N_DHS, replace=False, p=weights)
    dhs_seqs = syn["raw_sequence"].iloc[idx].tolist()

    seqs = diri_seqs + dhs_seqs
    rng.shuffle(seqs)
    assert len(seqs) == N_TOTAL
    for s in seqs:
        assert len(s) == 200
        assert set(s).issubset(set("ACGT"))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N_TOTAL} seqs ({N_DIRI} dirichlet + {N_DHS} DHS) to {OUT}")


if __name__ == "__main__":
    main()
