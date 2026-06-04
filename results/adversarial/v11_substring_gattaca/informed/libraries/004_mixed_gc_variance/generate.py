"""
E4: mixed_gc_variance

Half the library is pure IID random (uniform). Half has per-sequence
GC content drawn uniformly from [0.35, 0.65] - moderately wider than
random_uniform's natural binomial variance.

Tests v3 hypothesis: preserving and modestly widening compositional
variance helps the harder evals (esp. SK-N-SH) without hurting easy ones.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
N_RANDOM = 25_000
N_VARIED = N_SEQS - N_RANDOM

rng = np.random.default_rng(SEED)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    bases = np.array(["A", "C", "G", "T"])
    sequences = []

    # Half: pure IID random uniform
    for _ in range(N_RANDOM):
        idx = rng.integers(0, 4, size=SEQ_LEN)
        sequences.append("".join(bases[idx]))

    # Half: per-seq GC drawn uniformly from [0.35, 0.65]
    for _ in range(N_VARIED):
        gc = rng.uniform(0.35, 0.65)
        # split GC equally between G and C; AT equally between A and T
        # P(A)=P(T)=(1-gc)/2, P(C)=P(G)=gc/2
        p = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        sequences.append("".join(bases[idx]))

    # Shuffle to interleave
    rng.shuffle(sequences)

    # Validate
    assert len(sequences) == N_SEQS
    for s in sequences:
        assert len(s) == SEQ_LEN
        assert set(s).issubset({"A", "C", "G", "T"})

    with open(out_path, "w") as f:
        for s in sequences:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
