"""
E6: dinuc_markov_natural

50k 200bp sequences from a SYMMETRIC 1st-order Markov chain with
CpG depletion and mild AT depletion biases. Symmetric T => uniform
stationary distribution => exact 50% GC marginal.

Transition matrix (rows: current base, cols: next base):
       A      C      G      T
  A  0.25   0.275  0.275  0.20
  C  0.275  0.35   0.10   0.275
  G  0.275  0.10   0.35   0.275
  T  0.20   0.275  0.275  0.25

CpG (C→G and G→C) depleted to 0.10 (vs uniform 0.25).
AT (A→T and T→A) depleted to 0.20.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

BASES = np.array(["A", "C", "G", "T"])

T = np.array([
    [0.25, 0.275, 0.275, 0.20],
    [0.275, 0.35, 0.10, 0.275],
    [0.275, 0.10, 0.35, 0.275],
    [0.20, 0.275, 0.275, 0.25],
])

# Verify rows sum to 1
assert np.allclose(T.sum(axis=1), 1.0)
# Verify symmetric => uniform stationary
assert np.allclose(T, T.T)

# Stationary distribution
evals, evecs = np.linalg.eig(T.T)
i = np.argmin(np.abs(evals - 1.0))
pi = np.real(evecs[:, i]); pi = pi / pi.sum()
print(f"Stationary (A,C,G,T): {pi}")
print(f"GC fraction at stationarity: {pi[1] + pi[2]:.4f}")

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    # vectorize for speed: precompute cumulative T rows
    Tcum = np.cumsum(T, axis=1)
    with open(out_path, "w") as f:
        # initial bases from stationary (uniform)
        first = rng.integers(0, 4, size=N_SEQS, dtype=np.uint8)
        # for speed, generate all U values up front
        for n in range(N_SEQS):
            seq = np.empty(SEQ_LEN, dtype=np.uint8)
            cur = first[n]
            seq[0] = cur
            u = rng.random(SEQ_LEN - 1)
            for k in range(1, SEQ_LEN):
                cur = int(np.searchsorted(Tcum[cur], u[k - 1]))
                seq[k] = cur
            f.write("".join(BASES[seq]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
