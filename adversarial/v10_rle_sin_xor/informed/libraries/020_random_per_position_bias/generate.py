"""
Experiment 020: random with per-position base frequency bias.

Each position k has its own multinomial drawn from Dirichlet(1,1,1,1). Sequences
sampled IID per position from these distributions. Introduces "structured"
randomness: per-position bias rather than uniform.

Tests if positional structure (a la PWM/profile) helps the oracle.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))
# Sample per-position multinomial probabilities (Dirichlet(1,1,1,1) = uniform on simplex)
probs = rng.dirichlet(np.ones(4), size=200)  # (200, 4)

with open(OUT, "w") as f:
    for _ in range(50_000):
        seq = "".join(bases[rng.choice(4, p=probs[k])] for k in range(200))
        f.write(seq + "\n")
print("Done")
