"""Smoke test: 50,000 uniform random 200bp sequences."""
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent / "sequences_0.txt"
N, L = 50_000, 200
SEED = 0

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))
arr = bases[rng.integers(0, 4, size=(N, L))]
with OUT.open("w") as f:
    for row in arr:
        f.write("".join(row) + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
