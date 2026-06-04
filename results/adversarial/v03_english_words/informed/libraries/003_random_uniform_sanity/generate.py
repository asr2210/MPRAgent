"""
Experiment 003: Random uniform ACGT sanity check.

Purpose: validate that my prepare.py invocation is consistent with the baseline
table's `synth_oracle` strategy. Expected eval_01 ≈ 0.684 (the published
baseline). If I get something close, my pipeline is fine and exp 001/002
failures are real reflections of the synthseqs pool's unsuitability. If I get
something far off, I have a setup issue to debug.
"""

from pathlib import Path
import numpy as np

OUT = Path(__file__).parent / "sequences_0.txt"
N_TOTAL = 50_000
LEN = 200
SEED = 0


def main() -> None:
    rng = np.random.default_rng(SEED)
    bases = np.array(list("ACGT"))
    arr = rng.integers(0, 4, size=(N_TOTAL, LEN))
    seqs = ["".join(bases[row]) for row in arr]

    assert len(seqs) == N_TOTAL
    for s in seqs:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")

    print(f"Wrote {len(seqs)} random sequences to {OUT}")


if __name__ == "__main__":
    main()
