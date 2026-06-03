"""
Experiment 029: Manual composition design (4 corners + 6 edges + interior).

Hand-designed composition distribution mimicking Dirichlet(0.5)'s shape
but with explicit control over how much mass goes to each archetype:
- 4 "corners": one base dominant (60% one, 13.3% each other) — 10k
- 6 "edges": two bases dominant (40% each, 10% other two) — 15k
- interior: random Dirichlet(1.0) (~uniform) — 15k
- random Dirichlet(0.5) — 10k

Then for each composition, jitter with small Dirichlet(50) noise (~iid
sample around the target) to add diversity without losing the archetype.

Tests whether explicit composition design beats Dirichlet's natural
mixture of these archetypes.
"""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
BASES = np.array(list("ACGT"))


def jitter(p_target, rng, concentration=50.0):
    """Sample around p_target with given concentration (high = tight)."""
    alpha = p_target * concentration
    return rng.dirichlet(alpha)


def main():
    rng = np.random.default_rng(SEED)
    probs = []

    # 4 corners: 10k total, 2500 each
    for i in range(4):
        p = np.full(4, 0.4 / 3)
        p[i] = 0.6
        for _ in range(2_500):
            probs.append(jitter(p, rng))

    # 6 edges: 15k total, 2500 each
    edge_pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for i, j in edge_pairs:
        p = np.full(4, 0.1)
        p[i], p[j] = 0.4, 0.4
        for _ in range(2_500):
            probs.append(jitter(p, rng))

    # interior near-uniform: 15k
    for _ in range(15_000):
        probs.append(rng.dirichlet((1.0,)*4))

    # random Dirichlet(0.5): 10k
    for _ in range(10_000):
        probs.append(rng.dirichlet((0.5,)*4))

    probs = np.array(probs)
    rng.shuffle(probs)
    assert len(probs) == N

    seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])]) for i in range(N)]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} manual-design composition seqs to {OUT}")


if __name__ == "__main__":
    main()
