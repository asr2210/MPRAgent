"""
Experiment 012: SynthSeqs with topic-UNIFORM sampling.

Exp 001 sampled synthseqs weighted by natural NMF-component frequency
in the 3.59M DHS index (over-weighting Primitive/embryonic 17.4% and
Neural 12.9%, under-weighting minority components like Stromal A 1.6%).
Result: eval_01 = 0.1319.

This experiment samples UNIFORMLY across components: ~3125 sequences
from each of 16 components. Hypothesis: lifts performance on minority-
component cell types, especially K562 (erythroid; component 15 is
Myeloid/erythroid). Across all libraries, K562 head correlation is
consistently ~0.04 — the bottleneck on mean_r. Even modest K562 lift
(0.04 -> 0.08) raises mean by ~0.013.

Cleaner test of real-DHS sequences than 001 because topic balance
removes the over/under-sampling confounder.
"""
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
N_COMPONENTS = 16


def main():
    rng = np.random.default_rng(SEED)

    syn = pd.read_csv(DATA / "train_synthseqs.csv.gz", sep="\t",
                      compression="gzip")
    assert syn["raw_sequence"].str.len().min() == 200

    per_component = N // N_COMPONENTS  # 3125
    extra = N - per_component * N_COMPONENTS  # 0

    selected = []
    for c in range(1, N_COMPONENTS + 1):
        pool = syn.index[syn["component"] == c].to_numpy()
        take = per_component + (1 if (c - 1) < extra else 0)
        # synthseqs has 10k per component, so 3125 << 10k, no replacement issues
        chosen = rng.choice(pool, size=take, replace=False)
        selected.append(chosen)
    selected = np.concatenate(selected)
    rng.shuffle(selected)
    seqs = syn["raw_sequence"].iloc[selected].tolist()

    assert len(seqs) == N
    for s in seqs:
        assert len(s) == 200
        assert set(s).issubset(set("ACGT"))

    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} topic-uniform synthseqs to {OUT}")


if __name__ == "__main__":
    main()
