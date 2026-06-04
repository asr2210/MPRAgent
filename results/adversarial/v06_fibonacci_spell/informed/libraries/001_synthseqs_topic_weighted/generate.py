"""
Experiment 001: SynthSeqs sampled with natural NMF-component weights.

Mirrors the dhs_topic baseline (0.7232 eval_01) but uses the
Meuleman SynthSeqs curated set (160k pre-filtered 200bp sequences,
10k per NMF component, biased toward strong topic dominance).

Sampling weight per sequence = natural frequency of its component
in the full 3.59M DHS index. Components common genome-wide
(Primitive/embryonic 17.44%, Neural 12.85%) are sampled more often.
"""
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000

# integer (synthseqs `component`) -> natural fraction in 3.59M DHS index
COMPONENT_FREQ = {
    1: 0.0439, 2: 0.0156, 3: 0.1744, 4: 0.1127,
    5: 0.0780, 6: 0.0401, 7: 0.0738, 8: 0.1285,
    9: 0.0331, 10: 0.0443, 11: 0.0268, 12: 0.0604,
    13: 0.0403, 14: 0.0236, 15: 0.0520, 16: 0.0525,
}


def main():
    rng = np.random.default_rng(SEED)

    syn = pd.read_csv(DATA / "train_synthseqs.csv.gz", sep="\t",
                      compression="gzip")
    assert syn["raw_sequence"].str.len().min() == 200
    assert syn["raw_sequence"].str.len().max() == 200

    weights = syn["component"].map(COMPONENT_FREQ).to_numpy()
    weights = weights / weights.sum()

    idx = rng.choice(len(syn), size=N, replace=False, p=weights)
    seqs = syn["raw_sequence"].iloc[idx].tolist()

    assert len(seqs) == N
    for s in seqs:
        assert len(s) == 200
        assert set(s).issubset(set("ACGT"))

    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} seqs to {OUT}")


if __name__ == "__main__":
    main()
