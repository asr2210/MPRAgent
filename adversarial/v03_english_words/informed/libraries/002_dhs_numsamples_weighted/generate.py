"""
Experiment 002: DHS sequences weighted by `numsamples`.

Hypothesis (T1): The 001 failure was that the synthseqs pool is biased toward
cell-type-SPECIFIC elements (most active in only a few of 733 biosamples).
By sampling with weight proportional to `numsamples`, we shift toward
constitutive elements that are more likely to be active in K562/HepG2/SK-N-SH,
giving the model real labels to learn from.

Comparison: dhs_topic baseline = 0.7232 eval_01. dhs_random (uniform from full
DHS pool) = 0.7089. If `numsamples`-weighting on the 160k subset can recover
≥0.65 on eval_01, that confirms "elements active in our cell types" is the
dominant axis. Anything below ~0.55 means the pool itself is too narrow even
after re-weighting, and we need the full DHS Index.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).parent / "sequences_0.txt"
DATA = Path(__file__).resolve().parents[2] / "data" / "train_all_classifier_light.csv.gz"
N_TOTAL = 50_000
SEED = 0


def main() -> None:
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA, sep="\t")
    assert (df["raw_sequence"].str.len() == 200).all()

    # Linear weight by numsamples. Add 1 to give zero-numsamples nonzero prob
    # (none should be zero in this data, but defensive).
    w = df["numsamples"].astype(float).values + 0.0
    p = w / w.sum()

    # Sample 50k WITHOUT replacement (we have 160k pool, 50k is OK).
    idx = rng.choice(len(df), size=N_TOTAL, replace=False, p=p)
    chosen = df.iloc[idx].reset_index(drop=True)

    seqs = chosen["raw_sequence"].tolist()

    assert len(seqs) == N_TOTAL
    for s in seqs:
        assert len(s) == 200
        assert set(s).issubset(set("ACGT"))

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")

    # diagnostics
    print(f"Wrote {len(seqs)} to {OUT}")
    print(f"numsamples sampled: mean={chosen.numsamples.mean():.1f} "
          f"median={chosen.numsamples.median():.1f} "
          f"min={chosen.numsamples.min()} max={chosen.numsamples.max()}")
    print(f"proportion (NMF) sampled: mean={chosen.proportion.mean():.3f}")
    print(f"total_signal sampled: mean={chosen.total_signal.mean():.1f}")
    print(f"component spread: {chosen.component.value_counts().sort_index().tolist()}")


if __name__ == "__main__":
    main()
