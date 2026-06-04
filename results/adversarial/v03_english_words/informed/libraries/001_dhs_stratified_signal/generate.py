"""
Experiment 001: DHS sequences stratified by BOTH NMF component AND signal quartile.

Hypothesis: A library that ensures both cell-type-program coverage (16 NMF
components) AND activity dynamic-range coverage (4 signal quartiles per component)
will outperform dhs_topic. The model needs to learn the full sequence->activity
mapping, which requires seeing strong, moderate, and weak elements per regulatory
program — not just the most-cell-type-specific ones.

Source: meuleman.org pre-extracted DHS dataset, 160k 200bp sequences with NMF
component assignments and DHS signal intensity.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).parent / "sequences_0.txt"
DATA = Path(__file__).resolve().parents[2] / "data" / "train_all_classifier_light.csv.gz"
N_TOTAL = 50_000
N_PER_COMPONENT = N_TOTAL // 16  # 3125
N_PER_QUARTILE = N_PER_COMPONENT // 4  # 781
SEED = 0


def main() -> None:
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA, sep="\t")
    assert (df["raw_sequence"].str.len() == 200).all()

    parts = []
    for comp in sorted(df["component"].unique()):
        sub = df[df["component"] == comp].copy()
        sub["q"] = pd.qcut(sub["total_signal"], q=4, labels=False, duplicates="drop")
        for q in sorted(sub["q"].dropna().unique()):
            pool = sub[sub["q"] == q]
            take = min(N_PER_QUARTILE, len(pool))
            picked = pool.sample(n=take, random_state=int(rng.integers(0, 2**31 - 1)))
            parts.append(picked)

    chosen = pd.concat(parts, ignore_index=True)

    # If we're short due to integer division (50000 - 16*4*781 = 16),
    # top up by drawing more from the highest-signal pool of random components.
    deficit = N_TOTAL - len(chosen)
    if deficit > 0:
        remaining = df.drop(index=chosen.index, errors="ignore")
        topup = remaining.sample(n=deficit, random_state=SEED)
        chosen = pd.concat([chosen, topup], ignore_index=True)

    # Shuffle so order doesn't encode stratification (just hygiene).
    chosen = chosen.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    seqs = chosen["raw_sequence"].tolist()

    assert len(seqs) == N_TOTAL, f"Expected {N_TOTAL}, got {len(seqs)}"
    for s in seqs:
        assert len(s) == 200, f"Sequence not 200bp: len={len(s)}"
        assert set(s).issubset(set("ACGT")), f"Non-ACGT chars: {set(s) - set('ACGT')}"

    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")

    print(f"Wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
