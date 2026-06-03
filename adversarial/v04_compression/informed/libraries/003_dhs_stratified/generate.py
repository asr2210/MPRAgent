#!/usr/bin/env python3
"""
003_dhs_stratified — DHS sequences balanced across 16 NMF components.

Source: Meuleman SynthSeqs training set (160k sequences, 10k per NMF
component, chr3–chrY). I sample 3,125 per component (50k total) without
replacement, seed=0.

This is the published `dhs_stratified` strategy: each chromatin
accessibility program contributes equally regardless of its size in
the original DHS index.

Compare to 002_genome_random (0.4992): if E3 substantially beats it,
the DHS accessibility annotation is adding value on top of "real DNA".
If E3 only matches, real DNA per se is most of the signal.

Generalization argument: the 16 NMF components are a low-rank
factorization of chromatin accessibility across 733 biosamples. By
sampling equally across components, the library forces equal exposure
to all cell-lineage-defined regulatory programs — not just the
programs active in K562/HepG2/SKNSH. A model trained on this should
be better equipped to predict activity in unseen cell types whose
accessibility profile resembles any of these 16 programs.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N = 50_000
PER_COMPONENT = N // 16  # 3125; 16 * 3125 = 50000
DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv(DATA_DIR / "synthseqs_train.csv.gz", sep="\t")
    assert (df["raw_sequence"].str.len() == 200).all()
    components = sorted(df["component"].unique())
    assert len(components) == 16, components
    assert N == PER_COMPONENT * 16

    picks = []
    for c in components:
        sub = df[df["component"] == c]
        idx = rng.choice(len(sub), size=PER_COMPONENT, replace=False)
        picks.append(sub.iloc[idx]["raw_sequence"].values)
    seqs = np.concatenate(picks)
    rng.shuffle(seqs)
    assert len(seqs) == N

    valid = set("ACGT")
    bad = sum(1 for s in seqs if not set(s).issubset(valid))
    print(f"non-ACGT sequences: {bad} (will need fallback if >0)")
    assert bad == 0, f"need to handle masked bases ({bad} sequences)"

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s)
            f.write("\n")
    print(f"wrote {len(seqs)} sequences → {out_path}")


if __name__ == "__main__":
    main()
