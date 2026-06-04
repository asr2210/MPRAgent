#!/usr/bin/env python3
"""
Experiment 002 — DHS biological grounding + Dirichlet compositional diversity.

50,000 sequences = 25k stratified DHS (Meuleman 160k subset, 1562/component +
8 extras to round up) + 25k Dirichlet-composition synthetic (each sequence
draws per-base frequencies from Dirichlet(1,1,1,1), then samples 200 bp i.i.d.).

Hypothesis: compositional diversity (GC variance across sequences) is what
dirichlet_composition exploits in this harness. Biological grounding (DHS)
might be additive. If so, this should beat both individually.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_DHS = 25_000
N_DIRICHLET = 25_000
N_COMPONENTS = 16
SEQ_LEN = 200
ALPHABET = np.array(list("ACGT"))


def sample_dhs(rng):
    df = pd.read_csv(DATA, sep="\t")
    per = N_DHS // N_COMPONENTS  # 1562
    extras = N_DHS - per * N_COMPONENTS  # 8
    out = []
    for comp in range(1, N_COMPONENTS + 1):
        pool = df[df["component"] == comp]
        idx = rng.choice(len(pool), size=per, replace=False)
        out.append(pool.iloc[idx]["raw_sequence"].values)
    # add extras from random components
    extra_components = rng.choice(np.arange(1, N_COMPONENTS + 1), size=extras, replace=False)
    for comp in extra_components:
        pool = df[df["component"] == comp]
        idx = rng.choice(len(pool), size=1, replace=False)
        out.append(pool.iloc[idx]["raw_sequence"].values)
    seqs = np.concatenate(out)
    assert len(seqs) == N_DHS
    return seqs


def sample_dirichlet(rng):
    # Each of N_DIRICHLET sequences: draw a base-frequency vector from
    # Dirichlet(1,1,1,1), then sample SEQ_LEN i.i.d. bases.
    out = np.empty(N_DIRICHLET, dtype=object)
    for i in range(N_DIRICHLET):
        p = rng.dirichlet([1.0, 1.0, 1.0, 1.0])
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        out[i] = "".join(ALPHABET[idx])
    return out


def main():
    rng = np.random.default_rng(SEED)
    dhs = sample_dhs(rng)
    dir_ = sample_dirichlet(rng)
    seqs = np.concatenate([dhs, dir_])
    rng.shuffle(seqs)

    assert len(seqs) == N_TOTAL
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"

    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
