#!/usr/bin/env python3
"""
010_dhs_random_full — 50k DHS sites sampled uniformly at random from
the full Meuleman DHS index (~3.6M elements). Extract 200bp centered
on each summit. Matches the published `dhs_random` strategy.

Tests narrowness penalty more cleanly than E3 (which was NMF-stratified):
- E3 dhs_stratified: 0.4378 (NMF-balanced, narrow)
- E10 dhs_random_full: ? (uniform from full index, broader)
- E2 random genome: 0.4992

If E10 ≈ E3: DHS-ness per se hurts (stratification was not the issue).
If E10 between E3 and E2: stratification was a small part of the
narrowness penalty.
If E10 ≈ E2: at uniform sampling rates, DHS is no worse than random
genome (the stratification of E3 / SynthSeqs over-balances rare
programs and is what hurts).
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
N = 50_000
L = 200
DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    # Load DHS index
    print("loading DHS index...")
    df = pd.read_csv(DATA_DIR / "dhs_index.txt.gz", sep="\t")
    print(f"  {len(df):,} DHS records")
    print(f"  columns: {list(df.columns)}")
    # Restrict to standard chroms
    std = set(["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"])
    df = df[df["seqname"].isin(std)].reset_index(drop=True)
    print(f"  {len(df):,} on standard chroms")

    # Uniform random sample of 60k (oversample to absorb edge rejections)
    sample_idx = rng.choice(len(df), size=60_000, replace=False)
    sample = df.iloc[sample_idx].reset_index(drop=True)

    # Load chromosomes
    print("loading chromosomes...")
    chroms_needed = sample["seqname"].unique()
    chrom_seqs = {c: load_chrom(DATA_DIR / f"{c}.fa.gz") for c in chroms_needed}
    print(f"  loaded {len(chrom_seqs)} chromosomes")

    valid = set("ACGT")
    out = []
    for _, row in sample.iterrows():
        if len(out) >= N:
            break
        summit = int(row["summit"])
        start = summit - L // 2
        end = start + L
        seq = chrom_seqs[row["seqname"]]
        if start < 0 or end > len(seq):
            continue
        w = seq[start:end]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    assert len(out) == N, len(out)

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")
    print(f"wrote {N} sequences → {out_path}")


if __name__ == "__main__":
    main()
