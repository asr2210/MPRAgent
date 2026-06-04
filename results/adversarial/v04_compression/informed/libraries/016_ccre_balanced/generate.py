#!/usr/bin/env python3
"""
016_ccre_balanced — 200bp windows centered on ENCODE cCREs, balanced
across 4 element classes:
- 12.5k Promoter-like (PLS, PLS+CTCF)
- 12.5k Proximal enhancer-like (pELS, pELS+CTCF)
- 12.5k Distal enhancer-like (dELS, dELS+CTCF)
- 12.5k Insulator-like (CTCF-only, DNase-H3K4me3)

cCRE V3/GRCh38: 1.06M elements, 150-350bp each. Centered 200bp.
seed=0. Random sample within each class.

Hypothesis: explicit element-class balance differs from NMF
stratification (E3 = 0.4378) and DHS uniform (E10 = 0.4699).
If element-class balance is more "informative" than NMF-class
balance, this could land near or above DHS-uniform.

Predicted: 0.45-0.48 (between DHS strategies and random genome).
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 0
PER_CLASS = 12_500
L = 200
DATA = Path(__file__).resolve().parents[2] / "data"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

CLASS_GROUPS = {
    "PLS":   ["PLS", "PLS,CTCF-bound"],
    "pELS":  ["pELS", "pELS,CTCF-bound"],
    "dELS":  ["dELS", "dELS,CTCF-bound"],
    "other": ["CTCF-only,CTCF-bound", "DNase-H3K4me3",
              "DNase-H3K4me3,CTCF-bound"],
}


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    print("loading cCRE BED...")
    df = pd.read_csv(DATA / "cCRE_hg38.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "dccid", "sccid", "class"])
    print(f"  total cCREs: {len(df):,}")
    df = df[df["chrom"].isin(set(CHROMS))].reset_index(drop=True)
    print(f"  on standard chroms: {len(df):,}")

    # Assign each cCRE to a group
    cls_to_group = {c: g for g, cs in CLASS_GROUPS.items() for c in cs}
    df["group"] = df["class"].map(cls_to_group)
    df = df[df["group"].notna()].reset_index(drop=True)
    print(f"  group counts:")
    for g in CLASS_GROUPS:
        n = (df["group"] == g).sum()
        print(f"    {g}: {n:,}")

    # Sample
    pieces = []
    for g in CLASS_GROUPS:
        sub = df[df["group"] == g]
        n_sample = min(PER_CLASS * 2, len(sub))  # oversample for edges
        idx = rng.choice(len(sub), size=n_sample, replace=False)
        pieces.append(sub.iloc[idx].reset_index(drop=True))
    sampled = pd.concat(pieces, ignore_index=True)
    print(f"  sampled: {len(sampled):,} (will trim to {PER_CLASS*4})")

    chroms_needed = set(sampled["chrom"].unique())
    print(f"loading {len(chroms_needed)} chromosomes...")
    chrom_seqs = {c: load_chrom(DATA / f"{c}.fa.gz") for c in chroms_needed}

    valid = set("ACGT")
    out_by_group = {g: [] for g in CLASS_GROUPS}
    for _, row in sampled.iterrows():
        g = row["group"]
        if len(out_by_group[g]) >= PER_CLASS:
            continue
        mid = (int(row["start"]) + int(row["end"])) // 2
        start = mid - L // 2
        end = start + L
        seq = chrom_seqs[row["chrom"]]
        if start < 0 or end > len(seq):
            continue
        w = seq[start:end]
        if not set(w).issubset(valid):
            continue
        out_by_group[g].append(w)

    out = []
    for g in CLASS_GROUPS:
        n = len(out_by_group[g])
        print(f"  {g}: {n}")
        assert n == PER_CLASS, f"only got {n} of {PER_CLASS} for {g}"
        out.extend(out_by_group[g])
    rng.shuffle(out)

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")
    print(f"wrote {len(out)} → {out_path}")


if __name__ == "__main__":
    main()
