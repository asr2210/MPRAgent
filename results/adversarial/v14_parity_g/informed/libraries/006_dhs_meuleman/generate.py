#!/usr/bin/env python3
"""
006_dhs_meuleman — Direct port of dhs_random / dhs_stratified from the
informed baselines, using the Meuleman 2020 DHS Index (3.59M elements,
hg38, 16 NMF components for tissue/cell-type signal).

Goal: a faithful test of the dhs_topic / dhs_stratified baselines. If the
instructions.md table claim (r ≈ 0.71 on eval_01) is accurate, this library
should produce a measurable signal. If the score stays at v14 noise floor,
the informed baselines are not actually applicable to v14's evaluator.

Strategy: stratified sample across all 16 NMF components (n/16 per
component), 200bp centered on the DHS summit.
"""
import gzip
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
DHS = DATA / "DHS_Index_hg38.txt.gz"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

N_SEQ = 50_000
L = 200
SEED = 0
ACGT = set("ACGT")


def main():
    rng = random.Random(SEED)
    np.random.seed(SEED)

    # Group DHSs by component
    groups = defaultdict(list)
    with gzip.open(DHS, "rt") as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            summit = int(parts[6])
            component = parts[9]
            groups[component].append((chrom, summit))

    n_comp = len(groups)
    per_comp = N_SEQ // n_comp
    extra = N_SEQ - per_comp * n_comp
    print(f"Found {n_comp} NMF components, {per_comp} per (+{extra} extras to first comps)")

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

    sequences = []
    comp_keys = sorted(groups.keys())
    for i, comp in enumerate(comp_keys):
        target = per_comp + (1 if i < extra else 0)
        pool = groups[comp]
        rng.shuffle(pool)
        taken = 0
        for chrom, summit in pool:
            if taken >= target:
                break
            s = summit - L // 2
            e = s + L
            if s < 0:
                continue
            try:
                clen = len(fa[chrom])
            except KeyError:
                continue
            if e > clen:
                continue
            seq = str(fa[chrom][s:e]).upper()
            if len(seq) != L or not set(seq).issubset(ACGT):
                continue
            sequences.append(seq)
            taken += 1
        print(f"  {comp}: took {taken}/{target}")

    if len(sequences) < N_SEQ:
        # Top up from any pool
        all_remaining = [(c, s) for comp in groups for (c, s) in groups[comp]]
        rng.shuffle(all_remaining)
        for chrom, summit in all_remaining:
            if len(sequences) >= N_SEQ:
                break
            s = summit - L // 2
            e = s + L
            if s < 0:
                continue
            try:
                clen = len(fa[chrom])
            except KeyError:
                continue
            if e > clen:
                continue
            seq = str(fa[chrom][s:e]).upper()
            if len(seq) != L or not set(seq).issubset(ACGT):
                continue
            sequences.append(seq)

    rng.shuffle(sequences)  # mix components in stream
    assert len(sequences) == N_SEQ
    with open(OUT, "w") as f:
        for s in sequences:
            f.write(s + "\n")
    print(f"Wrote {N_SEQ} DHS-stratified sequences -> {OUT}")


if __name__ == "__main__":
    main()
