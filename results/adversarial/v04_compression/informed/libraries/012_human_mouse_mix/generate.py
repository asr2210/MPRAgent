#!/usr/bin/env python3
"""
012_human_mouse_mix — 25k random 200bp tiles from human chr1/17/19/22
+ 25k from mouse chr1/11/17/19. Tests whether cross-species mixing
adds informativeness beyond pure human DNA.

Context:
- E2 (pure human): 0.4992
- E11 (pure mouse): 0.4485
- Linear average expectation: 0.4739
- This experiment: > 0.4992 if diversity ADDS; ~ 0.49 if it's neutral;
  < 0.4992 if dilution dominates.

seed=0. Sequences shuffled within file. Each tile from species-specific
length-weighted chromosome choice.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N_HUMAN = 25_000
N_MOUSE = 25_000
L = 200
HUMAN_CHROMS = ["chr1", "chr17", "chr19", "chr22"]
MOUSE_CHROMS = ["chr1", "chr11", "chr17", "chr19"]
DATA = Path(__file__).resolve().parents[2] / "data"


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def sample_tiles(rng, chrom_seqs, n, L):
    lengths = np.array([len(chrom_seqs[c]) for c in chrom_seqs], dtype=float)
    weights = lengths / lengths.sum()
    chroms_list = list(chrom_seqs.keys())
    valid = set("ACGT")
    out = []
    attempts = 0
    while len(out) < n and attempts < n * 10:
        attempts += 1
        c = chroms_list[rng.choice(len(chroms_list), p=weights)]
        seq = chrom_seqs[c]
        start = int(rng.integers(0, len(seq) - L))
        w = seq[start : start + L]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    return out, attempts


def main():
    rng = np.random.default_rng(SEED)
    print("loading human chroms...")
    human = {c: load_chrom(DATA / f"{c}.fa.gz") for c in HUMAN_CHROMS}
    print(f"  {sum(len(s) for s in human.values()):,} bp")
    print("loading mouse chroms...")
    mouse = {c: load_chrom(DATA / "mm10" / f"{c}.fa.gz") for c in MOUSE_CHROMS}
    print(f"  {sum(len(s) for s in mouse.values()):,} bp")

    human_tiles, ha = sample_tiles(rng, human, N_HUMAN, L)
    print(f"human: {len(human_tiles)} of {ha} attempts")
    mouse_tiles, ma = sample_tiles(rng, mouse, N_MOUSE, L)
    print(f"mouse: {len(mouse_tiles)} of {ma} attempts")

    out = human_tiles + mouse_tiles
    assert len(out) == N_HUMAN + N_MOUSE
    rng.shuffle(out)

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")
    print(f"wrote {len(out)} → {out_path}")


if __name__ == "__main__":
    main()
