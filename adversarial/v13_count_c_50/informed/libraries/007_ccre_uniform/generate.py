"""
Experiment 007 — ccre_uniform

Uniform random 50 000 from the SCREEN cCRE pool (2.35 M elements), no
class balancing. Diagnostic for exp 006: tells me whether the +0.031
mean gain came from the cCRE source or from forcing class balance.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SEED = 0
N = 50_000
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def main():
    rng = np.random.default_rng(SEED)
    seqs = np.load(DATA / "ccre_seqs.npy", mmap_mode="r")
    idx = rng.choice(seqs.shape[0], size=N, replace=False)
    chosen = seqs[idx]
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} sequences from {seqs.shape[0]:,} cCRE pool")


if __name__ == "__main__":
    main()
