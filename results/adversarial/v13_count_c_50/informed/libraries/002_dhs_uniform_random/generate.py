"""
Experiment 002 — dhs_uniform_random

Pure uniform-random sample of 50,000 DHS elements, 200 bp centered on
the summit. No weighting.

Purpose: calibration. The published `dhs_random` baseline is 0.7089 on
eval_01. If I land near that, my data extraction + harness are sound and
the experiment 001 miss was a weighting bug. If I land much lower again,
something deeper is wrong with my pipeline.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SEED = 0
N = 50_000
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
DATA = ROOT / "data"


def main():
    rng = np.random.default_rng(SEED)
    seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")
    total = seqs.shape[0]
    print(f"Pool: {total:,} elements; sampling {N:,} without replacement", flush=True)
    idx = rng.choice(total, size=N, replace=False)
    chosen = seqs[idx]
    assert chosen.shape == (N, 200)
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
