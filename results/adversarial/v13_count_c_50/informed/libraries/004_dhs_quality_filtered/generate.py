"""
Experiment 004 — dhs_quality_filtered

Uniform-random DHS, but restricted to elements with mean_signal > 0.5
AND numsamples >= 3. Filters out the noisy single-biosample peaks that
dominate the raw pool's lower tail (1.18 M of 3.59 M elements survive).

Rationale: exp 001 showed that biasing toward 1-biosample peaks tanks
the library. The complement (excluding them) should improve information
per slot. If quality DHS noticeably beats uniform DHS, then noisy peaks
in the pool are actively hurting training.
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
    meta = pd.read_parquet(DATA / "dhs_meta.parquet")
    seqs = np.load(DATA / "dhs_seqs.npy", mmap_mode="r")

    keep = (meta["mean_signal"].to_numpy() > 0.5) & (meta["numsamples"].to_numpy() >= 3)
    pool = np.where(keep)[0]
    print(f"Pool: {len(pool):,} / {len(meta):,} after quality filter", flush=True)

    sel = rng.choice(pool, size=N, replace=False)
    chosen = seqs[sel]
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")
    print(f"Wrote {N:,} sequences", flush=True)


if __name__ == "__main__":
    main()
