"""
Experiment 001 — dhs_signal_specificity

Sample 50,000 DHS elements with weight ∝ mean_signal / sqrt(numsamples).
Extract 200 bp centered on the summit.

Rationale: the published `dhs_topic` baseline scores eval_01 = 0.7232 and
its description ("upweights elements with strong cell-type-specific
accessibility signal") matches this weighting closely. Without the NMF
topic-loading matrix on hand, this is the cleanest proxy we can build
from the index alone, and it anchors my deltas relative to the baseline
table for the rest of the run.
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
SEQS = DATA / "dhs_seqs.npy"
META = DATA / "dhs_meta.parquet"


def main():
    rng = np.random.default_rng(SEED)
    print("Loading DHS metadata...", flush=True)
    meta = pd.read_parquet(META)
    print(f"  {len(meta):,} candidate elements", flush=True)

    # weight ∝ mean_signal / sqrt(numsamples)
    w = meta["mean_signal"].to_numpy(dtype=np.float64)
    ns = meta["numsamples"].to_numpy(dtype=np.float64)
    w = w / np.sqrt(np.maximum(ns, 1.0))
    w = w / w.sum()

    print(f"Sampling {N:,} rows without replacement...", flush=True)
    idx = rng.choice(len(meta), size=N, replace=False, p=w)

    print("Loading sequence bytes...", flush=True)
    seqs = np.load(SEQS, mmap_mode="r")
    chosen = seqs[idx]
    assert chosen.shape == (N, 200)

    print(f"Writing {OUT}...", flush=True)
    with open(OUT, "wb") as f:
        for row in chosen:
            f.write(row.tobytes())
            f.write(b"\n")
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
