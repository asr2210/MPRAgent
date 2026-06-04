"""Experiment 028 — Per-seq GC sampled from natural cCRE distribution.

For each sequence, draw target_gc ~ N(0.485, 0.098) (cCRE empirical
mean/std), clip to [0.20, 0.80], then sample i.i.d. bases at that GC.

Tests composition variance separately from natural motif content:
- Tight-GC random (007): mean 0.50, std 0.035, plateau 0.386
- cCRE GC std 0.098, mean 0.485 (019): plateau 0.384
- 028: random with cCRE-LIKE per-seq GC variance, no motifs

Prediction: ≈ plateau (composition variance isn't the lever).
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)

from utils.seqlib import write_sequences, SEQLEN

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N_TARGET = 50_000
GC_MEAN = 0.485
GC_STD = 0.098
GC_MIN, GC_MAX = 0.20, 0.80


def gen_seq_at_gc(target_gc, rng):
    pA = (1 - target_gc) / 2
    pT = pA
    pG = target_gc / 2
    pC = pG
    probs = [pA, pC, pG, pT]
    idx = rng.choice(4, size=SEQLEN, p=probs)
    return "".join("ACGT"[i] for i in idx)


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    # Vector sampling
    print(f"Generating {N_TARGET:,} sequences with cCRE-like GC distribution...")
    gc_targets = np.clip(rng.normal(GC_MEAN, GC_STD, size=N_TARGET), GC_MIN, GC_MAX)

    seqs = []
    seen = set()
    for i, gc in enumerate(gc_targets):
        for _try in range(5):
            s = gen_seq_at_gc(float(gc), rng)
            if s not in seen:
                seen.add(s)
                seqs.append(s)
                break
        if (i + 1) % 10000 == 0:
            print(f"  ...{i+1:,}")

    assert len(seqs) == N_TARGET, f"got {len(seqs)}, want {N_TARGET}"
    rng.shuffle(seqs)
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}, "
          f"P10={np.percentile(gcs, 10):.3f}, P90={np.percentile(gcs, 90):.3f}")
    print(f"Unique: {len(set(seqs)):,}")
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
