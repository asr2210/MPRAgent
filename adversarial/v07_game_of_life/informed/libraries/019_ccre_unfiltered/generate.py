"""Experiment 019 — Real cCREs, NO GC filter (natural composition variation).

All cCRE-based libraries to date filter to GC ∈ [0.45, 0.55]. This one lets
the natural cCRE composition flow through to test whether the per-sequence
composition variation in REAL enhancers helps or hurts.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.seqlib import extract, write_sequences, fasta, SEQLEN

CCRE_BED = os.path.join(REPO, "data", "ENCODE_cCREs_v3.bed")
OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N_TARGET = 50_000


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    _ = fasta()  # warm up

    print(f"Loading cCREs from {CCRE_BED}...")
    ccres = []
    with open(CCRE_BED) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 3:
                continue
            try:
                ccres.append((p[0], int(p[1]), int(p[2])))
            except ValueError:
                continue
    rng.shuffle(ccres)
    print(f"Have {len(ccres):,} cCREs")

    seqs = []
    seen = set()
    for chrom, s, e in ccres:
        if len(seqs) >= N_TARGET:
            break
        center = (s + e) // 2
        # try center and small offsets
        for off in (0, -50, 50, -100, 100, -150, 150):
            seq = extract(chrom, center + off)
            if seq is None or seq in seen:
                continue
            seen.add(seq)
            seqs.append(seq)
            break

    print(f"Got {len(seqs):,} sequences")
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}, "
          f"min={gcs.min():.3f}, max={gcs.max():.3f}, "
          f"P10={np.percentile(gcs, 10):.3f}, P90={np.percentile(gcs, 90):.3f}")

    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
