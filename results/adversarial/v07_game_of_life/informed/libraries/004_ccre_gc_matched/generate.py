"""Experiment 004 — cCRE elements with GC content matched to 0.50 ± 0.05.

Sample cCRE-derived 200bp windows but only keep those with GC ∈ [0.45, 0.55].
This isolates "biological content with neutral composition" — tests whether
biology adds any lift above pure random when composition is matched.

Hypothesis: if composition is the dominant lever (theory v3 from exp 003),
this library scores ≈ gc_50 (~0.40). If biology adds value when composition
is controlled, score > 0.40.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.seqlib import extract, write_sequences, SEQLEN

CCRE_BED = os.path.join(REPO, "data", "ENCODE_cCREs_v3.bed")
OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N = 50000
GC_MIN = 0.45
GC_MAX = 0.55


def load_all_ccres():
    ccres = []
    with open(CCRE_BED) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 6:
                continue
            chrom, start, end, _, _, _ = p[:6]
            ccres.append((chrom, int(start), int(end)))
    return ccres


def gc(seq):
    return sum(1 for c in seq if c in "GC") / len(seq)


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    ccres = load_all_ccres()
    print(f"All cCREs: {len(ccres):,}")

    seqs = []
    seen = set()
    tried = 0
    rejected_gc = 0

    # Shuffle and walk; for each cCRE, try center-window, then jitters
    order = rng.permutation(len(ccres))
    for jitter_round in range(5):
        for idx in order:
            if len(seqs) >= N:
                break
            chrom, s, e = ccres[idx]
            center = (s + e) // 2
            jitter = 0 if jitter_round == 0 else int(rng.integers(-100, 101))
            seq = extract(chrom, center + jitter)
            tried += 1
            if seq is None or seq in seen:
                continue
            g = gc(seq)
            if g < GC_MIN or g > GC_MAX:
                rejected_gc += 1
                continue
            seen.add(seq)
            seqs.append(seq)
        if len(seqs) >= N:
            break

    print(f"Tried {tried:,} extractions, rejected {rejected_gc:,} for GC, collected {len(seqs):,}")
    assert len(seqs) == N, f"Need {N}, got {len(seqs)}"
    rng.shuffle(seqs)
    write_sequences(seqs, OUT)
    # Quick stats
    gcs = [gc(s) for s in seqs]
    print(f"Final GC: mean={np.mean(gcs):.3f}, std={np.std(gcs):.3f}")
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
