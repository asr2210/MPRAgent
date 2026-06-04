"""Experiment 003 — Pure promoter (PLS) library.

50k 200bp windows from ENCODE cCRE PLS-class elements (Promoter-Like Signature).
PLS elements include both 'PLS' and 'PLS,CTCF-bound' (~41k total). To reach 50k,
we generate multiple windows per element where the element is long enough
(median cCRE length 286bp), with random jitter.

Hypothesis: promoter-rich training data with richer MPRA activity signal will
let the model learn better features than diluted/random sequences. If true,
score > cCRE-balanced (0.392); if not, biology truly isn't the lever.
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


def load_pls():
    pls = []
    with open(CCRE_BED) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 6:
                continue
            chrom, start, end, _, _, label = p[:6]
            if label in ("PLS", "PLS,CTCF-bound"):
                pls.append((chrom, int(start), int(end)))
    return pls


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    pls = load_pls()
    print(f"PLS elements: {len(pls):,}")

    seqs = []
    seen = set()
    order = rng.permutation(len(pls))

    # Pass 1: extract one center-window per element
    for i in order:
        chrom, s, e = pls[i]
        seq = extract(chrom, (s + e) // 2)
        if seq is None or seq in seen:
            continue
        seen.add(seq)
        seqs.append(seq)
        if len(seqs) >= N:
            break

    # Pass 2: if still under target, generate jittered windows
    if len(seqs) < N:
        for _ in range(50):  # multiple rounds with different jitter
            if len(seqs) >= N:
                break
            for i in order:
                if len(seqs) >= N:
                    break
                chrom, s, e = pls[i]
                # jitter ±50bp around center if room
                center = (s + e) // 2
                jitter = int(rng.integers(-50, 51))
                seq = extract(chrom, center + jitter)
                if seq is None or seq in seen:
                    continue
                seen.add(seq)
                seqs.append(seq)

    print(f"Collected {len(seqs)} unique sequences")
    assert len(seqs) == N, f"Need {N}, got {len(seqs)}"
    rng.shuffle(seqs)
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
