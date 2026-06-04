"""Experiment 024 — 2-way hybrid: 25k natural cCREs + 25k random GC-50.

Test whether NATURAL info from cCREs combines constructively with synthetic
random sequences. 012's 4-way mix lost; cleaner 2-way blend may differ.

Prediction (theory v15): plateau, ~0.384.
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
N_CCRE = 25_000
N_RAND = 25_000


def load_ccre_seqs(n_target, rng):
    _ = fasta()
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
        if len(seqs) >= n_target:
            break
        center = (s + e) // 2
        for off in (0, -50, 50, -100, 100, -150, 150):
            seq = extract(chrom, center + off)
            if seq is None or seq in seen:
                continue
            seen.add(seq)
            seqs.append(seq)
            break
    return seqs, seen


def gen_random(n, rng, seen):
    lookup = np.array(list("ACGT"))
    out = []
    while len(out) < n:
        batch = n - len(out)
        idx = rng.integers(0, 4, size=(batch * 2, SEQLEN), dtype=np.int8)
        for row in lookup[idx]:
            s = "".join(row)
            if s in seen:
                continue
            seen.add(s)
            out.append(s)
            if len(out) >= n:
                break
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    print(f"Loading {N_CCRE:,} natural cCREs...")
    ccre_seqs, seen = load_ccre_seqs(N_CCRE, rng)
    print(f"Got {len(ccre_seqs):,} cCRE sequences")
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in ccre_seqs])
    print(f"  cCRE GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")

    print(f"Generating {N_RAND:,} GC-50 random sequences...")
    rand_seqs = gen_random(N_RAND, rng, seen)
    print(f"Got {len(rand_seqs):,} random sequences")
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in rand_seqs])
    print(f"  random GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")

    all_seqs = ccre_seqs + rand_seqs
    rng.shuffle(all_seqs)
    print(f"Total: {len(all_seqs):,} (unique={len(set(all_seqs)):,})")
    assert len(all_seqs) == 50_000

    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in all_seqs])
    print(f"Library GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")
    write_sequences(all_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
