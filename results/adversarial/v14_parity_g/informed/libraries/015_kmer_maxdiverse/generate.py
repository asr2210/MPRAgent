#!/usr/bin/env python3
"""
015_kmer_maxdiverse — algorithmically construct 50k sequences that maximize
unique 6-mer coverage. Greedy: each sequence is built to add the most
new 6-mers to the running set; once all 4^6 = 4,096 6-mers are seen,
revert to uniform random.

Rationale: if v14's model is k-mer based, maximum k-mer diversity should
give it the broadest feature exposure. Worst case it degenerates to random.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
K = 6
SEED = 0
ACGT = "ACGT"


def all_kmers(seq, k):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}


def main():
    rng = random.Random(SEED)
    seen = set()
    out = []
    saturated = False
    target = 4 ** K
    with open(OUT, "w") as f:
        for i in range(N_SEQ):
            if saturated:
                seq = "".join(rng.choices(ACGT, k=L))
            else:
                # Greedy: try 5 candidate seeds, keep best by new-kmer count
                best = None
                best_new = -1
                for _ in range(5):
                    cand = "".join(rng.choices(ACGT, k=L))
                    new = len(all_kmers(cand, K) - seen)
                    if new > best_new:
                        best_new = new
                        best = cand
                seq = best
                seen.update(all_kmers(seq, K))
                if len(seen) >= target:
                    saturated = True
            f.write(seq + "\n")
            out.append(seq)
    print(f"Wrote {N_SEQ}, k-mer coverage {len(seen)}/{target} -> {OUT}")


if __name__ == "__main__":
    main()
