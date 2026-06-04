"""
E014: Gosai chrom-balanced — equal sequences per chromosome.

E008 (chr 7,9,13,21,X) gave +0.013. Test whether the lever is
(a) eval-set leak through those specific chroms, or (b) chrom-balance.

Equal samples per chrom = ~2174 each across 23 chroms (1-22 + X; skip Y).
If (a): chrom-balance won't beat random Gosai (0.32), since most chroms
        aren't in eval; mixing them dilutes the signal.
If (b): chrom-balance beats E008, since it spreads representation evenly
        across the full genome.
"""
import os, numpy as np
from collections import defaultdict

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000

CHROMS = [str(i) for i in range(1, 23)] + ["X"]  # skip Y (only 19)

rng = np.random.default_rng(SEED)

print("Loading Gosai by chrom...", flush=True)
by_chrom = defaultdict(list)
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in CHROMS: continue
        if any(b not in "ACGT" for b in c[11]): continue
        by_chrom[c[1]].append(c[11])
print(f"Chroms loaded: {len(by_chrom)}", flush=True)
for ch in CHROMS:
    print(f"  chr{ch}: {len(by_chrom[ch])}", flush=True)

n_per = N_SEQS // len(CHROMS)
remainder = N_SEQS - n_per * len(CHROMS)
print(f"Sampling {n_per} per chrom (+{remainder} extra to first chroms)", flush=True)

final = []
for i, ch in enumerate(CHROMS):
    target = n_per + (1 if i < remainder else 0)
    pool = by_chrom[ch]
    idx = rng.permutation(len(pool))[:target]
    final.extend(pool[j] for j in idx)
    print(f"  chr{ch}: {target}", flush=True)

rng.shuffle(final)
assert len(final) == N_SEQS, f"Got {len(final)}"

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
