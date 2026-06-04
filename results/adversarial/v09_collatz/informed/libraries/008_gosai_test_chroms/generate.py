"""
E008: Gosai sequences from Malinois TEST chromosomes (7, 9, 13, 21, X).

Diagnostic for eval distribution: if prepare.py's eval uses sequences
from the Malinois test chromosomes (a natural choice — the published
held-out set), training on these same chromosomes should yield data
leakage and a measurable Pearson jump vs random Gosai (0.323).

If E008 ≈ E004 (random Gosai), the eval distribution is independent of
Malinois chromosome split, so the ceiling is structural.

50K random samples from chroms 7, 9, 13, 21, X (117K available).
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
TEST_CHROMS = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai (test chroms only)...", flush=True)
seqs = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in TEST_CHROMS: continue
        if any(b not in "ACGT" for b in c[11]): continue
        seqs.append(c[11])
print(f"Total in test chroms: {len(seqs)}", flush=True)

idx = rng.permutation(len(seqs))[:N_SEQS]
final = [seqs[i] for i in idx]
print(f"Selected {len(final)}", flush=True)
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {out}", flush=True)
