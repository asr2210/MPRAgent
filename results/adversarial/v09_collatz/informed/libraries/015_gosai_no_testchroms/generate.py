"""
E015: Gosai EXCLUDING test chroms (7, 9, 13, 21, X).

Direct converse of E008. If E008's +0.013 boost came from eval-set leak
through chr 7/9/13/21/X, then EXCLUDING those chroms should
*underperform* random Gosai (~0.32).

If score stays near 0.32, the chrom effect is symmetric (chr-specific
sequence statistics, not literal eval overlap).
If score drops to ~0.30 or lower, eval-set leak confirmed.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
EXCLUDE = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai (excluding test chroms)...", flush=True)
seqs = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] in EXCLUDE: continue
        if any(b not in "ACGT" for b in c[11]): continue
        seqs.append(c[11])
print(f"Pool: {len(seqs)}", flush=True)

idx = rng.permutation(len(seqs))[:N_SEQS]
final = [seqs[i] for i in idx]
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
