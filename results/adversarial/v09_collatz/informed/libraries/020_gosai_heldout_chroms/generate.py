"""
E020: Gosai chr 7, 9, 13, 19, 21, X (add Malinois val chr19 to test set).

E008's chr 7,9,13,21,X boost (+0.013) — if it's because these chroms are
Malinois held-out (test) and Gosai's labels on them are cleaner, then
adding chr19 (typical Malinois validation chrom) should give an additional
small boost.

If E020 > E008 (0.336): held-out-chrom hypothesis supported.
If E020 ≈ E008: chr19 is in Malinois training, not held-out — E008 boost
                stays at +0.013.
If E020 < E008: chr19 dilutes the chr 7,9,13,21,X signal — argues for
                eval overlap on those specific 5 chroms only.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
CHROMS = {"7", "9", "13", "19", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai (chr 7/9/13/19/21/X)...", flush=True)
seqs = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in CHROMS: continue
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
