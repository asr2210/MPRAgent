"""
E021: E008 replicate with seed=2 — test seed variance.

E008 (Gosai chr 7/9/13/21/X, seed=1) gave 0.336. Was that boost real
or single-seed noise? Re-run with seed=2 to estimate the noise floor.

If E021 ≈ 0.336: chr-quirk is reproducible.
If E021 ≈ 0.323 (random Gosai): single-seed noise inflated E008.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 2  # CHANGED from E008
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
print(f"Pool: {len(seqs)}", flush=True)

idx = rng.permutation(len(seqs))[:N_SEQS]
final = [seqs[i] for i in idx]
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
