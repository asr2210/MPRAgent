"""
E010: Reverse-complement augmentation — 25K Gosai + 25K RC of same.

Enhancer activity is largely strand-symmetric (TF binding works in both
orientations for most factors). Including RC versions of the same sequences
should give the model two independent measurements per unique site AND
encourage strand-invariant representations.

If RC augmentation helps (~0.34+ vs 0.32 random Gosai), strand-symmetry
prior is valuable. If it doesn't help, model already learns invariance or
training capacity isn't the limit.

25K unique sequences + their 25K reverse complements = 50K rows.
Each unique appears twice (forward + RC) with separate noisy labels.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_UNIQUE = 25_000

rng = np.random.default_rng(SEED)

COMP = str.maketrans("ACGT", "TGCA")

def rc(s):
    return s.translate(COMP)[::-1]

print("Loading Gosai...", flush=True)
seqs = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if any(b not in "ACGT" for b in c[11]): continue
        seqs.append(c[11])
print(f"Pool: {len(seqs)}", flush=True)

idx = rng.permutation(len(seqs))[:N_UNIQUE]
unique = [seqs[i] for i in idx]
rc_seqs = [rc(s) for s in unique]

# Verify all 200bp
for s in rc_seqs:
    assert len(s) == 200

final = unique + rc_seqs
rng.shuffle(final)
assert len(final) == 50_000

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {out}", flush=True)
