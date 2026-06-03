"""
E016: Gosai shuffled — k-mer scrambled Gosai sequences.

Shuffles bases within each Gosai sequence (preserves single-base
composition, destroys motifs/grammar).

Hypothesis test: does Gosai beat random because of motif content
(real grammar) or just because of composition statistics (GC%, etc.)?

If E016 ≈ random (E001 = 0.244): the boost from Gosai (+0.08) is mostly
from real motif grammar.
If E016 ≈ Gosai (E004 = 0.323): the boost is just composition statistics —
the model learns from k-mer marginals, not from sequence-level grammar.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000

rng = np.random.default_rng(SEED)

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

idx = rng.permutation(len(seqs))[:N_SEQS]
selected = [seqs[i] for i in idx]

# Shuffle each sequence in place (preserves base composition)
final = []
for s in selected:
    arr = np.array(list(s))
    rng.shuffle(arr)
    final.append("".join(arr))
assert len(final) == N_SEQS
for s in final[:3]:
    assert len(s) == 200 and all(b in "ACGT" for b in s)

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
