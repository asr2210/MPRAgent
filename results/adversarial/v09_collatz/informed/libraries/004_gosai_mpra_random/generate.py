"""
E004: 50k sequences from the Gosai et al. 2024 Malinois MPRA dataset
sampled uniformly at random.

The Gosai dataset (~798K sequences) was measured in K562, HepG2, SK-N-SH
— exactly the three cell types prepare.py uses. These sequences have
proven measurable activity. Hypothesis: this is the "right kind" of
training distribution for real-label learning, matching baseline
`mpra_real` (0.6026 on eval_01).

Filter to exactly 200bp sequences (763,684 available). Sample 50k random.
"""
import os, gzip, sys
import numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
LEN = 200

rng = np.random.default_rng(SEED)

# Load all 200bp sequences
print("Loading Gosai dataset...", flush=True)
seqs = []
with open(SRC) as f:
    header = f.readline()
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 12:
            continue
        seq = cols[11]
        if len(seq) != LEN:
            continue
        # Skip if any non-ACGT
        if any(b not in "ACGT" for b in seq):
            continue
        seqs.append(seq)
print(f"Loaded {len(seqs)} 200bp clean sequences", flush=True)

# Sample 50K uniformly without replacement
idx = rng.permutation(len(seqs))[:N_SEQS]
selected = [seqs[i] for i in idx]
print(f"Selected {len(selected)}", flush=True)
assert len(selected) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(selected))
    f.write("\n")
print(f"Wrote {out}", flush=True)
