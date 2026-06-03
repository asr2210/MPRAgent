"""
Experiment 008: real MPRA sequences from Gosai et al. 2024 (Nature 634:1211).

This dataset has 798k sequences (~763k at 200bp) measured in K562/HepG2/SK-N-SH.
It is almost certainly the source for the `mpra_oracle` baseline in instructions.md
Table 1, which reported eval_01=0.6643. If my pipeline gives similar, Table 1
applies and biology IS the key. If my pipeline gives ~0.52, Table 1 is irrelevant.

Either way, this is the highest-value diagnostic I can run.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "data", "evaluator_data", "41586_2024_8070_MOESM4_ESM.txt")
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

# Read all 200bp sequences
seqs = []
with open(SRC) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= seq_idx:
            continue
        s = p[seq_idx].upper()
        if len(s) == 200 and all(c in "ACGT" for c in s):
            seqs.append(s)
print(f"Loaded {len(seqs):,} valid 200bp sequences")

rng = np.random.default_rng(SEED)
idxs = rng.choice(len(seqs), size=50_000, replace=False)
with open(OUT, "w") as f:
    for i in idxs:
        f.write(seqs[i] + "\n")
print(f"Wrote 50,000 sequences")
