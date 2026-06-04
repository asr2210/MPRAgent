"""
Experiment 003 — Gosai MPRA dataset, train chromosomes.

Sample 50,000 sequences from the Gosai et al. 2024 MPRA dataset (Table_S2),
restricted to chromosomes typically used as training in Gosai's model splits.
Avoid chromosomes commonly held out: 7, 9, 13, 21, X, Y.

Hypothesis: v08's MPRA oracle and eval distribution come from the Gosai
K562/HepG2/SK-N-SH MPRA dataset (the only public MPRA dataset using exactly
these three cell types at 200bp). If true, training on Gosai sequences (from
non-held-out chromosomes) should produce a model that generalizes to eval.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200

# Chromosomes commonly held out in Gosai-style splits (Boda2/Malinois validation
# + test = {7, 9, 13, 21, X}). Also drop Y (only 19 sequences) for cleanliness.
HOLD_OUT = {"7", "9", "13", "21", "X", "Y"}

candidates = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    chr_idx = header.index("chr")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom = parts[chr_idx]
        if chrom in HOLD_OUT:
            continue
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN:
            continue
        if not set(seq).issubset(set("ACGT")):
            continue
        candidates.append(seq)

print(f"{len(candidates)} valid candidate sequences after filtering")
assert len(candidates) >= TARGET_N

random.shuffle(candidates)
seqs = candidates[:TARGET_N]

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
