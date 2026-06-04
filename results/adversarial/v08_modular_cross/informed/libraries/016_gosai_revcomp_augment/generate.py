"""
Experiment 016 — Reverse-complement augmentation.

Take 25K Gosai (lfcSE<0.5, random) and pair each with its reverse complement.
Final library = 25K real + 25K revcomps (50K total).

Test:
- If Malinois oracle is RC-invariant → revcomps get same labels → 2x training
  signal per real measurement. Model overfits less.
- If oracle is RC-aware → revcomps get different labels → model learns
  RC-invariance explicitly.

Either way, this asks: is the bottleneck the model not seeing enough
sequence-activity examples, or is it the oracle's label quality?
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200
N_REAL = 25_000

COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(COMP)[::-1]

seqs_pool = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    cols = {c: header.index(c) for c in
            ["K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k_se = float(parts[cols["K562_lfcSE"]])
            h_se = float(parts[cols["HepG2_lfcSE"]])
            s_se = float(parts[cols["SKNSH_lfcSE"]])
        except (ValueError, IndexError):
            continue
        if (k_se + h_se + s_se) / 3 >= 0.5:
            continue
        seqs_pool.append(seq)

print(f"{len(seqs_pool)} sequences after lfcSE<0.5")
random.shuffle(seqs_pool)
real = seqs_pool[:N_REAL]
rc = [revcomp(s) for s in real]
combined = real + rc
random.shuffle(combined)
assert len(combined) == TARGET_N

# Sanity: ensure all are 200bp ACGT
for s in combined:
    assert len(s) == SEQ_LEN and set(s).issubset(set("ACGT"))

with open(OUT, "w") as f:
    for s in combined:
        f.write(s + "\n")
print(f"wrote {len(combined)} (25K real + 25K revcomp) to {OUT}")
