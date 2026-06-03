"""
Experiment 009: 50k Gosai sequences with largest |SKNSH_log2FC|.

Hypothesis: SK-N-SH r jumped 0.000 → 0.026 with random Gosai sample. If we filter
to extreme SK-N-SH effect sequences (top 50k by |SKNSH_log2FC|), the model gets
maximum SK-N-SH dynamic range to learn from. Should push SK-N-SH r higher.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "data", "evaluator_data", "41586_2024_8070_MOESM4_ESM.txt")
OUT = os.path.join(HERE, "sequences_0.txt")

records = []
with open(SRC) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    sknsh_idx = header.index("SKNSH_log2FC")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= seq_idx:
            continue
        s = p[seq_idx].upper()
        if len(s) != 200 or not all(c in "ACGT" for c in s):
            continue
        try:
            sk = float(p[sknsh_idx])
        except ValueError:
            continue
        records.append((abs(sk), s))
print(f"Valid records: {len(records):,}")

# Sort by |SKNSH| descending, take top 50k
records.sort(reverse=True)
chosen = [s for _, s in records[:50_000]]
print(f"Chosen: {len(chosen)}")
print(f"Top |SKNSH| value: {records[0][0]:.3f}, 50k-th: {records[49999][0]:.3f}")
print(f"50k-th |SKNSH|: {records[49999][0]:.3f}, median: {records[len(records)//2][0]:.3f}")

with open(OUT, "w") as f:
    for s in chosen:
        f.write(s + "\n")
