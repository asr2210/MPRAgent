"""
Experiment 025 — GTEX lfcSE<0.7 + |mean activity| > 0.3.

Adds a deterministic effect-magnitude filter to remove near-null variants
from the GTEX-loose pool. Hypothesis: near-null variants add noise to the
training (label ≈ 0 with random sequence) — removing them reduces sampling
variance while keeping high-information sequences.

Sub-source: GTEX
Quality: lfcSE<0.7
Effect: |mean activity| > 0.3
Sample: random within filtered pool
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200
SE_CUT = 0.7
EFF_CUT = 0.3

seqs_pool = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    proj_idx = header.index("data_project")
    cols = {c: header.index(c) for c in
            ["K562_log2FC", "HepG2_log2FC", "SKNSH_log2FC",
             "K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[proj_idx] != "GTEX":
            continue
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k562 = float(parts[cols["K562_log2FC"]])
            hepg2 = float(parts[cols["HepG2_log2FC"]])
            sknsh = float(parts[cols["SKNSH_log2FC"]])
            k_se = float(parts[cols["K562_lfcSE"]])
            h_se = float(parts[cols["HepG2_lfcSE"]])
            s_se = float(parts[cols["SKNSH_lfcSE"]])
        except (ValueError, IndexError):
            continue
        if (k_se + h_se + s_se) / 3 >= SE_CUT:
            continue
        mean_act = (k562 + hepg2 + sknsh) / 3
        if abs(mean_act) <= EFF_CUT:
            continue
        seqs_pool.append(seq)

print(f"{len(seqs_pool)} GTEX seqs after lfcSE<{SE_CUT} AND |mean|>{EFF_CUT}")
random.shuffle(seqs_pool)
seqs = seqs_pool[:TARGET_N]
assert len(seqs) == TARGET_N
with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
