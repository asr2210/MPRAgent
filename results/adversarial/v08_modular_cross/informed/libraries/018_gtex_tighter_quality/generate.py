"""
Experiment 018 — GTEX with tighter quality (lfcSE<0.3) random.

Tightens quality on GTEX-only. exp 011 was lfcSE<0.3 stratified (0.0148).
This is lfcSE<0.3 random — parallel to 017 (lfcSE<0.5 random = 0.0190).
Tests whether more aggressive quality on GTEX beats the breakthrough.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200

seqs_pool = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    proj_idx = header.index("data_project")
    cols = {c: header.index(c) for c in
            ["K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[proj_idx] != "GTEX":
            continue
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k_se = float(parts[cols["K562_lfcSE"]])
            h_se = float(parts[cols["HepG2_lfcSE"]])
            s_se = float(parts[cols["SKNSH_lfcSE"]])
        except (ValueError, IndexError):
            continue
        if (k_se + h_se + s_se) / 3 >= 0.3:
            continue
        seqs_pool.append(seq)

print(f"{len(seqs_pool)} GTEX sequences after lfcSE<0.3")
random.shuffle(seqs_pool)
seqs = seqs_pool[:TARGET_N]
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
