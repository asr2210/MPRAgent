"""
Experiment 014 — Random Gosai with quality filter (no stratification).

Isolates the contribution of stratification vs quality.
- exp 003 (random Gosai, no filter): eval_01 = 0.0000
- exp 004 (quality + stratified): eval_01 = 0.0181
Two changes between them. This experiment runs random sampling AFTER quality
filter — same pool as exp 004 but no stratification. Compares:
  - if eval_01 ~ 0.018 → quality alone gives the lift; stratification doesn't help
  - if eval_01 << 0.018 → stratification is essential
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
        mean_se = (k_se + h_se + s_se) / 3
        if mean_se >= 0.5:
            continue
        seqs_pool.append(seq)

print(f"{len(seqs_pool)} sequences after lfcSE<0.5")
random.shuffle(seqs_pool)
seqs = seqs_pool[:TARGET_N]
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
