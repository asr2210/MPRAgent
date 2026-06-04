"""
Experiment 015 — Stricter quality (max-SE filter).

Previous filters used mean(lfcSE) < threshold. This still admits sequences with
ONE noisy cell type. Now require max(k_se, h_se, s_se) < 0.3 — all three cells
must have clean measurements.

Random sample (per exp 014 finding that stratification doesn't add lift).
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200
MAX_SE = 0.3

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
        if max(k_se, h_se, s_se) >= MAX_SE:
            continue
        seqs_pool.append(seq)

print(f"{len(seqs_pool)} sequences after max(lfcSE)<{MAX_SE}")
random.shuffle(seqs_pool)
seqs = seqs_pool[:TARGET_N]
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
