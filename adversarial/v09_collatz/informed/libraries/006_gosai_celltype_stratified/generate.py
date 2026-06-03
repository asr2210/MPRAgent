"""
E006: Gosai cell-type-stratified — equal balance of sequences dominated
by each cell type's activity.

Rationale: my per-cell-type ceiling is SKNSH=0.62 (high), HepG2=0.20, K562=0.15.
The model seems undertrained for K562/HepG2-specific patterns. If I force the
library to contain equal numbers of K562-dominant, HepG2-dominant, and
SKNSH-dominant sequences, the model should learn each cell type more uniformly.

Selection:
- For each cell C in {K562, HepG2, SKNSH}: select sequences where C's log2FC
  is the maximum across all three cells AND C's log2FC >= 1.0 (clearly active).
- Take 16,666 from each cell-dominant pool, plus 50K - 50000 = 0 extras.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
PER_CELL = N_SEQS // 3  # 16666 each; +2 extras handled below

rng = np.random.default_rng(SEED)

print("Loading Gosai...", flush=True)
rows = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12: continue
        seq = c[11]
        if len(seq) != 200: continue
        if any(b not in "ACGT" for b in seq): continue
        try:
            rows.append((seq, float(c[5]), float(c[6]), float(c[7])))
        except: continue
print(f"Total clean 200bp: {len(rows)}", flush=True)

arr = np.array([r[1:] for r in rows])
seqs = [r[0] for r in rows]

k_dom = (arr[:,0] > arr[:,1]) & (arr[:,0] > arr[:,2]) & (arr[:,0] >= 1.0)
h_dom = (arr[:,1] > arr[:,0]) & (arr[:,1] > arr[:,2]) & (arr[:,1] >= 1.0)
s_dom = (arr[:,2] > arr[:,0]) & (arr[:,2] > arr[:,1]) & (arr[:,2] >= 1.0)

k_idx = np.where(k_dom)[0]
h_idx = np.where(h_dom)[0]
s_idx = np.where(s_dom)[0]
print(f"K-dom pool: {len(k_idx)} | H-dom pool: {len(h_idx)} | S-dom pool: {len(s_idx)}", flush=True)

assert min(len(k_idx), len(h_idx), len(s_idx)) >= PER_CELL + 1

k_sel = rng.choice(k_idx, PER_CELL, replace=False)
h_sel = rng.choice(h_idx, PER_CELL, replace=False)
s_sel = rng.choice(s_idx, PER_CELL, replace=False)

# 16666*3 = 49998; need 2 more
extra_pool = np.concatenate([
    np.setdiff1d(k_idx, k_sel),
    np.setdiff1d(h_idx, h_sel),
    np.setdiff1d(s_idx, s_sel),
])
extra = rng.choice(extra_pool, 2, replace=False)

all_sel = np.concatenate([k_sel, h_sel, s_sel, extra])
rng.shuffle(all_sel)
final = [seqs[i] for i in all_sel]
print(f"Selected {len(final)}", flush=True)
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {out}", flush=True)
