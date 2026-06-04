"""
E005: Gosai sequences filtered for HIGH activity AND LOW noise.

Selection criteria (within Gosai 200bp sequences):
1. max(|K562_log2FC|, |HepG2_log2FC|, |SKNSH_log2FC|) >= 1.0 (clearly active in ≥1 cell)
2. lfcSE for all three cells < median (confident labels)

This gives a "loud and clean" subset. Hypothesis: under real-label
training, sequences with high-signal/low-noise labels train the model
more effectively than randomly sampled MPRA sequences.

If E005 substantially beats E004 (Gosai random = 0.32), the bottleneck
is sequence quality. If similar, the bottleneck is something else.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000

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
            k = float(c[5]); h = float(c[6]); s = float(c[7])
            kse = float(c[8]); hse = float(c[9]); sse = float(c[10])
        except: continue
        rows.append((seq, k, h, s, kse, hse, sse))
print(f"Total: {len(rows)}", flush=True)

arr = np.array([r[1:] for r in rows])
seqs = [r[0] for r in rows]

# Compute median lfcSE per cell
med_kse = np.median(arr[:,3]); med_hse = np.median(arr[:,4]); med_sse = np.median(arr[:,5])
print(f"Median lfcSE K562={med_kse:.3f} HepG2={med_hse:.3f} SKNSH={med_sse:.3f}", flush=True)

max_abs = np.maximum.reduce([np.abs(arr[:,0]), np.abs(arr[:,1]), np.abs(arr[:,2])])

mask = (
    (max_abs >= 1.0) &
    (arr[:,3] < med_kse) &
    (arr[:,4] < med_hse) &
    (arr[:,5] < med_sse)
)
print(f"Mask passes: {mask.sum()}", flush=True)

valid_idx = np.where(mask)[0]
if len(valid_idx) < N_SEQS:
    # Relax: only require |log2FC| >= 1
    print("Not enough — relaxing SE filter to any 2 of 3 < median")
    se_ok = ((arr[:,3] < med_kse).astype(int) + (arr[:,4] < med_hse).astype(int) + (arr[:,5] < med_sse).astype(int)) >= 2
    mask = (max_abs >= 1.0) & se_ok
    valid_idx = np.where(mask)[0]
    print(f"Relaxed mask passes: {mask.sum()}", flush=True)

assert len(valid_idx) >= N_SEQS, f"Only {len(valid_idx)}"

# Sample N_SEQS uniformly from valid
sel_idx = rng.permutation(len(valid_idx))[:N_SEQS]
final = [seqs[valid_idx[i]] for i in sel_idx]

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} sequences to {out}", flush=True)
