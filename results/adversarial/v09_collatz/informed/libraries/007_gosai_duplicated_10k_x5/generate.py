"""
E007: Diagnostic — 10K unique sequences each repeated 5x = 50K rows.

Hypothesis: my pipeline's per-cell-type ceiling (K562 0.15, HepG2 0.20)
might reflect HIGH measurement noise. If each row gets one noisy measurement,
duplicating a sequence 5x gives the training process 5 independent noise
samples to average over → cleaner per-sequence label → better model.

If E007 substantially beats E004 (Gosai random = 0.32), the bottleneck IS
label noise on individual measurements. If E007 ≈ E004, then the bottleneck
is elsewhere (model capacity, eval distribution mismatch, etc.).

Take 10K highest-|log2FC|-magnitude Gosai sequences (most informative
labels), then duplicate each 5x to fill 50K rows.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_UNIQUE = 10_000
N_COPIES = 5
N_SEQS = N_UNIQUE * N_COPIES

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
print(f"Total: {len(rows)}", flush=True)

arr = np.array([r[1:] for r in rows])
seqs = [r[0] for r in rows]

# Score by max |log2FC| across cells
score = np.max(np.abs(arr), axis=1)
top = np.argsort(-score)[:N_UNIQUE]
unique_seqs = [seqs[i] for i in top]
print(f"Top {N_UNIQUE} max|log2FC|: range {np.min(score[top]):.2f} - {np.max(score[top]):.2f}", flush=True)

# Replicate each 5 times, then shuffle
final = unique_seqs * N_COPIES  # 10K * 5 = 50K
rng.shuffle(final)
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {out}", flush=True)
