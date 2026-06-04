"""
E017: GC-content-matched random sequences.

Compute the GC distribution of Gosai sequences, then sample 50K
fully random sequences with per-sequence GC drawn from that distribution.

Tests: is the +0.06 lift of shuffled Gosai over uniform random just
about GC content, or do A-vs-T and C-vs-G ratios matter beyond GC?

- If E017 ≈ shuffled Gosai (0.305): GC alone is the signal.
- If E017 < shuffled Gosai (0.27-0.29): full mononucleotide marginals
  (incl. A/T and C/G asymmetry) carry extra info beyond GC.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
LEN = 200

rng = np.random.default_rng(SEED)

# Compute Gosai GC distribution
print("Computing Gosai GC distribution...", flush=True)
gcs = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if any(b not in "ACGT" for b in c[11]): continue
        s = c[11]
        gcs.append((s.count("G") + s.count("C")) / LEN)
gcs = np.array(gcs)
print(f"Gosai N={len(gcs)} GC mean={gcs.mean():.3f} std={gcs.std():.3f} "
      f"min={gcs.min():.3f} max={gcs.max():.3f}", flush=True)

# Sample 50K GC values from Gosai distribution
gc_samples = rng.choice(gcs, N_SEQS, replace=True)

# For each, generate a random sequence with that GC content
final = []
for gc in gc_samples:
    # Random binary: each position is GC with prob gc, else AT
    is_gc = rng.random(LEN) < gc
    # Within GC pick G or C uniformly; within AT pick A or T uniformly
    g_or_c = rng.integers(0, 2, LEN)  # 0=G, 1=C
    a_or_t = rng.integers(0, 2, LEN)  # 0=A, 1=T
    seq = []
    for i in range(LEN):
        if is_gc[i]:
            seq.append("G" if g_or_c[i] == 0 else "C")
        else:
            seq.append("A" if a_or_t[i] == 0 else "T")
    final.append("".join(seq))

assert len(final) == N_SEQS
for s in final[:5]:
    assert len(s) == LEN and all(b in "ACGT" for b in s)

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
