"""
E024: E022 replicate with seed=2.

Test seed stability of the new best library (test-chr + top SNR = 0.341).
If E024 ≈ 0.341: combo is robust.
If E024 drops significantly: chr+SNR has high variance and 0.341 is lucky.

Note: top-SNR ranking is deterministic, so the only variable is the
final shuffle order. Should give identical result if pipeline is order-invariant
to within seed variance.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 2
N_SEQS = 50_000
TEST_CHROMS = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai test chroms with quality scores...", flush=True)
rows = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in TEST_CHROMS: continue
        if any(b not in "ACGT" for b in c[11]): continue
        try:
            kfc, hfc, sfc = float(c[5]), float(c[6]), float(c[7])
            kse, hse, sse = float(c[8]), float(c[9]), float(c[10])
        except: continue
        snr = max(abs(kfc)/max(kse, 0.05),
                  abs(hfc)/max(hse, 0.05),
                  abs(sfc)/max(sse, 0.05))
        rows.append((c[11], snr))
print(f"Pool: {len(rows)}", flush=True)

rows.sort(key=lambda r: -r[1])
top = rows[:N_SEQS]
final = [r[0] for r in top]
rng.shuffle(final)
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
