"""
E027: Top 50K Gosai by SNR, restricted to Malinois held-out chroms
(chr 7, 9, 13, 19, 21, X) — broader hold-out set than E022.

E020 (random chr 7/9/13/19/21/X) was 0.333 — diluted by chr19.
E022 (top SNR chr 7/9/13/21/X) was 0.341.

Now combine: SNR ranking + chr19 inclusion. If chr19 has high-SNR
sequences that help, score should rise. If chr19 adds dilution that
SNR can't compensate, score stays ~0.341 or below.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000
CHROMS = {"7", "9", "13", "19", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai held-out chroms with SNR...", flush=True)
rows = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if c[1] not in CHROMS: continue
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
print(f"Top SNR range: {top[0][1]:.2f} to {top[-1][1]:.2f}", flush=True)

final = [r[0] for r in top]
rng.shuffle(final)
assert len(final) == N_SEQS

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
