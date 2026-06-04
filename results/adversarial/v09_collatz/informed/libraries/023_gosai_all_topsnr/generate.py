"""
E023: Top 50K Gosai by SNR (no chr restriction).

E022 stacked chr-restriction + SNR for 0.341. Does SNR alone (no chr
filter) also beat random Gosai (0.323)?

If E023 > 0.323: SNR filter has its own boost (~+0.005 like in E022).
If E023 ≈ 0.323: SNR alone does nothing; only the chr+SNR combination works.

SNR = max across cells of |log2FC| / lfcSE.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_SEQS = 50_000

rng = np.random.default_rng(SEED)

print("Loading Gosai (all chroms) with quality...", flush=True)
rows = []
with open(SRC) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
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
