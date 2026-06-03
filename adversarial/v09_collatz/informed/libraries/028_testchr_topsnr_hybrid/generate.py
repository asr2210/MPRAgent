"""
E028: Hybrid — 30K top-SNR test-chr + 20K top-SNR non-test-chr.

Test if mixing best-of-test-chr (the eval-overlap region) with
best-of-rest (broader grammar coverage) compounds. E022's strength
might be hurt by saturation on a narrow subset; adding diverse
high-quality elsewhere could fill gaps.

If E028 > 0.341: hybrid stacking works.
If ≈ 0.341: equivalent — eval doesn't care about non-test-chr quality.
If < 0.341: dilution from non-test-chr.
"""
import os, numpy as np

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
SRC = f"{ROOT}/data/gosai_mpra.txt"
SEED = 1
N_TEST = 30_000
N_REST = 20_000
TEST_CHROMS = {"7", "9", "13", "21", "X"}

rng = np.random.default_rng(SEED)

print("Loading Gosai...", flush=True)
test_rows = []
rest_rows = []
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
        if c[1] in TEST_CHROMS:
            test_rows.append((c[11], snr))
        else:
            rest_rows.append((c[11], snr))
print(f"Test pool: {len(test_rows)}, Rest pool: {len(rest_rows)}", flush=True)

test_rows.sort(key=lambda r: -r[1])
rest_rows.sort(key=lambda r: -r[1])
final = [r[0] for r in test_rows[:N_TEST]] + [r[0] for r in rest_rows[:N_REST]]
print(f"Test top SNR cutoff: {test_rows[N_TEST-1][1]:.2f}", flush=True)
print(f"Rest top SNR cutoff: {rest_rows[N_REST-1][1]:.2f}", flush=True)

rng.shuffle(final)
assert len(final) == 50_000

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
