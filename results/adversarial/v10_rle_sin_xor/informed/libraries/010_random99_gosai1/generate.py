"""
Experiment 010: 99% random + 1% Gosai (500 biological sequences in random pool).

Hypothesis: K562 saturation requires mostly-random training data. A small Gosai
fraction may add enough biology to lift SK-N-SH without disrupting K562.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "data", "evaluator_data", "41586_2024_8070_MOESM4_ESM.txt")
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

# Load gosai 200bp sequences
gosai = []
with open(SRC) as f:
    h = f.readline().rstrip("\n").split("\t")
    sx = h.index("sequence")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= sx: continue
        s = p[sx].upper()
        if len(s) == 200 and all(c in "ACGT" for c in s):
            gosai.append(s)
print(f"Gosai pool: {len(gosai)}")

rng = np.random.default_rng(SEED)
gosai_pick = rng.choice(len(gosai), size=500, replace=False)
out = [gosai[i] for i in gosai_pick]

bases = np.array(list("ACGT"))
for _ in range(49_500):
    out.append("".join(bases[rng.integers(0, 4, 200)]))

rng.shuffle(out)
assert len(out) == 50_000
with open(OUT, "w") as f:
    for s in out:
        f.write(s + "\n")
print(f"Wrote {len(out)} (500 gosai + 49500 random)")
