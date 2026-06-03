"""
Experiment 011: composition-mode-stratified Dirichlet.

Generate 500,000 Dirichlet(0.5) compositions. Categorize each by
(dominant base, second-dominant base) — 4*3 = 12 modes. Sample equal
numbers from each mode (50000/12 ≈ 4166 per mode) to ensure even
representation across composition geometries.

Tests whether Dirichlet(0.5) at N=50k undersamples some composition
modes. Hypothesis: dirichlet's natural draw may over-represent the
"AT-dominant" or "GC-dominant" axis depending on luck; stratifying
forces uniform mode coverage.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
ALPHA = (0.5, 0.5, 0.5, 0.5)
CANDIDATES = 500_000
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)

    cand = rng.dirichlet(ALPHA, size=CANDIDATES)               # (M, 4)
    sort = np.argsort(-cand, axis=1)                            # (M, 4)
    dom = sort[:, 0]                                            # (M,)
    sec = sort[:, 1]                                            # (M,)
    mode = dom * 4 + sec                                        # 0..15 (skips dom==sec)
    uniq, inv = np.unique(mode, return_inverse=True)
    n_modes = len(uniq)
    per_mode = N // n_modes
    extra = N - per_mode * n_modes

    selected_idx = []
    for k, m in enumerate(uniq):
        idx_in_mode = np.where(inv == k)[0]
        take = per_mode + (1 if k < extra else 0)
        if len(idx_in_mode) < take:
            # shouldn't happen but be safe
            take = len(idx_in_mode)
        chosen = rng.choice(idx_in_mode, size=take, replace=False)
        selected_idx.append(chosen)
    selected_idx = np.concatenate(selected_idx)
    rng.shuffle(selected_idx)
    selected_idx = selected_idx[:N]
    probs = cand[selected_idx]

    print(f"modes used: {n_modes}, per-mode counts (first 5): "
          f"{[len(rng.choice(np.where(inv==k)[0], size=per_mode, replace=False)) for k in range(min(5, n_modes))]}")
    print(f"selected: {len(probs)}")

    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} mode-stratified dirichlet seqs to {OUT}")


if __name__ == "__main__":
    main()
