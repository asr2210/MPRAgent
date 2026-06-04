"""Experiment 017 — Multi-seed oracle pool (5 × 10k from 5 independent seeds).

Tests whether averaging the oracle selection across multiple seeds reduces
the noise we see between 007 (mean_r=0.3861) and 014 (mean_r=0.3827, same
recipe, different seed).

Design: 5 independent oracle runs, each with 100k GC-50 random pool + Malinois
scoring + top-10k selection. Concatenate the 50k for the final library.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np
import torch

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
sys.path.insert(0, "/data/users/arao/.private/mpra_exp/boda2")

from utils.seqlib import write_sequences, SEQLEN
from boda.model.basset import BassetBranched
from boda.common.constants import MPRA_UPSTREAM, MPRA_DOWNSTREAM

OUT = os.path.join(THIS, "sequences_0.txt")
N_SEEDS = 5
POOL_PER_SEED = 100_000
TOP_PER_SEED = 10_000
N_TARGET = N_SEEDS * TOP_PER_SEED
BATCH = 1024
ART_DIR = "/data/users/arao/.private/mpra_exp/data/malinois_artifacts/artifacts"


def load_malinois(device):
    ckpt = torch.load(os.path.join(ART_DIR, "torch_checkpoint.pt"),
                       map_location="cpu", weights_only=False)
    hp = vars(ckpt["model_hparams"])
    model = BassetBranched(**hp)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model.to(device)


_UP, _DN = None, None
def get_flanks():
    global _UP, _DN
    if _UP is None:
        m = {"A": 0, "C": 1, "G": 2, "T": 3}
        _UP = np.array([m[c] for c in MPRA_UPSTREAM[-200:]], dtype=np.int8)
        _DN = np.array([m[c] for c in MPRA_DOWNSTREAM[:200]], dtype=np.int8)
    return _UP, _DN


def onehot_batch_from_idx(core_idx, device):
    up, dn = get_flanks()
    B = core_idx.shape[0]
    full = np.empty((B, 600), dtype=np.int64)
    full[:, :200] = up[None, :]
    full[:, 200:400] = core_idx
    full[:, 400:600] = dn[None, :]
    x = torch.zeros(B, 4, 600, device=device)
    x.scatter_(1, torch.from_numpy(full).to(device).unsqueeze(1), 1.0)
    return x


def score_all(model, all_idx, device):
    N = all_idx.shape[0]
    preds = np.zeros((N, 3), dtype=np.float32)
    for i in range(0, N, BATCH):
        b = all_idx[i:i + BATCH]
        x = onehot_batch_from_idx(b, device)
        with torch.no_grad():
            y = model(x)
        preds[i:i + BATCH] = y.cpu().numpy()
    return preds


def idx_to_seqs(idx_arr):
    lookup = np.array(list("ACGT"))
    return ["".join(row) for row in lookup[idx_arr]]


def main():
    t0 = time.time()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    print("Loading Malinois...")
    model = load_malinois(device)

    all_selected = []
    for seed in range(1, N_SEEDS + 1):
        t_s = time.time()
        rng = np.random.default_rng(seed * 1000 + 7)
        print(f"\n[seed={seed}] Generating {POOL_PER_SEED:,} pool...")
        pool = rng.integers(0, 4, size=(POOL_PER_SEED, SEQLEN), dtype=np.int8)
        preds = score_all(model, pool, device)
        score = preds.max(axis=1)
        sel = np.argpartition(-score, TOP_PER_SEED)[:TOP_PER_SEED]
        selected = pool[sel]
        all_selected.append(selected)
        print(f"  selected top-{TOP_PER_SEED} | "
              f"score min={score[sel].min():.2f}, max={score[sel].max():.2f} | "
              f"{time.time()-t_s:.1f}s")

    combined = np.concatenate(all_selected)
    assert len(combined) == N_TARGET
    # Shuffle
    rng = np.random.default_rng(99)
    perm = rng.permutation(N_TARGET)
    combined = combined[perm]

    gc = ((combined == 1) | (combined == 2)).mean(axis=1)
    print(f"\nFinal lib GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    seqs = idx_to_seqs(combined)
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
