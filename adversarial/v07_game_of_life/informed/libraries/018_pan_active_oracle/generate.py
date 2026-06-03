"""Experiment 018 — Pan-active Malinois oracle (top by MIN cross-cell).

Hypothesis: sequences that activate ALL THREE training cell types should be
more likely to activate held-out cell types in the evaluator. Standard 007
picks by max-cell (cell-specific extremes); 018 picks by min-cell
(pan-active, transferable).

Design: 500k GC-50 random, score with Malinois, top 50k by
MIN(K562, HepG2, SKNSH).
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
SCORES_OUT = os.path.join(THIS, "scores.npz")
SEED = 1
N_CAND = 500_000
N_TARGET = 50_000
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
    rng = np.random.default_rng(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    print("Loading Malinois...")
    model = load_malinois(device)

    print(f"Generating {N_CAND:,} GC-50 random candidates...")
    pool = rng.integers(0, 4, size=(N_CAND, SEQLEN), dtype=np.int8)
    print(f"Scoring with Malinois...")
    preds = score_all(model, pool, device)
    print(f"All preds: K562 [{preds[:,0].min():.2f}, {preds[:,0].max():.2f}] mean {preds[:,0].mean():.3f}")
    print(f"           HepG2 [{preds[:,1].min():.2f}, {preds[:,1].max():.2f}] mean {preds[:,1].mean():.3f}")
    print(f"           SKNSH [{preds[:,2].min():.2f}, {preds[:,2].max():.2f}] mean {preds[:,2].mean():.3f}")

    # Pan-active: MIN across cells
    score = preds.min(axis=1)
    print(f"min-cell score: min={score.min():.3f}, max={score.max():.3f}, "
          f"median={np.median(score):.3f}")
    sel = np.argpartition(-score, N_TARGET)[:N_TARGET]
    print(f"Selected min-score: min={score[sel].min():.3f}, "
          f"max={score[sel].max():.3f}, mean={score[sel].mean():.3f}")
    selected_preds = preds[sel]
    print(f"Selected preds: K562 mean {selected_preds[:,0].mean():.3f}, "
          f"HepG2 {selected_preds[:,1].mean():.3f}, SKNSH {selected_preds[:,2].mean():.3f}")
    gc = ((pool[sel] == 1) | (pool[sel] == 2)).mean(axis=1)
    print(f"Selected GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    rng.shuffle(sel)
    seqs = idx_to_seqs(pool[sel])
    np.savez(SCORES_OUT,
              selected_preds=preds[sel].astype(np.float16))
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
