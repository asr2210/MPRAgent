"""Experiment 008 — Malinois CELL-TYPE-DISCRIMINATING selection.

500k GC-50 random candidates scored with Malinois. Select 50k by composite
score = magnitude + 2*variance across (K562, HepG2, SKNSH) predictions.
Half from highest composite, weighted toward high cross-cell variance.

Hypothesis: cell-type-DISCRIMINATING sequences (those that activate strongly
in one cell type but not others) expose the model to the grammar of
cell-type specificity. This grammar should transfer to held-out cell types
because TFs are shared — only their cell-type expression patterns differ.

Generalization justification: a model that learns motif → cell-type-pattern
correspondences from our 3 cells should be able to predict ANY held-out
cell type's response to a sequence given the held-out cell's TF expression.
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


def gen_gc50_indices(rng, n):
    return rng.integers(0, 4, size=(n, SEQLEN), dtype=np.int8)


_UP_IDX, _DN_IDX = None, None
def get_flanks():
    global _UP_IDX, _DN_IDX
    if _UP_IDX is None:
        m = {"A": 0, "C": 1, "G": 2, "T": 3}
        _UP_IDX = np.array([m[c] for c in MPRA_UPSTREAM[-200:]], dtype=np.int8)
        _DN_IDX = np.array([m[c] for c in MPRA_DOWNSTREAM[:200]], dtype=np.int8)
    return _UP_IDX, _DN_IDX


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
    t0 = time.time()
    for i in range(0, N, BATCH):
        batch_idx = all_idx[i:i + BATCH]
        x = onehot_batch_from_idx(batch_idx, device)
        with torch.no_grad():
            y = model(x)
        preds[i:i + BATCH] = y.cpu().numpy()
        if (i // BATCH) % 50 == 0:
            dt = time.time() - t0
            done = i + batch_idx.shape[0]
            rate = done / max(dt, 1e-6)
            eta = (N - done) / max(rate, 1e-6)
            print(f"  scored {done:,}/{N:,} | {rate:.0f}/s | ETA {eta:.0f}s", flush=True)
    return preds


def idx_to_seqs(idx_arr):
    lookup = np.array(list("ACGT"))
    chars = lookup[idx_arr]
    return ["".join(row) for row in chars]


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print("Loading Malinois...")
    model = load_malinois(device)

    print(f"Generating {N_CAND:,} GC-50 random candidates...", flush=True)
    all_idx = gen_gc50_indices(rng, N_CAND)
    print(f"  Generated shape {all_idx.shape} in {time.time()-t0:.1f}s", flush=True)

    print(f"Scoring with Malinois (batch={BATCH})...", flush=True)
    preds = score_all(model, all_idx, device)

    # Compute per-sequence summary stats
    magnitude = preds.max(axis=1)
    mean_act = preds.mean(axis=1)
    variance = preds.std(axis=1)
    print(f"Magnitude: min={magnitude.min():.2f}, median={np.median(magnitude):.2f}, max={magnitude.max():.2f}")
    print(f"Variance:  min={variance.min():.3f}, median={np.median(variance):.3f}, max={variance.max():.3f}")

    # Composite score: rewards both magnitude AND cell-type discrimination.
    # Use rank-normalized to make them comparable.
    def rank_norm(x):
        order = np.argsort(x)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(len(x))
        return ranks / (len(x) - 1)
    score = rank_norm(magnitude) + 2.0 * rank_norm(variance)
    print(f"Composite score range: [{score.min():.3f}, {score.max():.3f}]")

    # Select top N_TARGET by composite
    sel = np.argpartition(-score, N_TARGET)[:N_TARGET]
    print(f"  Selected {len(sel)} sequences")

    selected_preds = preds[sel]
    print(f"Selected stats:")
    print(f"  K562  : mean {selected_preds[:,0].mean():.2f}, std {selected_preds[:,0].std():.2f}, range [{selected_preds[:,0].min():.2f}, {selected_preds[:,0].max():.2f}]")
    print(f"  HepG2 : mean {selected_preds[:,1].mean():.2f}, std {selected_preds[:,1].std():.2f}, range [{selected_preds[:,1].min():.2f}, {selected_preds[:,1].max():.2f}]")
    print(f"  SKNSH : mean {selected_preds[:,2].mean():.2f}, std {selected_preds[:,2].std():.2f}, range [{selected_preds[:,2].min():.2f}, {selected_preds[:,2].max():.2f}]")
    print(f"  per-seq cross-cell variance: mean {selected_preds.std(axis=1).mean():.3f}")

    # GC stats
    sel_idx = all_idx[sel]
    gc = ((sel_idx == 1) | (sel_idx == 2)).mean(axis=1)
    print(f"  GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    rng.shuffle(sel)
    selected_seqs = idx_to_seqs(all_idx[sel])
    np.savez(SCORES_OUT,
              preds_all=preds.astype(np.float16),
              selected_idx=sel.astype(np.int32),
              selected_preds=preds[sel].astype(np.float16))
    write_sequences(selected_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
