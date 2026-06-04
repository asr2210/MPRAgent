"""Experiment 013 — Malinois TOP from 4M candidate pool (8× larger than 007).

Tests whether the oracle benefit scales with candidate pool size, or whether
007's top-50k from 500k already saturates the oracle's discrimination ability.

Prediction (theory v6):
- If 013 > 007 by >0.003: oracle scales → more pool = better.
- If 013 ≈ 007: oracle saturates around 500k.
- If 013 < 007 (like 011 hyperactive): too-extreme sequences hurt downstream.

Design: 4M GC=50 random ACGT, score with Malinois, top 50k by
max(K562, HepG2, SKNSH). Stream in chunks to manage memory.
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
N_CAND = 4_000_000
CHUNK = 200_000          # generate + score in chunks of 200k
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


def score_chunk(model, idx_chunk, device):
    N = idx_chunk.shape[0]
    preds = np.zeros((N, 3), dtype=np.float32)
    for i in range(0, N, BATCH):
        b = idx_chunk[i:i + BATCH]
        x = onehot_batch_from_idx(b, device)
        with torch.no_grad():
            y = model(x)
        preds[i:i + BATCH] = y.cpu().numpy()
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

    # We need to remember the top N_TARGET by max-cell-score globally,
    # but cannot hold 4M×200 idx in memory comfortably (~800MB). Stream chunks.
    # Keep a running heap of (score, idx_array).
    # Simpler: keep top-K idx + scores; for each chunk, score and merge.
    top_scores = np.full(N_TARGET, -np.inf, dtype=np.float32)
    top_idx = np.zeros((N_TARGET, SEQLEN), dtype=np.int8)
    top_preds = np.zeros((N_TARGET, 3), dtype=np.float32)

    n_total = 0
    n_chunks = (N_CAND + CHUNK - 1) // CHUNK
    for c in range(n_chunks):
        bs = min(CHUNK, N_CAND - n_total)
        t_c = time.time()
        idx_chunk = rng.integers(0, 4, size=(bs, SEQLEN), dtype=np.int8)
        preds = score_chunk(model, idx_chunk, device)
        chunk_score = preds.max(axis=1)

        # Merge with running top-K: concatenate and re-select top
        merged_scores = np.concatenate([top_scores, chunk_score])
        merged_idx = np.concatenate([top_idx, idx_chunk])
        merged_preds = np.concatenate([top_preds, preds])
        sel = np.argpartition(-merged_scores, N_TARGET)[:N_TARGET]
        top_scores = merged_scores[sel]
        top_idx = merged_idx[sel]
        top_preds = merged_preds[sel]

        n_total += bs
        dt_c = time.time() - t_c
        dt_total = time.time() - t0
        rate = n_total / max(dt_total, 1e-6)
        eta = (N_CAND - n_total) / max(rate, 1e-6)
        thresh = top_scores.min()
        print(f"  chunk {c+1}/{n_chunks}: total {n_total:,}/{N_CAND:,} | "
              f"{dt_c:.0f}s | overall {rate:.0f}/s | ETA {eta:.0f}s | "
              f"top-K min={thresh:.3f}, max={top_scores.max():.3f}", flush=True)

    # Final report
    print(f"\nSelection complete. Top {N_TARGET:,} chosen from {N_CAND:,}.")
    print(f"Selected score: min={top_scores.min():.3f}, "
          f"mean={top_scores.mean():.3f}, max={top_scores.max():.3f}")
    print(f"Selected pred K562: mean {top_preds[:,0].mean():.3f}")
    print(f"Selected pred HepG2: mean {top_preds[:,1].mean():.3f}")
    print(f"Selected pred SKNSH: mean {top_preds[:,2].mean():.3f}")
    # GC
    gc = ((top_idx == 1) | (top_idx == 2)).mean(axis=1)
    print(f"Selected GC: mean={gc.mean():.3f}, std={gc.std():.3f}")

    # Shuffle and write
    order = rng.permutation(N_TARGET)
    top_idx = top_idx[order]
    top_scores = top_scores[order]
    top_preds = top_preds[order]
    seqs = idx_to_seqs(top_idx)
    np.savez(SCORES_OUT,
              top_scores=top_scores.astype(np.float16),
              top_preds=top_preds.astype(np.float16))
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
