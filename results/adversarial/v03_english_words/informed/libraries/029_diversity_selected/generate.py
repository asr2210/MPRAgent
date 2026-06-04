"""
Experiment 029: Oversample + 6-mer diversity selection.

Strategy:
  1. Generate 150k candidates from a MIX of best-known designs:
     - 50k via exp 020 design (dense overlap, uniform bg)
     - 50k via exp 022 design (dense10, uniform bg — SK-N-SH-friendly)
     - 50k via exp 026 design (dinuc bg, dense 8 — SK-N-SH-friendliest)
  2. Score each candidate's NOVELTY by counting how many of its 6-mers
     are still rare in the accumulating selected set.
  3. Greedily pick 50k highest-novelty seqs.

Bypasses per-cell-type opt-conflict by letting selection criterion
(diversity) determine library composition. If the eval rewards
training-data diversity, this should improve over single-design libs.
"""

import re
from pathlib import Path
import numpy as np
from collections import Counter

ROOT = Path(__file__).resolve().parents[2]
JASPAR = ROOT / "data" / "JASPAR2024_vertebrates.txt"
OUT = Path(__file__).parent / "sequences_0.txt"

N_FINAL = 50_000
N_PER_POOL = 50_000
LEN = 200
SEED = 0
K = 6  # 6-mer for diversity scoring
BASES = np.array(list("ACGT"))

DINUC_PROBS = np.array([
    [0.327, 0.176, 0.240, 0.257],
    [0.354, 0.252, 0.048, 0.346],
    [0.282, 0.193, 0.247, 0.278],
    [0.213, 0.211, 0.252, 0.324],
])
BASE_FREQS = np.array([0.295, 0.205, 0.205, 0.295])

TARGET_TOKENS = [
    "GATA1", "GATA2", "GATA3", "KLF1", "KLF4", "KLF15", "NFE2", "MAF::NFE2",
    "MYB", "TAL1", "RUNX1", "SPI1",
    "HNF4A", "HNF4G", "HNF1A", "HNF1B", "FOXA1", "FOXA2", "FOXA3",
    "CEBPA", "CEBPB", "CEBPD", "CEBPG", "ONECUT1", "ONECUT2", "ONECUT3",
    "RXRA", "NR1H4", "PPARA",
    "NEUROD1", "NEUROG1", "NEUROG2", "ASCL1", "ASCL2", "OLIG1", "OLIG2", "OLIG3",
    "SOX2", "SOX10", "SOX21", "POU3F1", "POU3F2", "POU3F3", "POU3F4",
    "ISL2", "MEF2C", "MYCN", "PHOX2A", "PHOX2B", "REST", "RFX1", "RFX3", "RFX5",
    "PAX3", "PAX6",
    "SP1", "SP2", "SP3", "CTCF", "TBP", "NFYA", "NFYB", "NFYC",
    "FOS", "JUN", "FOSL", "JUND", "BACH",
    "ELK1", "ELK4", "ETV", "GABPA",
    "E2F1", "E2F4", "E2F6",
    "MYC", "MAX", "MAZ",
    "ATF", "CREB", "USF1", "USF2",
    "YY1", "EGR1", "NRF1",
]


def parse_jaspar_filtered(path: Path) -> list[tuple[str, np.ndarray]]:
    text = path.read_text()
    pfms = []
    blocks = re.split(r"\n>", "\n" + text.lstrip())
    targets_upper = [t.upper() for t in TARGET_TOKENS]
    for blk in blocks:
        blk = blk.strip()
        if not blk:
            continue
        lines = blk.splitlines()
        m = re.match(r"(\S+)\s+(\S+)", lines[0].strip())
        if not m:
            continue
        name = m.group(2).upper()
        if not any(t in name for t in targets_upper):
            continue
        rows = {}
        for ln in lines[1:5]:
            mb = re.match(r"\s*([ACGT])\s*\[\s*(.*?)\s*\]\s*$", ln)
            if not mb:
                break
            rows[mb.group(1)] = [float(x) for x in mb.group(2).split()]
        if set(rows.keys()) != set("ACGT"):
            continue
        L = len(rows["A"])
        if L < 4 or L > 30:
            continue
        mat = np.stack([rows["A"], rows["C"], rows["G"], rows["T"]], axis=0) + 0.1
        mat = mat / mat.sum(axis=0, keepdims=True)
        pfms.append((name, mat))
    return pfms


def sample_from_pfm(pfm: np.ndarray, rng: np.random.Generator) -> str:
    cols = [rng.choice(4, p=pfm[:, j]) for j in range(pfm.shape[1])]
    return "".join(BASES[c] for c in cols)


def random_bg(rng: np.random.Generator) -> list:
    return list(BASES[rng.integers(0, 4, size=LEN)])


def dinuc_bg(rng: np.random.Generator) -> list:
    out = [int(rng.choice(4, p=BASE_FREQS))]
    for _ in range(LEN - 1):
        out.append(int(rng.choice(4, p=DINUC_PROBS[out[-1]])))
    return [BASES[c] for c in out]


def insert_motifs(seq: list, pfms: list, rng: np.random.Generator, n: int) -> None:
    for _ in range(n):
        _, pfm = pfms[int(rng.integers(0, len(pfms)))]
        mlen = pfm.shape[1]
        if mlen > LEN:
            continue
        pos = int(rng.integers(0, LEN - mlen + 1))
        motif = sample_from_pfm(pfm, rng)
        for j, ch in enumerate(motif):
            seq[pos + j] = ch


def build(rng, pfms, bg_fn, n_motifs):
    seq = bg_fn(rng)
    insert_motifs(seq, pfms, rng, n_motifs)
    return "".join(seq)


def kmer_set(seq: str, k: int) -> set:
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}


def select_diverse(candidates: list[str], n_pick: int, k: int) -> list[int]:
    """Greedy selection: pick seq whose k-mer set adds the most NEW k-mers,
    weighted inversely by frequency in selected set."""
    n = len(candidates)
    print(f"Pre-computing {k}-mer sets for {n} candidates...")
    kmers_per = [kmer_set(s, k) for s in candidates]

    # Track how often each k-mer is "covered" in selection
    covered = Counter()
    picked = []
    picked_set = set()

    # Initial scores: novelty = number of unique k-mers in seq
    # Random initial choice
    rng = np.random.default_rng(SEED + 1)
    first = int(rng.integers(0, n))
    picked.append(first)
    picked_set.add(first)
    for km in kmers_per[first]:
        covered[km] += 1

    # For efficiency: precompute per-candidate kmer arrays as lists of int hashes
    # but simple approach: each round, score a SAMPLE of remaining candidates
    # and pick the best. This is O(n_pick * sample_size).
    sample_size = 200
    print(f"Greedy selecting {n_pick} of {n} (sample size {sample_size})...")

    for it in range(1, n_pick):
        if it % 5000 == 0:
            print(f"  {it}/{n_pick} picked, covered {len(covered)} unique {k}-mers")
        # sample candidates not yet picked
        cand_idx = rng.integers(0, n, size=sample_size)
        best_score = -1.0
        best_i = -1
        for ci in cand_idx:
            if int(ci) in picked_set:
                continue
            # score: sum of (1 / (1 + covered_count)) for each k-mer in seq
            kms = kmers_per[ci]
            s = 0.0
            for km in kms:
                s += 1.0 / (1.0 + covered[km])
            if s > best_score:
                best_score = s
                best_i = int(ci)
        if best_i < 0:
            # fallback: pick any unpicked
            for ci in range(n):
                if ci not in picked_set:
                    best_i = ci
                    break
        picked.append(best_i)
        picked_set.add(best_i)
        for km in kmers_per[best_i]:
            covered[km] += 1

    return picked


def main() -> None:
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")

    candidates = []
    print("Generating pool A: dense 8 overlap, random bg (020 design)...")
    for _ in range(N_PER_POOL):
        candidates.append(build(rng, pfms, random_bg, 8))
    print("Generating pool B: dense 10 overlap, random bg (022 design)...")
    for _ in range(N_PER_POOL):
        candidates.append(build(rng, pfms, random_bg, 10))
    print("Generating pool C: dense 8 overlap, mammalian dinuc bg (026 design)...")
    for _ in range(N_PER_POOL):
        candidates.append(build(rng, pfms, dinuc_bg, 8))

    print(f"Total candidates: {len(candidates)}")
    picked_idx = select_diverse(candidates, N_FINAL, K)
    seqs = [candidates[i] for i in picked_idx]
    assert len(seqs) == N_FINAL
    for s in seqs[:50]:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    gc = np.array([(s.count("G") + s.count("C")) / LEN for s in seqs[:5000]])
    print(f"GC: mean={gc.mean():.3f}, std={gc.std():.3f}")
    # Composition breakdown
    n_a = sum(1 for i in picked_idx if i < N_PER_POOL)
    n_b = sum(1 for i in picked_idx if N_PER_POOL <= i < 2*N_PER_POOL)
    n_c = sum(1 for i in picked_idx if i >= 2*N_PER_POOL)
    print(f"Picked: pool A (020-style)={n_a}, B (022-style)={n_b}, C (026-style)={n_c}")
    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")


if __name__ == "__main__":
    main()
