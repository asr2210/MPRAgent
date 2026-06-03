"""
Experiment 030: Refined diversity selection (200k from 4 pools).

Builds on exp 029 breakthrough (mean_r 0.4288 vs prior plateau 0.4284).
Adds a 4th candidate pool (008-style: 5 motifs no-overlap, random bg)
to provide another short-context distribution. Uses larger candidate
pool (200k vs 150k) for more selection headroom.

Greedy 6-mer-novelty selection of 50k. Same k=6 and sample_size=200
as 029 (proven to work).

If this matches or exceeds 029 (eval_01 > 0.428), confirms T20 (diversity
selection robustly breaks plateau). The 4-pool variant should give the
selector more diversity to work with.
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
K = 6
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


def random_bg(rng):
    return list(BASES[rng.integers(0, 4, size=LEN)])


def dinuc_bg(rng):
    out = [int(rng.choice(4, p=BASE_FREQS))]
    for _ in range(LEN - 1):
        out.append(int(rng.choice(4, p=DINUC_PROBS[out[-1]])))
    return [BASES[c] for c in out]


def insert_overlap(seq, pfms, rng, n):
    for _ in range(n):
        _, pfm = pfms[int(rng.integers(0, len(pfms)))]
        mlen = pfm.shape[1]
        if mlen > LEN:
            continue
        pos = int(rng.integers(0, LEN - mlen + 1))
        motif = sample_from_pfm(pfm, rng)
        for j, ch in enumerate(motif):
            seq[pos + j] = ch


def insert_nooverlap(seq, pfms, rng, n):
    """008-style: try to place without overlap with previously placed motifs."""
    occupied = np.zeros(LEN, dtype=bool)
    for _ in range(n):
        _, pfm = pfms[int(rng.integers(0, len(pfms)))]
        mlen = pfm.shape[1]
        if mlen > LEN:
            continue
        for _try in range(20):
            pos = int(rng.integers(0, LEN - mlen + 1))
            if not occupied[pos:pos+mlen].any():
                motif = sample_from_pfm(pfm, rng)
                for j, ch in enumerate(motif):
                    seq[pos + j] = ch
                occupied[pos:pos+mlen] = True
                break


def kmer_set(seq, k):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}


def select_diverse(candidates, n_pick, k):
    n = len(candidates)
    print(f"Pre-computing {k}-mer sets for {n} candidates...")
    kmers_per = [kmer_set(s, k) for s in candidates]

    covered = Counter()
    picked = []
    picked_set = set()

    rng = np.random.default_rng(SEED + 1)
    first = int(rng.integers(0, n))
    picked.append(first)
    picked_set.add(first)
    for km in kmers_per[first]:
        covered[km] += 1

    sample_size = 200
    print(f"Greedy selecting {n_pick} of {n} (sample size {sample_size})...")
    for it in range(1, n_pick):
        if it % 5000 == 0:
            print(f"  {it}/{n_pick} picked, covered {len(covered)} unique {k}-mers")
        cand_idx = rng.integers(0, n, size=sample_size)
        best_score = -1.0
        best_i = -1
        for ci in cand_idx:
            if int(ci) in picked_set:
                continue
            kms = kmers_per[ci]
            s = 0.0
            for km in kms:
                s += 1.0 / (1.0 + covered[km])
            if s > best_score:
                best_score = s
                best_i = int(ci)
        if best_i < 0:
            for ci in range(n):
                if ci not in picked_set:
                    best_i = ci
                    break
        picked.append(best_i)
        picked_set.add(best_i)
        for km in kmers_per[best_i]:
            covered[km] += 1
    return picked


def main():
    rng = np.random.default_rng(SEED)
    pfms = parse_jaspar_filtered(JASPAR)
    print(f"Loaded {len(pfms)} targeted PFMs")

    candidates = []
    pool_starts = [0]

    print("Pool A: dense 8 overlap, random bg (020 design)...")
    for _ in range(N_PER_POOL):
        seq = random_bg(rng)
        insert_overlap(seq, pfms, rng, 8)
        candidates.append("".join(seq))
    pool_starts.append(len(candidates))

    print("Pool B: dense 10 overlap, random bg (022 design)...")
    for _ in range(N_PER_POOL):
        seq = random_bg(rng)
        insert_overlap(seq, pfms, rng, 10)
        candidates.append("".join(seq))
    pool_starts.append(len(candidates))

    print("Pool C: dense 8 overlap, mammalian dinuc bg (026 design)...")
    for _ in range(N_PER_POOL):
        seq = dinuc_bg(rng)
        insert_overlap(seq, pfms, rng, 8)
        candidates.append("".join(seq))
    pool_starts.append(len(candidates))

    print("Pool D: 5 motifs NO-overlap, random bg (008 design)...")
    for _ in range(N_PER_POOL):
        seq = random_bg(rng)
        insert_nooverlap(seq, pfms, rng, 5)
        candidates.append("".join(seq))
    pool_starts.append(len(candidates))

    print(f"Total candidates: {len(candidates)} from {len(pool_starts)-1} pools")
    picked_idx = select_diverse(candidates, N_FINAL, K)
    seqs = [candidates[i] for i in picked_idx]
    assert len(seqs) == N_FINAL
    for s in seqs[:50]:
        assert len(s) == LEN
        assert set(s).issubset(set("ACGT"))

    # Composition
    names = ["A_020", "B_022", "C_026", "D_008"]
    for i in range(len(pool_starts) - 1):
        lo, hi = pool_starts[i], pool_starts[i+1]
        cnt = sum(1 for j in picked_idx if lo <= j < hi)
        print(f"  Pool {names[i]}: picked {cnt}")

    gc = np.array([(s.count("G") + s.count("C")) / LEN for s in seqs[:5000]])
    print(f"GC: mean={gc.mean():.3f}, std={gc.std():.3f}")
    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {len(seqs)} to {OUT}")


if __name__ == "__main__":
    main()
