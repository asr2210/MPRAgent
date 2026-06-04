#!/usr/bin/env python3
import gzip
import os
import random
from collections import defaultdict

try:
    from twobitreader import TwoBitFile
except Exception:
    TwoBitFile = None


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "library", "sequences.txt")
CCRE_BED = os.path.join(DATA, "GRCh38-cCREs-V4.bed.gz")
HG38_2BIT = os.path.join(DATA, "hg38.2bit")
SEED = 20260522
LENGTH = 200
TARGET = 50000

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

MOTIFS = {
    "AP1": "TGASTCA",
    "ETS": "CCGGAAGT",
    "GATA": "WGATAR",
    "FOX": "TRTTTRY",
    "RUNX": "TGTGGT",
    "NFKB": "GGGRNNYYCC",
    "STAT": "TTCNNNGAA",
    "IRF": "GAAA",
    "ISRE": "GAAANNGAAA",
    "SMAD": "GTCTAGAC",
    "TEAD": "CATTCCA",
    "SOX": "AACAAT",
    "OCT": "ATGCAAAT",
    "CREB": "TGACGTCA",
    "CEBP": "TTGCNNAA",
    "KLFSP": "GGGCGG",
    "YY1": "CCATNTT",
    "NFY": "CCAAT",
    "EBOX": "CACGTG",
    "MYCMAX": "CACGTG",
    "HIF": "RCGTG",
    "NR": "RGKTCA",
    "ER": "GGTCANNNTGACC",
    "GR": "AGAACANNNTGTTCT",
    "HNF4": "RGGTCA",
    "HNF1": "GTTAATNATTAAC",
    "P53": "RRRCWWGYYY",
    "MEF2": "YTAWWWWTAR",
    "RFX": "GTNRCCNNRGYAAC",
    "REST": "TTCAGCACCACGGACAG",
    "CTCF": "CCASYAGRKGGCRS",
    "TATA": "TATAWAWR",
    "INR": "YYANWYY",
}

FAMILIES = {
    "immune": ["NFKB", "IRF", "ISRE", "STAT", "AP1", "RUNX"],
    "development": ["SOX", "OCT", "TEAD", "SMAD", "FOX", "GATA"],
    "metabolic": ["HNF4", "HNF1", "CEBP", "NR", "GR", "FOX"],
    "proliferation": ["ETS", "AP1", "MYCMAX", "EBOX", "CREB", "KLFSP"],
    "neural": ["REST", "SOX", "MEF2", "RFX", "CREB", "FOX"],
    "promoter": ["KLFSP", "NFY", "ETS", "YY1", "TATA", "INR"],
    "boundary": ["CTCF", "YY1", "RFX", "KLFSP"],
}

CLASS_COUNTS = {
    "dELS": 10000,
    "pELS": 5000,
    "PLS": 5000,
    "CA": 5000,
    "CA-CTCF": 4000,
    "TF": 3000,
    "CA-H3K4me3": 2500,
    "CA-TF": 1500,
}


def gc_random(rng, n, gc):
    seq = []
    for _ in range(n):
        if rng.random() < gc:
            seq.append("G" if rng.random() < 0.5 else "C")
        else:
            seq.append("A" if rng.random() < 0.5 else "T")
    return "".join(seq)


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def instantiate(consensus, rng):
    return "".join(rng.choice(IUPAC.get(ch, "ACGT")) for ch in consensus.upper())


def mutate(seq, rng, rate=0.18):
    bases = "ACGT"
    out = []
    for ch in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in bases if b != ch]))
        else:
            out.append(ch)
    return "".join(out)


def place(seq, insert, pos):
    return seq[:pos] + insert + seq[pos + len(insert):]


def valid(seq):
    return len(seq) == LENGTH and set(seq) <= set("ACGT")


def open_genome():
    if TwoBitFile is None or not os.path.exists(HG38_2BIT):
        return None
    return TwoBitFile(HG38_2BIT)


def chrom_sizes(genome):
    sizes = {}
    for chrom in genome.keys():
        if chrom.startswith("chr") and "_" not in chrom and chrom not in {"chrM"}:
            try:
                sizes[chrom] = len(genome[chrom])
            except Exception:
                pass
    return sizes


def fetch_window(genome, chrom, center, rng, jitter=90):
    if chrom not in genome:
        return None
    chrom_len = len(genome[chrom])
    c = center + rng.randint(-jitter, jitter)
    start = max(0, min(chrom_len - LENGTH, c - LENGTH // 2))
    try:
        seq = genome[chrom][start:start + LENGTH].upper()
    except Exception:
        return None
    if valid(seq):
        return seq
    return None


def load_ccres():
    by_class = defaultdict(list)
    if not os.path.exists(CCRE_BED):
        return by_class
    with gzip.open(CCRE_BED, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, start, end, cls = fields[0], int(fields[1]), int(fields[2]), fields[5]
            if chrom.startswith("chr") and "_" not in chrom and chrom != "chrM":
                by_class[cls].append((chrom, start, end))
    return by_class


def genomic_ccre_sequences(rng, genome, count):
    by_class = load_ccres()
    if not by_class or genome is None:
        return []
    seqs = []
    for cls, n in CLASS_COUNTS.items():
        pool = by_class.get(cls, [])
        if not pool:
            continue
        tries = 0
        while n > 0 and tries < n * 40:
            tries += 1
            chrom, start, end = rng.choice(pool)
            seq = fetch_window(genome, chrom, (start + end) // 2, rng)
            if seq:
                if rng.random() < 0.5:
                    seq = revcomp(seq)
                seqs.append(seq)
                n -= 1
    rng.shuffle(seqs)
    return seqs[:count]


def random_genomic_sequences(rng, genome, count):
    if genome is None:
        return []
    sizes = chrom_sizes(genome)
    chroms = list(sizes)
    weights = [sizes[c] for c in chroms]
    seqs = []
    tries = 0
    while len(seqs) < count and tries < count * 100:
        tries += 1
        chrom = rng.choices(chroms, weights=weights, k=1)[0]
        start = rng.randint(0, sizes[chrom] - LENGTH)
        seq = genome[chrom][start:start + LENGTH].upper()
        if valid(seq):
            if rng.random() < 0.5:
                seq = revcomp(seq)
            seqs.append(seq)
    return seqs


def synthetic_enhancer(rng, idx):
    gc = rng.choice([0.32, 0.38, 0.45, 0.52, 0.60, 0.68])
    seq = gc_random(rng, LENGTH, gc)
    family = rng.choice(list(FAMILIES))
    motifs = list(FAMILIES[family])
    rng.shuffle(motifs)
    copies = rng.randint(2, 7)
    occupied = []
    for i in range(copies):
        name = motifs[i % len(motifs)] if rng.random() < 0.7 else rng.choice(list(MOTIFS))
        motif = instantiate(MOTIFS[name], rng)
        if rng.random() < 0.5:
            motif = revcomp(motif)
        if idx % 11 == 0:
            motif = mutate(motif, rng, 0.22)
        for _ in range(80):
            pos = rng.randint(8, LENGTH - len(motif) - 8)
            if all(abs(pos - old) > 5 for old in occupied):
                seq = place(seq, motif, pos)
                occupied.append(pos)
                break
    return seq


def synthetic_promoter(rng, idx):
    gc = rng.choice([0.48, 0.56, 0.64, 0.72])
    seq = gc_random(rng, LENGTH, gc)
    if rng.random() < 0.75:
        tata = instantiate(MOTIFS["TATA"], rng)
        seq = place(seq, tata, rng.randint(60, 85))
    inr = instantiate(MOTIFS["INR"], rng)
    seq = place(seq, inr, rng.randint(92, 103))
    for name in rng.sample(["KLFSP", "NFY", "ETS", "YY1", "CREB", "EBOX"], rng.randint(2, 5)):
        motif = instantiate(MOTIFS[name], rng)
        if rng.random() < 0.35:
            motif = revcomp(motif)
        pos = rng.choice([rng.randint(10, 55), rng.randint(120, 185 - len(motif))])
        if idx % 13 == 0:
            motif = mutate(motif, rng, 0.18)
        seq = place(seq, motif, pos)
    return seq


def synthetic_ctcf(rng, idx):
    seq = gc_random(rng, LENGTH, rng.choice([0.38, 0.45, 0.52, 0.60]))
    motif = instantiate(MOTIFS["CTCF"], rng)
    if rng.random() < 0.5:
        motif = revcomp(motif)
    if idx % 7 == 0:
        motif = mutate(motif, rng, 0.20)
    seq = place(seq, motif, rng.randint(75, 110))
    for name in rng.sample(["YY1", "RFX", "KLFSP", "ETS"], rng.randint(0, 2)):
        m = instantiate(MOTIFS[name], rng)
        if rng.random() < 0.5:
            m = revcomp(m)
        seq = place(seq, m, rng.choice([rng.randint(10, 55), rng.randint(135, 185 - len(m))]))
    return seq


def synthetic_background(rng, idx):
    seq = gc_random(rng, LENGTH, rng.choice([0.25, 0.32, 0.40, 0.50, 0.60, 0.70]))
    if idx % 3 == 0:
        tract = rng.choice(["A", "T", "C", "G"]) * rng.randint(6, 14)
        seq = place(seq, tract, rng.randint(0, LENGTH - len(tract)))
    if idx % 5 == 0:
        motif = mutate(instantiate(MOTIFS[rng.choice(list(MOTIFS))], rng), rng, 0.35)
        seq = place(seq, motif, rng.randint(0, LENGTH - len(motif)))
    return seq


def synthetic_sequences(rng, count):
    seqs = []
    for i in range(count):
        bucket = i % 8
        if bucket in {0, 1, 2, 3}:
            seq = synthetic_enhancer(rng, i)
        elif bucket in {4, 5}:
            seq = synthetic_promoter(rng, i)
        elif bucket == 6:
            seq = synthetic_ctcf(rng, i)
        else:
            seq = synthetic_background(rng, i)
        seqs.append(seq)
    return seqs


def dedupe_fill(seqs, rng):
    seen = set()
    out = []
    for seq in seqs:
        if valid(seq) and seq not in seen:
            out.append(seq)
            seen.add(seq)
    while len(out) < TARGET:
        seq = synthetic_enhancer(rng, len(out))
        if seq not in seen:
            out.append(seq)
            seen.add(seq)
    return out[:TARGET]


def main():
    rng = random.Random(SEED)
    genome = open_genome()
    seqs = []
    seqs.extend(genomic_ccre_sequences(rng, genome, 36000))
    seqs.extend(random_genomic_sequences(rng, genome, 6000))
    seqs.extend(synthetic_sequences(rng, 8000))
    if len(seqs) < TARGET:
        seqs.extend(synthetic_sequences(rng, TARGET - len(seqs)))
    rng.shuffle(seqs)
    seqs = dedupe_fill(seqs, rng)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        for seq in seqs:
            fh.write(seq + "\n")
    assert len(seqs) == TARGET
    assert all(valid(seq) for seq in seqs)


if __name__ == "__main__":
    main()
