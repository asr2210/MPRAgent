"""
Experiment 001 — cCRE class-balanced library.

Sample 50,000 ENCODE cCREs stratified across classes:
- dELS, dELS+CTCF (distal enhancers)
- pELS, pELS+CTCF (proximal enhancers)
- PLS, PLS+CTCF (promoters)
- CTCF-only, CTCF+CTCF (insulators)
- DNase-H3K4me3, DNase-H3K4me3+CTCF

Extract 200bp centered on each cCRE midpoint, hg38 sequence.

Rationale: cCREs are the V3 ENCODE catalog (1.06M elements, union of accessibility +
chromatin signatures across 1518 biosamples). Class-balanced sampling promotes equal
representation of regulatory programs (promoters, enhancers, insulators), which we
hypothesize should generalize better than frequency-weighted sampling because the
held-out cell types may use different mixes of these programs.
"""
import os, sys, random
import pyfaidx

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)

SEQ_LEN = 200
TARGET_N = 50_000

# Read cCRE bed file: chrom, start, end, id1, id2, classification
ccre_by_class = {}
with open(os.path.join(DATA, "GRCh38-cCREs.bed")) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, start, end, _, _, classification = parts
        # Skip non-standard contigs
        if "_" in chrom or chrom in ("chrM", "chrEBV"):
            continue
        ccre_by_class.setdefault(classification, []).append(
            (chrom, int(start), int(end))
        )

print("cCRE class counts:")
for cls, lst in sorted(ccre_by_class.items(), key=lambda x: -len(x[1])):
    print(f"  {cls}: {len(lst)}")

n_classes = len(ccre_by_class)
per_class = TARGET_N // n_classes
remainder = TARGET_N - per_class * n_classes
print(f"\n{n_classes} classes, targeting {per_class} per class (+{remainder} extras)")

fa = pyfaidx.Fasta(os.path.join(DATA, "hg38.fa"))

selected = []
# Sort classes deterministically (largest first so extras go to largest classes)
sorted_classes = sorted(ccre_by_class.keys(), key=lambda c: -len(ccre_by_class[c]))
for i, cls in enumerate(sorted_classes):
    lst = ccre_by_class[cls]
    target = per_class + (1 if i < remainder else 0)
    if len(lst) <= target:
        chosen = lst[:]
    else:
        chosen = random.sample(lst, target)
    selected.extend([(cls, c, s, e) for (c, s, e) in chosen])

print(f"\nSelected {len(selected)} cCREs")

def get_seq(chrom, start, end):
    """Extract 200bp centered on midpoint of cCRE."""
    mid = (start + end) // 2
    s = mid - SEQ_LEN // 2
    e = s + SEQ_LEN
    if s < 0 or e > len(fa[chrom]):
        return None
    seq = str(fa[chrom][s:e]).upper()
    if len(seq) != SEQ_LEN:
        return None
    if "N" in seq:
        return None
    if not set(seq).issubset(set("ACGT")):
        return None
    return seq

# Extract sequences, replacing failures with new draws from same class
seqs = []
backlog = list(range(len(selected)))
random.shuffle(backlog)
attempt_pool = {cls: [c for c in ccre_by_class[cls]] for cls in sorted_classes}
for cls in sorted_classes:
    random.shuffle(attempt_pool[cls])

failed_indices = []
for idx in backlog:
    cls, chrom, s, e = selected[idx]
    seq = get_seq(chrom, s, e)
    if seq is None:
        failed_indices.append(idx)
        continue
    seqs.append(seq)

# Refill any failed sequences with fresh draws from the same class
print(f"\n{len(failed_indices)} failed; refilling")
fill_iter = {cls: iter(attempt_pool[cls]) for cls in sorted_classes}
for idx in failed_indices:
    cls = selected[idx][0]
    while True:
        try:
            chrom, s, e = next(fill_iter[cls])
        except StopIteration:
            # exhausted this class, fall back to any class
            for fb_cls in sorted_classes:
                try:
                    chrom, s, e = next(fill_iter[fb_cls])
                    break
                except StopIteration:
                    continue
            else:
                raise RuntimeError("ran out of cCREs")
        seq = get_seq(chrom, s, e)
        if seq is not None:
            seqs.append(seq)
            break

assert len(seqs) == TARGET_N, f"got {len(seqs)}, expected {TARGET_N}"

# Sanity-check: all 200bp ACGT
for i, s in enumerate(seqs):
    assert len(s) == SEQ_LEN, f"line {i}: len {len(s)}"
    assert set(s).issubset(set("ACGT")), f"line {i}: bad chars"

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"\nWrote {len(seqs)} sequences to {OUT}")
