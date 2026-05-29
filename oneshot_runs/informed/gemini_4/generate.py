import os
import gzip
import numpy as np
import random
import time
from twobitreader import TwoBitFile

# Set random seed for perfect reproducibility
np.random.seed(42)
random.seed(42)

start_time = time.time()

# Ensure library/ directory exists
os.makedirs("library", exist_ok=True)

# 1. Load reference genome
print("Loading hg38.2bit genome...")
genome = TwoBitFile("data/hg38.2bit")
print(f"Loaded genome in {time.time() - start_time:.2f}s")

# 2. Load NMF loadings (Mixture matrix)
load_start = time.time()
print("Loading NMF loadings (Mixture)...")
with gzip.open("data/2018-06-08NC16_NNDSVD_Mixture.npy.gz", "rb") as f:
    nmf_loadings = np.load(f)
print(f"Loaded NMF loadings with shape {nmf_loadings.shape} in {time.time() - load_start:.2f}s")

# 3. Load DHS metadata (chrom, start, end, summit, mean_signal)
meta_start = time.time()
print("Loading DHS metadata...")
num_elements = nmf_loadings.shape[1]
chroms = []
summits = np.zeros(num_elements, dtype=np.int32)
mean_signals = np.zeros(num_elements, dtype=np.float32)

with gzip.open("data/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz", "rt") as f:
    header = f.readline()
    for idx, line in enumerate(f):
        fields = line.strip().split("\t")
        if len(fields) < 10:
            continue
        chroms.append(fields[0])
        mean_signals[idx] = float(fields[4])
        summits[idx] = int(fields[6])
print(f"Loaded {len(chroms)} DHS elements in {time.time() - meta_start:.2f}s")

# Verification
assert len(chroms) == num_elements, "Metadata length does not match loadings!"

# Sampling parameters
total_library_size = 50000
num_topic_seqs = 35000  # 70% of 50000 (topic-specific, high-signal active elements)
num_random_seqs = 7500  # 15% of 50000 (diverse genomic chromatin elements)
num_synth_seqs = 7500   # 15% of 50000 (fully random synthetic background)

seen_seqs = set()
selected_sequences = []

# Helper to extract a safe 200bp sequence centered at the summit
def extract_sequence(chrom, summit):
    if chrom not in genome:
        return None
    # Center 200bp on summit
    start_pos = max(0, summit - 100)
    end_pos = start_pos + 200
    try:
        seq = genome[chrom][start_pos:end_pos].upper()
        if len(seq) == 200 and 'N' not in seq:
            return seq
    except Exception:
        pass
    return None

# --- Part A: Topic-Specific DHS Sampling (Cream of the Crop) ---
print("\nSampling topic-specific high-signal DHS elements...")
# Sample 35,000 sequences distributed over the 16 tissue/cell-type components (topics)
for k in range(16):
    needed = 2188 if k < 8 else 2187  # Total = 8*2188 + 8*2187 = 35000
    
    # Calculate weight: topic continuous loading * overall accessibility signal
    loadings = nmf_loadings[k, :]
    weights = loadings * mean_signals
    
    nonzero_indices = np.where(weights > 0)[0]
    sub_weights = weights[nonzero_indices]
    sub_weights /= np.sum(sub_weights)
    
    # Pre-draw candidates to ensure we have enough after filtering out 'N' and duplicates
    candidates = np.random.choice(nonzero_indices, size=int(needed * 2.5), replace=False, p=sub_weights)
    
    count = 0
    for idx_choice in candidates:
        if count >= needed:
            break
        seq = extract_sequence(chroms[idx_choice], summits[idx_choice])
        if seq and seq not in seen_seqs:
            seen_seqs.add(seq)
            selected_sequences.append(seq)
            count += 1
    print(f"Topic {k:02d}: Sampled {count}/{needed} sequences.")

# --- Part B: Broad Genomic DHS Sampling (Chromatin Diversity) ---
print("\nSampling diverse chromatin background elements...")
# Draw candidate indices uniformly from the entire DHS index pool
candidate_indices = np.random.choice(num_elements, size=num_random_seqs * 5, replace=False)

count = 0
for idx_choice in candidate_indices:
    if count >= num_random_seqs:
        break
    seq = extract_sequence(chroms[idx_choice], summits[idx_choice])
    if seq and seq not in seen_seqs:
        seen_seqs.add(seq)
        selected_sequences.append(seq)
        count += 1
print(f"Sampled {count}/{num_random_seqs} random diverse DHS sequences.")

# --- Part C: Synthetic Sequence Generation (Sequence Space Diversity) ---
print("\nGenerating synthetic background sequences...")
bases = ['A', 'C', 'G', 'T']
count = 0
while count < num_synth_seqs:
    seq = "".join(random.choices(bases, k=200))
    if seq not in seen_seqs:
        seen_seqs.add(seq)
        selected_sequences.append(seq)
        count += 1
print(f"Generated {count}/{num_synth_seqs} synthetic sequences.")

# Shuffle the final list to mix the different components nicely for training stability
print("\nShuffling final sequence list...")
random.shuffle(selected_sequences)

# Ensure exact count and validity of sequences
assert len(selected_sequences) == total_library_size, f"Incorrect library size! Got {len(selected_sequences)}."
for idx, seq in enumerate(selected_sequences):
    assert len(seq) == 200, f"Sequence at index {idx} has length {len(seq)} instead of 200!"
    assert all(b in bases for b in seq), f"Sequence at index {idx} contains non-ACGT bases: {seq}"

# Write to output file
out_path = "library/sequences.txt"
print(f"Writing sequences to {out_path}...")
with open(out_path, "w") as f:
    for seq in selected_sequences:
        f.write(seq + "\n")

print(f"\nSuccessfully generated {total_library_size} sequences in {time.time() - start_time:.2f}s!")
