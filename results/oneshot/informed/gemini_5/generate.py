import os
import sys
import random
import bisect
import collections
import time
import pyfaidx

def main():
    # Set random seed for reproducibility
    random.seed(42)

    print("Starting library generation...")
    start_time = time.time()

    # Define primary chromosomes to use
    PRIMARY_CHROMS = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])

    # 1. Parse chromosome lengths from hg38.fa.fai
    chrom_lengths = {}
    with open("data/hg38.fa.fai", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] in PRIMARY_CHROMS:
                chrom_lengths[parts[0]] = int(parts[1])

    print(f"Loaded {len(chrom_lengths)} chromosome lengths.")

    # 2. Parse cCRE BED file and group by category and chromosome
    # Categories of interest: PLS, pELS, dELS, CTCF-only, DNase-H3K4me3
    # We will keep track of cCRE intervals for overlap checking in genomic negatives
    ccre_by_chrom = collections.defaultdict(list)
    ccre_by_category = collections.defaultdict(list)

    print("Parsing ENCODE cCRE BED file...")
    ccre_parse_start = time.time()
    with open("data/encodeCcreCombined.bed", "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            if chrom not in PRIMARY_CHROMS:
                continue
            
            start = int(parts[1])
            end = int(parts[2])
            category = parts[10] # e.g. PLS, pELS, dELS, CTCF-only, DNase-H3K4me3
            
            # Store coordinate for overlap checking
            ccre_by_chrom[chrom].append((start, end))
            
            # Store for category sampling
            ccre_by_category[category].append((chrom, start, end))

    print(f"Parsed BED file in {time.time() - ccre_parse_start:.2f} seconds.")
    for cat, items in ccre_by_category.items():
        print(f"  Category '{cat}': {len(items)} elements")

    # 3. Sort and merge cCRE intervals per chromosome for efficient overlap checking
    # Sorting is required for binary search / bisect
    print("Sorting cCRE intervals per chromosome...")
    for chrom in ccre_by_chrom:
        ccre_by_chrom[chrom].sort()

    # Define overlap check function using binary search
    def check_overlap(chrom, interval_start, interval_end):
        intervals = ccre_by_chrom.get(chrom, [])
        if not intervals:
            return False
        
        # Binary search
        idx = bisect.bisect_right(intervals, (interval_start, float('inf'))) - 1
        if idx >= 0:
            if intervals[idx][1] > interval_start:
                return True
        if idx + 1 < len(intervals):
            if intervals[idx + 1][0] < interval_end:
                return True
        return False

    # 4. Open genome FASTA
    print("Opening hg38.fa genome...")
    genome = pyfaidx.Fasta('data/hg38.fa')

    # Target counts for active elements
    target_counts = {
        "PLS": 10000,
        "pELS": 10000,
        "dELS": 10000,
        "CTCF-only": 3000,
        "DNase-H3K4me3": 2000
    }

    sequences = []
    metadata = [] # To document the origin of each sequence

    # Helper function to validate sequence (only {A, C, G, T} and length 200)
    def validate_seq(seq):
        if len(seq) != 200:
            return False
        return all(c in "ACGT" for c in seq)

    # 5. Sample active elements from each cCRE category
    for category, target_count in target_counts.items():
        print(f"Sampling {target_count} elements from category '{category}'...")
        candidates = ccre_by_category.get(category, [])
        
        # Shuffle candidates to get a diverse random selection
        random.shuffle(candidates)
        
        count = 0
        skipped_non_acgt = 0
        
        for chrom, start, end in candidates:
            if count >= target_count:
                break
                
            # Get center of cCRE and extract exactly 200bp
            center = (start + end) // 2
            seq_start = center - 100
            seq_end = center + 100
            
            # Check boundaries
            if seq_start < 0 or seq_end > chrom_lengths[chrom]:
                continue
                
            # Query sequence
            seq = str(genome[chrom][seq_start:seq_end]).upper()
            
            if validate_seq(seq):
                sequences.append(seq)
                metadata.append(f"cCRE|{category}|{chrom}:{seq_start}-{seq_end}")
                count += 1
            else:
                skipped_non_acgt += 1
                
        print(f"  Successfully sampled {count} elements (skipped {skipped_non_acgt} containing non-ACGT bases).")
        if count < target_count:
            print(f"WARNING: Could only sample {count} elements instead of {target_count} target!")

    # 6. Sample Inactive Elements (Genomic Negatives)
    print("Sampling 10,000 genomic negatives...")
    neg_target = 10000
    neg_count = 0
    neg_skipped_overlap = 0
    neg_skipped_non_acgt = 0

    # Build chromosome selection distribution based on chromosome length
    chrom_list = list(chrom_lengths.keys())
    chrom_weights = [chrom_lengths[c] for c in chrom_list]

    while neg_count < neg_target:
        # Sample chromosome proportional to its length
        chrom = random.choices(chrom_list, weights=chrom_weights, k=1)[0]
        length = chrom_lengths[chrom]
        
        # Random start coordinate
        start = random.randint(1000, length - 1200)
        end = start + 200
        
        # Check overlap with any cCRE on this chromosome
        if check_overlap(chrom, start, end):
            neg_skipped_overlap += 1
            continue
            
        # Extract sequence
        seq = str(genome[chrom][start:end]).upper()
        
        if validate_seq(seq):
            sequences.append(seq)
            metadata.append(f"genomic_negative|{chrom}:{start}-{end}")
            neg_count += 1
        else:
            neg_skipped_non_acgt += 1

    print(f"Successfully sampled {neg_count} genomic negatives.")
    print(f"  (Skipped {neg_skipped_overlap} due to cCRE overlap, {neg_skipped_non_acgt} due to non-ACGT bases).")

    # 7. Generate Synthetic Sequences
    # We want 5,000 synthetic sequences in total:
    # - 2,500 random synthetic
    # - 2,500 motif-inserted synthetic
    print("Generating 5,000 synthetic sequences...")
    
    MOTIFS = {
        "AP-1": ["TGAGTCA", "TGACTCA"],
        "CTCF": ["CCACCAGGGGGCGGC", "GCCGCCCCCTGGTGG"],  # Forward and RC
        "Sp1": ["GGGCGG", "CCGCCC"],
        "NF-kB": ["GGGAATTTCC", "GGAAATTCCC"],
        "TATA": ["TATAAA", "TTTATA"],
        "YY1": ["CCGCCATTTT", "AAAATGGCGG"],
        "CREB": ["TGACGTCA"]
    }

    def generate_synthetic(target_gc, insert_motifs):
        p_gc = target_gc / 2.0
        p_at = (1.0 - target_gc) / 2.0
        bases = ['A', 'C', 'G', 'T']
        weights = [p_at, p_gc, p_gc, p_at]
        
        seq_chars = random.choices(bases, weights=weights, k=200)
        
        if insert_motifs:
            # Number of motifs to insert (1 to 3)
            num_motifs = random.randint(1, 3)
            selected_motif_names = random.sample(list(MOTIFS.keys()), num_motifs)
            
            occupied = []
            for name in selected_motif_names:
                motif_seq = random.choice(MOTIFS[name])
                m_len = len(motif_seq)
                
                # Try to place safely without overlap and away from edges
                placed = False
                for _ in range(50):
                    pos = random.randint(10, 200 - m_len - 10)
                    overlap = False
                    for o_start, o_end in occupied:
                        if not (pos + m_len <= o_start or pos >= o_end):
                            overlap = True
                            break
                    if not overlap:
                        seq_chars[pos:pos+m_len] = list(motif_seq)
                        occupied.append((pos, pos+m_len))
                        placed = True
                        break
                        
        return "".join(seq_chars)

    # Generate 2,500 pure random synthetic
    syn_rand_count = 0
    while syn_rand_count < 2500:
        # Vary GC content between 35% and 65%
        gc = random.uniform(0.35, 0.65)
        seq = generate_synthetic(gc, insert_motifs=False)
        if validate_seq(seq):
            sequences.append(seq)
            metadata.append(f"synthetic|random|gc:{gc:.3f}")
            syn_rand_count += 1

    # Generate 2,500 motif-inserted synthetic
    syn_motif_count = 0
    while syn_motif_count < 2500:
        # Vary GC content between 35% and 65%
        gc = random.uniform(0.35, 0.65)
        seq = generate_synthetic(gc, insert_motifs=True)
        if validate_seq(seq):
            sequences.append(seq)
            metadata.append(f"synthetic|motif_inserted|gc:{gc:.3f}")
            syn_motif_count += 1

    print(f"Generated {syn_rand_count} random synthetic and {syn_motif_count} motif-inserted synthetic sequences.")

    # 8. Final verification
    print(f"Total sequences generated: {len(sequences)}")
    assert len(sequences) == 50000, f"Error: Generated {len(sequences)} sequences instead of 50,000!"

    # Write sequences to sequences.txt
    output_dir = "library"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "sequences.txt")
    
    with open(output_path, "w") as f_out:
        for seq in sequences:
            f_out.write(seq + "\n")
            
    # Write metadata map for tracking/analysis
    metadata_path = os.path.join(output_dir, "metadata.txt")
    with open(metadata_path, "w") as f_meta:
        for meta in metadata:
            f_meta.write(meta + "\n")

    print(f"Successfully wrote 50,000 sequences to {output_path}")
    print(f"Successfully wrote sequence metadata to {metadata_path}")
    print(f"Total processing completed in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()
