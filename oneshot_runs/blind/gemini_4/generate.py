#!/usr/bin/env python3
"""
generate.py
Author: Autonomous CLI Agent (Gemini CLI)
Date: May 27, 2026

Generates exactly 50,000 high-quality, diverse 200bp sequences designed
from first-principles to train the best possible general-purpose model
of gene regulatory activity across all human cell lines.

The library contains a balanced combination of:
1. Natural Genomic Sequences (16,000 WT)
   - 4,000 Promoter-like (PLS)
   - 4,000 Proximal Enhancer-like (pELS)
   - 4,000 Distal Enhancer-like (dELS)
   - 2,000 Chromatin Accessible CTCF-only (Insulators)
   - 2,000 Neutral Genomic backgrounds (non-functional regions)
2. Synthetic Motif Grammar (16,000 WT)
   - 5,333 Motif Scanning / Positional effect templates
   - 5,333 Homotypic cooperative clusters (1-4 copies, varying spacing)
   - 5,334 Heterotypic cooperative pairs (varying spacing/orientations)
3. Pure Random Backgrounds (3,000)
   - Uniform, continuous distribution of GC content (30% to 75%)
4. Mutation & Perturbation Pairs (15,000 Mutants)
   - 7,500 Natural WT mutants (point mutations & knockout scrambles)
   - 7,500 Synthetic WT mutants (point mutations & knockout scrambles)
   
Total = 50,000 sequences. Saved directly to library/sequences.txt.
Full metadata saved to library/metadata.tsv for documentation.
"""

import os
import json
import gzip
import random
import bisect
import collections

# Set fixed seed for perfect reproducibility
random.seed(42)

# --- Define Constants ---
OUTPUT_DIR = "library"
SEQUENCES_FILE = os.path.join(OUTPUT_DIR, "sequences.txt")
METADATA_FILE = os.path.join(OUTPUT_DIR, "metadata.tsv")

# Core Transcription Factor motifs from major families
MOTIFS = {
    'AP-1': 'TGAGTCA',           # Basic leucine zipper (ubiquitous, pioneer/helper)
    'GATA': 'TGATAG',            # GATA zinc finger (erythroid/K562 and lineage regulator)
    'HNF4A': 'TGACCTTGACT',      # Nuclear receptor (liver/HepG2 specific)
    'SOX9': 'CATTGTT',           # SOX/HMG box (neuronal/SK-N-SH and stem cell specific)
    'SP1': 'GGGCGG',             # GC-box, classical promoter element
    'NF-Y': 'CCAAT',             # CCAAT-box, classical promoter element
    'YY1': 'CCATATT',            # Initiator and ubiquitous chromatin loop driver
    'CTCF': 'CCACCAGGGGGGG',     # CCCTC-binding factor, insulator / loop anchor
    'TATA': 'TATAAA',            # TATA box promoter core
    'CREB': 'TGACGTCA'           # cAMP response element, activity-dependent factor
}

# --- Helper Functions ---

def rev_comp(seq):
    """Returns the reverse complement of a DNA sequence."""
    rc_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(rc_dict[x] for x in reversed(seq))

def load_chromosomes():
    """Loads FASTA sequences for chr20, chr21, and chr22 from data folder."""
    chroms = {}
    for chrom in ['chr20', 'chr21', 'chr22']:
        filepath = f"data/{chrom}.fa.gz"
        print(f"Loading sequence for {chrom} from {filepath}...")
        with gzip.open(filepath, "rt") as f:
            f.readline()  # Skip header
            seq = "".join(line.strip() for line in f).upper()
        print(f"Loaded {chrom} with length {len(seq)} bp.")
        chroms[chrom] = seq
    return chroms

def load_ccre_registry():
    """Loads pre-fetched ENCODE cCRE registry from data folder."""
    filepath = "data/ccres.json"
    print(f"Loading cCRE registry from {filepath}...")
    with open(filepath) as f:
        return json.load(f)

def generate_random_bg(length, gc_content):
    """Generates a random DNA sequence of target GC content."""
    p_gc = gc_content
    p_at = 1.0 - p_gc
    choices = ['A', 'C', 'G', 'T']
    weights = [p_at / 2, p_gc / 2, p_gc / 2, p_at / 2]
    return "".join(random.choices(choices, weights=weights, k=length))

def scramble_subseq(subseq):
    """Scrambles a sub-sequence to destroy motif signature while retaining nucleotide composition."""
    chars = list(subseq)
    random.shuffle(chars)
    scrambled = "".join(chars)
    if scrambled == subseq:
        # If scramble didn't change anything, mutate characters
        chars = list(subseq)
        for i in range(len(chars)):
            chars[i] = random.choice([x for x in 'ACGT' if x != chars[i]])
        scrambled = "".join(chars)
    return scrambled

def point_mutate_seq(seq, pos):
    """Changes a single nucleotide at a specific position to a different one."""
    chars = list(seq)
    orig = chars[pos]
    choices = [x for x in 'ACGT' if x != orig]
    chars[pos] = random.choice(choices)
    return "".join(chars)

def find_motif_matches(seq, motifs):
    """Finds all forward and reverse complement matches of known motifs in a sequence."""
    matches = []
    seq_upper = seq.upper()
    for name, consensus in motifs.items():
        # Forward search
        start = 0
        while True:
            pos = seq_upper.find(consensus, start)
            if pos == -1:
                break
            matches.append((pos, pos + len(consensus), name, False))
            start = pos + 1
        # Reverse complement search
        rc_consensus = rev_comp(consensus)
        if rc_consensus != consensus:  # avoid duplicate counts for palindromic motifs
            start = 0
            while True:
                pos = seq_upper.find(rc_consensus, start)
                if pos == -1:
                    break
                matches.append((pos, pos + len(rc_consensus), name, True))
                start = pos + 1
    return matches

# --- Sampling Logic ---

def sample_natural_ccres(ccres_data, chrom_seqs, ccre_class_filter, target_count):
    """Samples unique, non-overlapping, valid 200bp sequences for specific cCRE classes."""
    pool = []
    # Combine cCREs from all three chromosomes matching the class filter
    for chrom, items in ccres_data.items():
        for c in items:
            if c.get('cCRE_class') in ccre_class_filter:
                pool.append((chrom, c['chromStart'], c['chromEnd'], c['name']))
                
    print(f"Total candidate pool size for filter {ccre_class_filter}: {len(pool)}")
    random.shuffle(pool)
    
    extracted_seqs = []
    used_positions = collections.defaultdict(list) # prevent duplicate/heavy overlapping coordinates
    
    for chrom, start, end, name in pool:
        if len(extracted_seqs) >= target_count:
            break
            
        center = (start + end) // 2
        seq_start = center - 100
        seq_end = center + 100
        
        # Check boundary
        if seq_start < 0 or seq_end > len(chrom_seqs[chrom]):
            continue
            
        # Check if overlaps with an already selected sequence on the same chromosome by more than 50bp
        overlap = False
        for s, e in used_positions[chrom]:
            if max(s, seq_start) < min(e, seq_end):
                overlap = True
                break
        if overlap:
            continue
            
        seq = chrom_seqs[chrom][seq_start:seq_end]
        
        # Verify characters
        if len(seq) == 200 and all(x in 'ACGT' for x in seq):
            extracted_seqs.append({
                'sequence': seq,
                'class': ccre_class_filter[0],
                'details': f"chrom={chrom};start={seq_start};end={seq_end};ccre={name}"
            })
            used_positions[chrom].append((seq_start, seq_end))
            
    print(f"Successfully extracted {len(extracted_seqs)} valid genomic sequences.")
    assert len(extracted_seqs) == target_count, f"Could not extract enough sequences for filter {ccre_class_filter}"
    return extracted_seqs

def sample_neutral_genomic(ccres_data, chrom_seqs, target_count):
    """Samples non-functional genomic regions (at least 5kb away from any cCRE)."""
    print(f"Sampling {target_count} neutral genomic sequences...")
    # Build sorted forbidden intervals
    forbidden = collections.defaultdict(list)
    for chrom, items in ccres_data.items():
        for c in items:
            forbidden[chrom].append((c['chromStart'] - 5000, c['chromEnd'] + 5000))
            
    # Merge overlapping forbidden intervals for extremely fast binary search
    merged_forbidden = {}
    for chrom in forbidden:
        forbidden[chrom].sort()
        merged = []
        for start, end in forbidden[chrom]:
            if not merged or merged[-1][1] < start:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        merged_forbidden[chrom] = merged

    def is_forbidden_overlap(chrom, start, end):
        intervals = merged_forbidden[chrom]
        # Binary search
        idx = bisect.bisect_right(intervals, [start, float('inf')])
        if idx > 0 and intervals[idx-1][1] >= start:
            return True
        if idx < len(intervals) and intervals[idx][0] <= end:
            return True
        return False

    extracted_seqs = []
    chromosomes = ['chr20', 'chr21', 'chr22']
    attempts = 0
    
    while len(extracted_seqs) < target_count and attempts < 200000:
        attempts += 1
        chrom = random.choice(chromosomes)
        seq_len = len(chrom_seqs[chrom])
        start = random.randint(10000, seq_len - 10200)
        end = start + 200
        
        if not is_forbidden_overlap(chrom, start, end):
            seq = chrom_seqs[chrom][start:end]
            if len(seq) == 200 and all(x in 'ACGT' for x in seq):
                extracted_seqs.append({
                    'sequence': seq,
                    'class': 'NeutralGenomic',
                    'details': f"chrom={chrom};start={start};end={end}"
                })
                # Add to forbidden dynamically so we don't select the same region
                merged_forbidden[chrom].append([start - 2000, end + 2000])
                merged_forbidden[chrom].sort() # keep sorted
                
    print(f"Successfully extracted {len(extracted_seqs)} neutral genomic sequences in {attempts} attempts.")
    assert len(extracted_seqs) == target_count, "Could not extract enough neutral genomic sequences"
    return extracted_seqs

# --- Synthetic Library Generation ---

def create_scanning_seqs(count):
    """Creates synthetic scanning sequences: a single motif is systematically placed along varying GC backgrounds."""
    print(f"Creating {count} synthetic motif scanning sequences...")
    seqs = []
    motifs_list = list(MOTIFS.items())
    gc_choices = [0.35, 0.45, 0.55, 0.65]
    
    attempts = 0
    while len(seqs) < count and attempts < 100000:
        attempts += 1
        name, consensus = random.choice(motifs_list)
        motif_len = len(consensus)
        
        # Select start position with some buffer
        start_pos = random.randint(10, 200 - 10 - motif_len)
        gc = random.choice(gc_choices)
        
        # Decide orientation
        motif_instance = consensus if random.random() < 0.5 else rev_comp(consensus)
        
        # Generate background
        bg = generate_random_bg(200, gc)
        # Implant motif
        seq = bg[:start_pos] + motif_instance + bg[start_pos + motif_len:]
        
        seqs.append({
            'sequence': seq,
            'class': 'Scanning',
            'details': f"motif={name};start={start_pos};gc={gc};rc={motif_instance != consensus}"
        })
    return seqs

def create_homotypic_seqs(count):
    """Creates homotypic clusters: 2 to 4 copies of the same motif with varying spacing."""
    print(f"Creating {count} synthetic homotypic cooperative sequences...")
    seqs = []
    motifs_list = list(MOTIFS.items())
    gc_choices = [0.35, 0.45, 0.55, 0.65]
    spacing_choices = [5, 10, 15, 20, 25, 30]
    
    attempts = 0
    while len(seqs) < count and attempts < 100000:
        attempts += 1
        name, consensus = random.choice(motifs_list)
        motif_len = len(consensus)
        
        num_copies = random.choice([2, 3, 4])
        spacing = random.choice(spacing_choices)
        
        total_len = num_copies * motif_len + (num_copies - 1) * spacing
        if total_len > 180:
            continue
            
        start_pos = random.randint(10, 200 - 10 - total_len)
        gc = random.choice(gc_choices)
        bg = generate_random_bg(200, gc)
        
        segment_parts = []
        for i in range(num_copies):
            # alternating orientations possible
            motif_instance = consensus if random.random() < 0.5 else rev_comp(consensus)
            segment_parts.append(motif_instance)
            if i < num_copies - 1:
                segment_parts.append(generate_random_bg(spacing, gc))
                
        segment = "".join(segment_parts)
        seq = bg[:start_pos] + segment + bg[start_pos + len(segment):]
        # Pad or truncate to ensure exactly 200bp
        if len(seq) != 200:
            seq = (seq + generate_random_bg(200, gc))[:200]
            
        seqs.append({
            'sequence': seq,
            'class': 'Homotypic',
            'details': f"motif={name};copies={num_copies};spacing={spacing};gc={gc}"
        })
    return seqs

def create_heterotypic_seqs(count):
    """Creates heterotypic pairs: two different motifs with varying spacing and orientations."""
    print(f"Creating {count} synthetic heterotypic cooperative sequences...")
    seqs = []
    motifs_list = list(MOTIFS.items())
    gc_choices = [0.35, 0.45, 0.55, 0.65]
    spacing_choices = [5, 10, 15, 20, 25, 30]
    
    attempts = 0
    while len(seqs) < count and attempts < 100000:
        attempts += 1
        m1_name, m1_consensus = random.choice(motifs_list)
        m2_name, m2_consensus = random.choice([x for x in motifs_list if x[0] != m1_name])
        
        m1_instance = m1_consensus if random.random() < 0.5 else rev_comp(m1_consensus)
        m2_instance = m2_consensus if random.random() < 0.5 else rev_comp(m2_consensus)
        
        spacing = random.choice(spacing_choices)
        total_len = len(m1_instance) + spacing + len(m2_instance)
        if total_len > 180:
            continue
            
        start_pos = random.randint(10, 200 - 10 - total_len)
        gc = random.choice(gc_choices)
        bg = generate_random_bg(200, gc)
        
        segment = m1_instance + generate_random_bg(spacing, gc) + m2_instance
        seq = bg[:start_pos] + segment + bg[start_pos + len(segment):]
        if len(seq) != 200:
            seq = (seq + generate_random_bg(200, gc))[:200]
            
        seqs.append({
            'sequence': seq,
            'class': 'Heterotypic',
            'details': f"motifs={m1_name}+{m2_name};spacing={spacing};gc={gc}"
        })
    return seqs

# --- Mutational / Perturbation Processing ---

def mutate_single_seq(seq):
    """Mutates a sequence by either scrambling a matched motif (knockout) or introducing a point mutation."""
    matches = find_motif_matches(seq, MOTIFS)
    
    if matches:
        # Mutate a randomly matched motif
        match = random.choice(matches)
        start, end, name, is_rc = match
        
        if random.random() < 0.5:
            # Knockout / Scramble the motif
            scrambled = scramble_subseq(seq[start:end])
            mutated_seq = seq[:start] + scrambled + seq[end:]
            mtype = "knockout"
            details = f"target_motif={name};pos={start}-{end};type=knockout"
        else:
            # Point mutation at center of motif
            mid_pos = (start + end) // 2
            mutated_seq = point_mutate_seq(seq, mid_pos)
            mtype = "point"
            details = f"target_motif={name};pos={mid_pos};type=point;orig={seq[mid_pos]};new={mutated_seq[mid_pos]}"
    else:
        # No motif matched, mutate the center 100bp of the sequence
        if random.random() < 0.5:
            # Knockout/scramble a 10bp window at the center
            scrambled = scramble_subseq(seq[95:105])
            mutated_seq = seq[:95] + scrambled + seq[105:]
            mtype = "scramble_center"
            details = f"target_motif=None;pos=95-105;type=scramble_center"
        else:
            # Point mutation at a random center position
            mut_pos = random.randint(50, 150)
            mutated_seq = point_mutate_seq(seq, mut_pos)
            mtype = "point_center"
            details = f"target_motif=None;pos={mut_pos};type=point_center;orig={seq[mut_pos]};new={mutated_seq[mut_pos]}"
            
    assert len(mutated_seq) == 200, f"Mutation output length mismatch: {len(mutated_seq)}"
    return mutated_seq, mtype, details


# --- Main Orchestration ---

def main():
    print("="*60)
    print("      Massively Parallel Reporter Assay (MPRA) Library Generator")
    print("="*60)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load genomic datasets
    chrom_seqs = load_chromosomes()
    ccres_data = load_ccre_registry()
    
    library = []
    
    # -------------------------------------------------------------
    # PART 1: Natural Genomic Sequences (16,000 sequences)
    # -------------------------------------------------------------
    print("\n--- Generating PART 1: Natural Genomic Sequences ---")
    
    # 1.1 Promoter-like (PLS): 4,000 sequences
    pls_seqs = sample_natural_ccres(ccres_data, chrom_seqs, ['Promoter', 'CA-H3K4me3'], 4000)
    for s in pls_seqs:
        s['type'] = 'natural_wt'
    library.extend(pls_seqs)
    
    # 1.2 Proximal Enhancers (pELS): 4,000 sequences
    pels_seqs = sample_natural_ccres(ccres_data, chrom_seqs, ['Proximal enhancer'], 4000)
    for s in pels_seqs:
        s['type'] = 'natural_wt'
    library.extend(pels_seqs)
    
    # 1.3 Distal Enhancers (dELS): 4,000 sequences
    dels_seqs = sample_natural_ccres(ccres_data, chrom_seqs, ['Distal enhancer'], 4000)
    for s in dels_seqs:
        s['type'] = 'natural_wt'
    library.extend(dels_seqs)
    
    # 1.4 CTCF-only / Insulators (CTCF): 2,000 sequences
    ctcf_seqs = sample_natural_ccres(ccres_data, chrom_seqs, ['CA-CTCF'], 2000)
    for s in ctcf_seqs:
        s['type'] = 'natural_wt'
    library.extend(ctcf_seqs)
    
    # 1.5 Neutral Genomic backgrounds: 2,000 sequences
    neutral_seqs = sample_neutral_genomic(ccres_data, chrom_seqs, 2000)
    for s in neutral_seqs:
        s['type'] = 'natural_wt'
    library.extend(neutral_seqs)
    
    # Verify Natural total
    natural_wt_count = len(pls_seqs) + len(pels_seqs) + len(dels_seqs) + len(ctcf_seqs) + len(neutral_seqs)
    print(f"Total Natural Genomic WT sequences: {natural_wt_count}")
    
    # -------------------------------------------------------------
    # PART 2: Synthetic Motif Grammar (16,000 sequences)
    # -------------------------------------------------------------
    print("\n--- Generating PART 2: Synthetic Motif Grammar ---")
    
    # 2.1 Scanning: 5,333 sequences
    scanning_seqs = create_scanning_seqs(5333)
    for s in scanning_seqs:
        s['type'] = 'synthetic_wt'
    library.extend(scanning_seqs)
    
    # 2.2 Homotypic: 5,333 sequences
    homotypic_seqs = create_homotypic_seqs(5333)
    for s in homotypic_seqs:
        s['type'] = 'synthetic_wt'
    library.extend(homotypic_seqs)
    
    # 2.3 Heterotypic: 5,334 sequences
    heterotypic_seqs = create_heterotypic_seqs(5334)
    for s in heterotypic_seqs:
        s['type'] = 'synthetic_wt'
    library.extend(heterotypic_seqs)
    
    synthetic_wt_count = len(scanning_seqs) + len(homotypic_seqs) + len(heterotypic_seqs)
    print(f"Total Synthetic Grammar WT sequences: {synthetic_wt_count}")
    
    # -------------------------------------------------------------
    # PART 3: Pure Random Backgrounds (3,000 sequences)
    # -------------------------------------------------------------
    print("\n--- Generating PART 3: Pure Random Backgrounds ---")
    random_seqs = []
    for i in range(3000):
        # Vary GC content continuously from 30% to 75%
        gc = 0.30 + 0.45 * (i / 3000.0)
        seq = generate_random_bg(200, gc)
        random_seqs.append({
            'sequence': seq,
            'type': 'random_bg',
            'class': 'Random',
            'details': f"gc={gc:.4f}"
        })
    library.extend(random_seqs)
    print(f"Generated 3000 pure random backgrounds.")
    
    # -------------------------------------------------------------
    # PART 4: Mutants / Perturbation Pairs (15,000 sequences)
    # -------------------------------------------------------------
    print("\n--- Generating PART 4: Mutants / Perturbation Pairs ---")
    
    # Assign unique IDs to the current library entries so we can link mutants to their parents
    for idx, item in enumerate(library):
        item['id'] = f"SEQ_{idx+1:05d}"
        item['parent_id'] = "None"
        
    # We want exactly 7,500 mutants of natural WT sequences, and 7,500 mutants of synthetic WT sequences
    natural_wt_indices = [i for i, item in enumerate(library) if item['type'] == 'natural_wt']
    synthetic_wt_indices = [i for i, item in enumerate(library) if item['type'] == 'synthetic_wt']
    
    # Sample 7,500 of each
    mut_nat_indices = random.sample(natural_wt_indices, 7500)
    mut_syn_indices = random.sample(synthetic_wt_indices, 7500)
    
    mutants = []
    mut_idx = len(library) + 1
    
    print("Mutating selected natural sequences...")
    for idx in mut_nat_indices:
        parent = library[idx]
        mut_seq, mtype, mdetails = mutate_single_seq(parent['sequence'])
        mutants.append({
            'id': f"SEQ_{mut_idx:05d}",
            'sequence': mut_seq,
            'type': 'natural_mut',
            'class': 'Mutant',
            'details': f"parent_class={parent['class']};mtype={mtype};{mdetails}",
            'parent_id': parent['id']
        })
        mut_idx += 1
        
    print("Mutating selected synthetic sequences...")
    for idx in mut_syn_indices:
        parent = library[idx]
        mut_seq, mtype, mdetails = mutate_single_seq(parent['sequence'])
        mutants.append({
            'id': f"SEQ_{mut_idx:05d}",
            'sequence': mut_seq,
            'type': 'synthetic_mut',
            'class': 'Mutant',
            'details': f"parent_class={parent['class']};mtype={mtype};{mdetails}",
            'parent_id': parent['id']
        })
        mut_idx += 1
        
    library.extend(mutants)
    print(f"Generated {len(mutants)} mutational perturbation sequences.")
    
    # -------------------------------------------------------------
    # Verification & Saving
    # -------------------------------------------------------------
    print("\n--- Running Final Library Verifications ---")
    print(f"Total sequences in library: {len(library)}")
    assert len(library) == 50000, f"Library size mismatch! Expected 50000, got {len(library)}"
    
    # Verify that every sequence is exactly 200bp and contains only {A,C,G,T}
    print("Checking sequence lengths and character alphabets...")
    for item in library:
        seq = item['sequence']
        assert len(seq) == 200, f"Sequence {item['id']} is not 200bp! Length is {len(seq)}"
        assert all(x in 'ACGT' for x in seq), f"Sequence {item['id']} contains invalid characters: {set(seq) - set('ACGT')}"
    print("All 50,000 sequences successfully passed validation!")
    
    # Save sequences to sequences.txt
    print(f"Writing sequences to {SEQUENCES_FILE}...")
    with open(SEQUENCES_FILE, 'w') as f:
        for item in library:
            f.write(item['sequence'] + "\n")
            
    # Save metadata to metadata.tsv
    print(f"Writing rich metadata to {METADATA_FILE}...")
    with open(METADATA_FILE, 'w') as f:
        # Write header
        f.write("id\ttype\tclass\tparent_id\tdetails\tsequence\n")
        for item in library:
            f.write(f"{item['id']}\t{item['type']}\t{item['class']}\t{item['parent_id']}\t{item['details']}\t{item['sequence']}\n")
            
    print("\n" + "="*60)
    print("   LIBRARY GENERATION COMPLETED SUCCESSFULLY!")
    print("   Total Sequences: 50,000 lines (200bp)")
    print("   Saved to: library/sequences.txt")
    print("   Metadata: library/metadata.tsv")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
