import os
import random
import gzip
import pandas as pd
import numpy as np
from tqdm import tqdm
import twobitreader
from dinuc_shuffle import dinuc_shuffle

def main():
    print("Initializing hg38 genome from 2bit file...")
    genome_path = "data/hg38.2bit"
    if not os.path.exists(genome_path):
        raise FileNotFoundError(f"Genome file not found at {genome_path}")
    genome = twobitreader.TwoBitFile(genome_path)
    
    dhs_path = "data/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz"
    print(f"Loading DHS Index from {dhs_path}...")
    
    # Load DHS Index. Load only required columns to be efficient and avoid parsing issues with other columns.
    df = pd.read_csv(dhs_path, sep="\t", usecols=["seqname", "mean_signal", "summit", "component"], dtype={
        "seqname": str,
        "mean_signal": float,
        "summit": float,
        "component": str
    })
    
    # Drop rows with NA in our required columns and convert summit to integer
    df = df.dropna(subset=["seqname", "mean_signal", "summit", "component"])
    df["summit"] = df["summit"].astype(int)
    
    print(f"Total rows loaded: {len(df)}")
    
    # List of components
    components = sorted(df["component"].unique())
    print(f"Found {len(components)} distinct DHS components:")
    for i, comp in enumerate(components):
        print(f"  {i+1}. {comp}")
        
    # We want 35,000 biological sequences stratified across the 16 components.
    # 8 components will have 2,188 sequences, 8 components will have 2,187 sequences.
    num_seqs_per_component = [2188] * 8 + [2187] * 8
    # Shuffle this count list to be fair across components
    random.seed(42)
    random.shuffle(num_seqs_per_component)
    
    component_targets = dict(zip(components, num_seqs_per_component))
    
    # We will also shuffle the top 625 sequences from each component to make 10,000 shuffled sequences total
    shuf_target_per_component = 625
    
    biological_sequences = []
    shuffled_sequences = []
    
    print("\nExtracting biological and shuffled sequences...")
    for comp in tqdm(components):
        target_biol = component_targets[comp]
        
        # Filter for this component and sort by mean_signal descending
        comp_df = df[df["component"] == comp].sort_values(by="mean_signal", ascending=False)
        
        biol_extracted = []
        shuf_extracted = []
        
        # Iterate over sorted candidates to find valid ones (no Ns)
        for _, row in comp_df.iterrows():
            chrom = row["seqname"]
            summit = int(row["summit"])
            
            if chrom not in genome:
                continue
                
            start_pos = summit - 100
            end_pos = summit + 100
            
            # Check boundaries
            if start_pos < 0 or end_pos > len(genome[chrom]):
                continue
                
            seq = genome[chrom][start_pos:end_pos].upper()
            
            # Check for non-ACGT bases
            if any(char not in "ACGT" for char in seq):
                continue
                
            # Valid sequence found
            biol_extracted.append(seq)
            
            # If we haven't reached the shuffle target yet, generate dinuc shuffle
            if len(shuf_extracted) < shuf_target_per_component:
                shuf_seq = dinuc_shuffle(seq)
                shuf_extracted.append(shuf_seq)
                
            if len(biol_extracted) >= target_biol:
                break
                
        if len(biol_extracted) < target_biol:
            print(f"WARNING: Only extracted {len(biol_extracted)}/{target_biol} for component {comp}")
            
        biological_sequences.extend(biol_extracted)
        shuffled_sequences.extend(shuf_extracted)
        
    print(f"\nExtracted {len(biological_sequences)} biological sequences.")
    print(f"Extracted {len(shuffled_sequences)} dinucleotide-shuffled sequences.")
    
    # Generate 5,000 fully random synthetic sequences of length 200bp
    print("Generating 5,000 fully random synthetic sequences...")
    synthetic_sequences = []
    # Use standard seed for reproducible synthetic generation
    np.random.seed(42)
    bases = ["A", "C", "G", "T"]
    for _ in range(5000):
        seq = "".join(np.random.choice(bases, size=200))
        synthetic_sequences.append(seq)
        
    # Combine all sequences
    all_sequences = biological_sequences + shuffled_sequences + synthetic_sequences
    print(f"Combined total sequences: {len(all_sequences)}")
    
    # Verify exact constraints
    assert len(all_sequences) == 50000, f"Error: Total sequences is {len(all_sequences)}, expected 50000"
    for i, seq in enumerate(all_sequences):
        assert len(seq) == 200, f"Error: Sequence at index {i} has length {len(seq)}, expected 200"
        assert all(c in "ACGT" for c in seq), f"Error: Sequence at index {i} contains non-ACGT bases: {seq}"
        
    print("All sequences passed validation constraints!")
    
    # Shuffle the final list to remove sequential structure biases
    print("Randomly shuffling the combined library...")
    random.shuffle(all_sequences)
    
    # Create library/ directory and write sequences
    os.makedirs("library", exist_ok=True)
    out_path = "library/sequences.txt"
    print(f"Writing sequences to {out_path}...")
    with open(out_path, "w") as f_out:
        for seq in all_sequences:
            f_out.write(seq + "\n")
            
    print("Library successfully generated!")

if __name__ == "__main__":
    main()
