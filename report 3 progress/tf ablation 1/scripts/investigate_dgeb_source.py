"""
Investigate DGEB's data source and structure.
If DGEB got data from BioCyc (genome-ordered consecutive genes),
then mega-distance pairs that are NOT adjacent in DGEB are highly suspicious.
"""

print("=" * 80)
print("INVESTIGATING DGEB DATA SOURCE AND STRUCTURE")
print("=" * 80)

# Check the DGEB metadata
print("\n1. DGEB METADATA")
print("-" * 80)
try:
    if hasattr(task, 'metadata'):
        print(f"Task metadata available")
        # Try to print available attributes
        try:
            print(f"  Type: {task.metadata.type if hasattr(task.metadata, 'type') else 'Unknown'}")
        except:
            pass
        try:
            print(f"  Description: {task.metadata.description if hasattr(task.metadata, 'description') else 'Unknown'}")
        except:
            pass
    print(f"\nDataset info:")
    print(f"  Splits: {hf_ds.keys() if hf_ds else 'Unknown'}")
    print(f"  Total examples: {len(df)}")
except Exception as e:
    print(f"  Could not access metadata: {e}")

# Step 2: Check if DGEB df is in genomic order
print("\n" + "=" * 80)
print("2. IS DGEB ORDERED BY GENOMIC POSITION?")
print("=" * 80)

print("\nChecking first 20 consecutive entries...")

# For each consecutive pair in df, check genomic distance
genomic_ordered = True
non_consecutive_count = 0

for i in range(min(20, len(df)-1)):
    entry_a = df.iloc[i]['Entry']
    entry_b = df.iloc[i+1]['Entry']
    
    # Try to extract protein IDs
    import re
    match_a = re.search(r"U00096\.3_prot_([A-Z0-9]+\.\d+)", entry_a)
    match_b = re.search(r"U00096\.3_prot_([A-Z0-9]+\.\d+)", entry_b)
    
    if match_a and match_b:
        # These are protein IDs from DGEB
        pid_a_dgeb = match_a.group(1)
        pid_b_dgeb = match_b.group(1)
        
        # Now try to map to genomic positions
        # We need to use our original gene_annotations, but the DGEB IDs might be AAC format
        # Let's check genomic order
        
        # Try to find in our protein sequences by matching sequences
        seq_a = df.iloc[i]['Sequence']
        seq_b = df.iloc[i+1]['Sequence']
        
        # Find which proteins have these sequences
        matching_a = [pid for pid, seq in protein_sequences.items() if seq == seq_a]
        matching_b = [pid for pid, seq in protein_sequences.items() if seq == seq_b]
        
        if len(matching_a) == 1 and len(matching_b) == 1:
            # Unique matches - check if they're consecutive in genome
            pid_a = matching_a[0]
            pid_b = matching_b[0]
            
            if pid_a in gene_annotations and pid_b in gene_annotations:
                ann_a = gene_annotations[pid_a]
                ann_b = gene_annotations[pid_b]
                
                distance = abs(ann_b['start'] - ann_a['end'])
                
                print(f"\nRow {i}→{i+1}: {pid_a} → {pid_b}")
                print(f"  Distance: {distance:,} bp")
                print(f"  Label: {df.iloc[i]['Label']}")
                
                if distance > 10000:
                    print(f"  ⚠️  Large distance for consecutive DGEB rows!")
                    genomic_ordered = False
                    non_consecutive_count += 1

if genomic_ordered and non_consecutive_count == 0:
    print("\n✅ DGEB appears to be in genomic order (first 20 pairs)")
else:
    print(f"\n⚠️  DGEB has {non_consecutive_count} non-consecutive pairs in first 20")

# Step 3: Look at one of the problematic mega-pairs in detail
print("\n" + "=" * 80)
print("3. DETAILED ANALYSIS OF ONE MEGA-PAIR")
print("=" * 80)

# Get one mega-pair from pairs_df_fixed
mega_pair = pairs_df_fixed[pairs_df_fixed['intergenic_distance_fixed'] > 1000000].iloc[0]

print(f"\nMega-pair: {mega_pair['protein_id_A']} → {mega_pair['protein_id_B']}")
print(f"Distance: {mega_pair['intergenic_distance_fixed']:,.0f} bp")
print(f"Label: {mega_pair['Label']}")

# Find these sequences in DGEB
seq_a = mega_pair['Sequence_A']
seq_b = mega_pair['Sequence_B']

mask_a = df['Sequence'] == seq_a
mask_b = df['Sequence'] == seq_b

if mask_a.any() and mask_b.any():
    idx_a_dgeb = df[mask_a].index[0]
    idx_b_dgeb = df[mask_b].index[0]
    
    print(f"\nIn DGEB:")
    print(f"  Sequence A at row: {idx_a_dgeb}")
    print(f"    Entry: {df.iloc[idx_a_dgeb]['Entry']}")
    print(f"  Sequence B at row: {idx_b_dgeb}")
    print(f"    Entry: {df.iloc[idx_b_dgeb]['Entry']}")
    print(f"  Rows apart: {abs(idx_b_dgeb - idx_a_dgeb)}")
    
    # Now check: who are sequence A's TRUE genomic neighbors in DGEB?
    print(f"\n  Sequence A's DGEB neighbors:")
    if idx_a_dgeb > 0:
        prev_entry = df.iloc[idx_a_dgeb - 1]['Entry']
        prev_seq = df.iloc[idx_a_dgeb - 1]['Sequence']
        print(f"    Previous (row {idx_a_dgeb-1}): {prev_entry}")
    if idx_a_dgeb < len(df) - 1:
        next_entry = df.iloc[idx_a_dgeb + 1]['Entry']
        next_seq = df.iloc[idx_a_dgeb + 1]['Sequence']
        print(f"    Next (row {idx_a_dgeb+1}): {next_entry}")
    
    # Check if sequence A appears multiple times in DGEB
    all_matches_a = df[df['Sequence'] == seq_a]
    if len(all_matches_a) > 1:
        print(f"\n  ⚠️  Sequence A appears {len(all_matches_a)} times in DGEB!")
        print("  Rows:", all_matches_a.index.tolist())
        print("  → This is why our pairs get misaligned!")

# Step 4: Check for duplicate sequences in DGEB
print("\n" + "=" * 80)
print("4. DUPLICATE SEQUENCES IN DGEB")
print("=" * 80)

from collections import Counter
dgeb_seq_counts = Counter(df['Sequence'])
dgeb_duplicates = {seq: count for seq, count in dgeb_seq_counts.items() if count > 1}

print(f"\nUnique sequences in DGEB: {len(dgeb_seq_counts)}")
print(f"Duplicate sequences in DGEB: {len(dgeb_duplicates)}")
print(f"Total proteins in DGEB: {len(df)}")

if dgeb_duplicates:
    print(f"\nSample duplicates in DGEB:")
    for i, (seq, count) in enumerate(list(dgeb_duplicates.items())[:5]):
        print(f"\n  Duplicate {i+1}: {count} entries with same sequence")
        print(f"    Sequence length: {len(seq)} aa")
        # Find which entries have this sequence
        entries_with_seq = df[df['Sequence'] == seq]['Entry'].tolist()
        print(f"    Entries: {', '.join(entries_with_seq[:10])}")

print("\n" + "=" * 80)
print("💡 KEY INSIGHT")
print("=" * 80)

if dgeb_duplicates:
    print("""
DGEB HAS DUPLICATE SEQUENCES!

This means:
1. DGEB includes paralogous genes (same sequence, different locations)
2. When we use create_pairs(), we pair consecutive DGEB rows correctly
3. But when we map sequences → RefSeq IDs, we can pick the WRONG paralog
4. This creates artificial mega-distance pairs

EXAMPLE:
  DGEB row 100: Gene X (sequence S, at position 100kb)
  DGEB row 101: Gene Y (sequence T, at position 101kb)  ← correct pair
  
  But sequence S also exists at position 3Mbp (Gene X')
  
  When we map sequence S → RefSeq:
  - If we pick Gene X' instead of Gene X
  - We calculate distance between Gene X' (3Mbp) and Gene Y (101kb)
  - Result: 2.9 Mbp "intergenic distance" (wrong!)

THE PROBLEM:
  Our mapping method (sequence → RefSeq ID) cannot distinguish between paralogs!
  We need position-aware mapping or direct ID matching.
""")
else:
    print("""
DGEB doesn't have duplicate sequences (surprising!).
Need to investigate further why mega-pairs are not adjacent.
""")
