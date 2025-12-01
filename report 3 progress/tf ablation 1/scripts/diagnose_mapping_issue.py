"""
Diagnose the protein sequence → RefSeq ID mapping issue
"""

print("=" * 80)
print("DIAGNOSING SEQUENCE → REFSEQ MAPPING ISSUE")
print("=" * 80)

# Step 1: Check for duplicate sequences in protein_sequences
print("\n1. CHECKING FOR DUPLICATE SEQUENCES IN PROTEIN DATABASE")
print("-" * 80)

# Reload protein sequences
from Bio import SeqIO
protein_sequences = {}
for record in SeqIO.parse(gzip.open('data/GCF_000005845.2_ASM584v2_protein.faa.gz', 'rt'), "fasta"):
    protein_sequences[record.id] = str(record.seq)

print(f"Total proteins in .faa file: {len(protein_sequences)}")

# Check for duplicates
from collections import Counter
seq_counts = Counter(protein_sequences.values())
duplicate_seqs = {seq: count for seq, count in seq_counts.items() if count > 1}

print(f"Unique sequences: {len(seq_counts)}")
print(f"Duplicate sequences: {len(duplicate_seqs)}")

if duplicate_seqs:
    print(f"\n🚨 FOUND {len(duplicate_seqs)} DUPLICATE SEQUENCES!")
    print("Sample duplicates:")
    for i, (seq, count) in enumerate(list(duplicate_seqs.items())[:5]):
        print(f"\n  Duplicate {i+1}: {count} proteins with same sequence")
        print(f"    Sequence length: {len(seq)} aa")
        # Find which proteins have this sequence
        proteins_with_seq = [pid for pid, s in protein_sequences.items() if s == seq]
        print(f"    Proteins: {', '.join(proteins_with_seq[:10])}")
        if len(proteins_with_seq) > 10:
            print(f"    ... and {len(proteins_with_seq) - 10} more")

# Step 2: Check the seq_to_pid mapping
print("\n" + "=" * 80)
print("2. ANALYZING SEQ_TO_PID MAPPING")
print("=" * 80)

seq_to_pid = {seq: pid for pid, seq in protein_sequences.items()}
print(f"\nOriginal proteins: {len(protein_sequences)}")
print(f"Mapped sequences: {len(seq_to_pid)}")
print(f"Lost in mapping: {len(protein_sequences) - len(seq_to_pid)}")

if len(protein_sequences) != len(seq_to_pid):
    print(f"\n🚨 {len(protein_sequences) - len(seq_to_pid)} proteins LOST due to duplicate sequences!")
    print("   → When multiple proteins have same sequence, dict keeps only the last one")
    print("   → This causes incorrect mappings!")

# Step 3: Check mapping for the mega-distance pairs
print("\n" + "=" * 80)
print("3. CHECKING MEGA-DISTANCE PAIR SEQUENCES")
print("=" * 80)

# Get worst pairs
worst_pairs = pairs_df.nlargest(3, 'intergenic_distance')

for i, (idx, row) in enumerate(worst_pairs.iterrows()):
    print(f"\n{i+1}. Pair: {row['refseq_A']} → {row['refseq_B']}")
    print(f"   Distance: {row['intergenic_distance']:,} bp")
    
    # Check sequences
    seq_a = row['Sequence_A']
    seq_b = row['Sequence_B']
    
    print(f"   Sequence A length: {len(seq_a)} aa")
    print(f"   Sequence B length: {len(seq_b)} aa")
    
    # Check if these sequences map correctly
    mapped_a = seq_to_pid.get(seq_a, None)
    mapped_b = seq_to_pid.get(seq_b, None)
    
    print(f"   Mapped A: {mapped_a} (matches refseq_A: {mapped_a == row['refseq_A']})")
    print(f"   Mapped B: {mapped_b} (matches refseq_B: {mapped_b == row['refseq_B']})")
    
    # Check if sequences are duplicates
    if seq_a in duplicate_seqs:
        print(f"   ⚠️  Sequence A is a DUPLICATE! ({duplicate_seqs[seq_a]} proteins share it)")
    if seq_b in duplicate_seqs:
        print(f"   ⚠️  Sequence B is a DUPLICATE! ({duplicate_seqs[seq_b]} proteins share it)")

# Step 4: Trace back to original DGEB for these pairs
print("\n" + "=" * 80)
print("4. TRACING BACK TO ORIGINAL DGEB")
print("=" * 80)

for i, (idx, row) in enumerate(worst_pairs.iterrows()):
    print(f"\n{i+1}. Pair index {idx} in pairs_df:")
    print(f"   refseq_A: {row['refseq_A']}, refseq_B: {row['refseq_B']}")
    
    # These pairs were created from consecutive rows in df
    # pairs_df row i comes from df rows i and i+1
    # But we need to account for dropped rows...
    
    # Try to find by sequence
    seq_a = row['Sequence_A']
    seq_b = row['Sequence_B']
    
    mask_a = df['Sequence'] == seq_a
    mask_b = df['Sequence'] == seq_b
    
    if mask_a.any():
        idx_a_in_df = df[mask_a].index[0]
        print(f"   Sequence A found at df index: {idx_a_in_df}")
        print(f"     Entry: {df.iloc[idx_a_in_df]['Entry']}")
        print(f"     protein_id: {df.iloc[idx_a_in_df]['protein_id']}")
    else:
        print(f"   Sequence A not found in df")
    
    if mask_b.any():
        idx_b_in_df = df[mask_b].index[0]
        print(f"   Sequence B found at df index: {idx_b_in_df}")
        print(f"     Entry: {df.iloc[idx_b_in_df]['Entry']}")
        print(f"     protein_id: {df.iloc[idx_b_in_df]['protein_id']}")
        
        if mask_a.any():
            print(f"   Rows apart in df: {abs(idx_b_in_df - idx_a_in_df)}")
    else:
        print(f"   Sequence B not found in df")

print("\n" + "=" * 80)
print("🎯 ROOT CAUSE ANALYSIS")
print("=" * 80)

if duplicate_seqs:
    print(f"""
ROOT CAUSE: DUPLICATE SEQUENCES IN PROTEIN DATABASE

✓ {len(duplicate_seqs)} sequences are shared by multiple proteins
✓ seq_to_pid dict keeps only ONE protein per sequence (the last one)
✓ This causes INCORRECT MAPPINGS:
  - Sequence from gene X in DGEB
  - Maps to gene Y's RefSeq ID (because they share sequence)
  - Creates artificial "pairs" of distant genes

EXAMPLE:
  - df row i has gene with sequence S
  - df row i+1 has gene with sequence T
  - create_pairs() makes pair (S, T) - looks correct
  - map_protein_ids_to_pairs() maps:
    - S → NP_X (but S might actually belong to NP_A, NP_B, NP_X)
    - T → NP_Y (but T might actually belong to NP_C, NP_D, NP_Y)
  - If wrong mappings chosen → distant genes paired!

BIOLOGICAL NOTE:
Duplicate protein sequences are common in bacteria due to:
- Gene duplications
- Paralogs (related genes from duplication events)
- Conserved domains
- Short proteins with limited sequence space

SOLUTION NEEDED:
- Use UNIQUE identifiers instead of sequences for mapping
- Or: Track one-to-many sequence→protein relationships
- Or: Use the AAC IDs from DGEB directly with GenBank annotations
""")
else:
    print("""
ROOT CAUSE: UNKNOWN
- No duplicate sequences found
- Need more investigation
""")
