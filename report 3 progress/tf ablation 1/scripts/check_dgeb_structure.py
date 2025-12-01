"""
Check if DGEB dataset really has consecutive genes in adjacent rows
"""

print("=" * 80)
print("INVESTIGATING DGEB DATASET STRUCTURE")
print("=" * 80)

# Load the original DGEB dataframe
print("\n1. Original DGEB dataset structure:")
print(f"   Total rows: {len(df)}")
print(f"   Columns: {df.columns.tolist()}")
print("\nFirst 10 rows:")
display(df.head(10))

# Check if df is already sorted by genomic position
print("\n2. Checking if DGEB is sorted by genomic position...")

# We need to map sequences to genomic positions
# Let's check a few consecutive entries in df
for i in range(min(10, len(df)-1)):
    entry_a = df.iloc[i]['Entry']
    entry_b = df.iloc[i+1]['Entry']
    
    # Extract protein IDs
    import re
    match_a = re.search(r"U00096\.3_prot_([A-Z0-9]+\.\d+)", entry_a)
    match_b = re.search(r"U00096\.3_prot_([A-Z0-9]+\.\d+)", entry_b)
    
    if match_a and match_b:
        pid_a = match_a.group(1)
        pid_b = match_b.group(1)
        
        if pid_a in gene_annotations and pid_b in gene_annotations:
            ann_a = gene_annotations[pid_a]
            ann_b = gene_annotations[pid_b]
            
            distance = ann_b['start'] - ann_a['end']
            strand_match = "✓" if ann_a['strand'] == ann_b['strand'] else "✗"
            
            print(f"\nRow {i} → {i+1}:")
            print(f"  {pid_a} (pos: {ann_a['start']:,}) → {pid_b} (pos: {ann_b['start']:,})")
            print(f"  Distance: {distance:,} bp")
            print(f"  Same strand: {strand_match}")
            print(f"  Label: {df.iloc[i]['Label']} (1=same operon, 0=different)")

# Now check what create_pairs() actually produces
print("\n" + "=" * 80)
print("3. CHECKING CREATE_PAIRS() OUTPUT")
print("=" * 80)

# Recreate what happened
print("\nOriginal df shape:", df.shape)
print("After extracting protein_id and dropping NaN:", pairs_df.shape)

# Check first few pairs in pairs_df
print("\nFirst 5 pairs in pairs_df:")
for i in range(min(5, len(pairs_df))):
    row = pairs_df.iloc[i]
    refseq_a = row['refseq_A']
    refseq_b = row['refseq_B']
    
    if refseq_a in gene_annotations and refseq_b in gene_annotations:
        ann_a = gene_annotations[refseq_a]
        ann_b = gene_annotations[refseq_b]
        
        distance = ann_b['start'] - ann_a['end']
        
        print(f"\nPair {i}: {refseq_a} → {refseq_b}")
        print(f"  Position A: {ann_a['start']:,} - {ann_a['end']:,}")
        print(f"  Position B: {ann_b['start']:,} - {ann_b['end']:,}")
        print(f"  Distance: {distance:,} bp")
        print(f"  Same strand: {ann_a['strand'] == ann_b['strand']}")
        print(f"  Label: {row['Label']}")

# Check the problematic pairs
print("\n" + "=" * 80)
print("4. CHECKING THE SUSPICIOUS 4.1 Mbp PAIR")
print("=" * 80)

# Find this pair in pairs_df
suspicious = pairs_df[(pairs_df['refseq_A'] == 'NP_414895.1') & 
                      (pairs_df['refseq_B'] == 'NP_418697.1')]

if len(suspicious) > 0:
    print("\n🔍 Found the suspicious pair in pairs_df!")
    print(f"   Index: {suspicious.index[0]}")
    print(f"   Label: {suspicious['Label'].values[0]}")
    print(f"   Intergenic distance: {suspicious['intergenic_distance'].values[0]:,} bp")
    
    # Find where these appear in the ORIGINAL df
    print("\n🔍 Looking for these proteins in original DGEB df...")
    
    # Search for gene A
    mask_a = df['Entry'].str.contains('NP_414895.1', na=False)
    if mask_a.any():
        idx_a = df[mask_a].index[0]
        print(f"   Gene A (NP_414895.1) found at df index: {idx_a}")
    
    # Search for gene B  
    mask_b = df['Entry'].str.contains('NP_418697.1', na=False)
    if mask_b.any():
        idx_b = df[mask_b].index[0]
        print(f"   Gene B (NP_418697.1) found at df index: {idx_b}")
        print(f"   Distance in df indices: {abs(idx_b - idx_a)} rows apart")
        
        if abs(idx_b - idx_a) == 1:
            print("\n   ✓ These ARE adjacent rows in DGEB!")
            print("   → So DGEB does contain non-consecutive (genomically) pairs!")
        else:
            print("\n   ✗ These are NOT adjacent rows in DGEB")
            print("   → The pair was created some other way")
else:
    print("Pair not found in pairs_df")

# Statistical analysis
print("\n" + "=" * 80)
print("5. STATISTICAL ANALYSIS OF ALL PAIRS")
print("=" * 80)

# For all pairs with genomic annotations, calculate distances
distances = []
labels = []
for i, row in pairs_df.iterrows():
    refseq_a = row['refseq_A']
    refseq_b = row['refseq_B']
    
    if refseq_a in gene_annotations and refseq_b in gene_annotations:
        ann_a = gene_annotations[refseq_a]
        ann_b = gene_annotations[refseq_b]
        
        if ann_a['strand'] == ann_b['strand']:
            distance = abs(ann_b['start'] - ann_a['end'])
            distances.append(distance)
            labels.append(row['Label'])

distances = np.array(distances)
labels = np.array(labels)

print(f"\nAll same-strand pairs (N={len(distances)}):")
print(f"  Mean distance: {distances.mean():,.0f} bp")
print(f"  Median distance: {np.median(distances):,.0f} bp")
print(f"  Min: {distances.min():,} bp")
print(f"  Max: {distances.max():,} bp")

# Split by label
operon_pairs = distances[labels == 1]
non_operon_pairs = distances[labels == 0]

print(f"\nSame operon pairs (Label=1, N={len(operon_pairs)}):")
print(f"  Mean distance: {operon_pairs.mean():,.0f} bp")
print(f"  Median distance: {np.median(operon_pairs):,.0f} bp")
print(f"  Max: {operon_pairs.max():,} bp")

print(f"\nDifferent operon pairs (Label=0, N={len(non_operon_pairs)}):")
print(f"  Mean distance: {non_operon_pairs.mean():,.0f} bp")
print(f"  Median distance: {np.median(non_operon_pairs):,.0f} bp")
print(f"  Max: {non_operon_pairs.max():,} bp")

# Count mega-regions
mega = (distances > 1000000).sum()
print(f"\nPairs with >1 Mbp 'intergenic distance': {mega} ({mega/len(distances)*100:.1f}%)")

# Visualize
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

# Distribution by label
ax = axes[0]
ax.hist(np.log10(operon_pairs + 1), bins=50, alpha=0.6, label='Same operon (Label=1)', color='green')
ax.hist(np.log10(non_operon_pairs + 1), bins=50, alpha=0.6, label='Different operons (Label=0)', color='red')
ax.set_xlabel('Log10(Distance + 1)', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title('Distance Distribution by Operon Label\n(DGEB pairs)', fontsize=13, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)
ax.axvline(x=5, color='orange', linestyle='--', linewidth=2, label='100 kb')

# Box plot
ax = axes[1]
bp = ax.boxplot([np.log10(operon_pairs + 1), np.log10(non_operon_pairs + 1)],
                labels=['Same operon\n(Label=1)', 'Different operons\n(Label=0)'],
                patch_artist=True)
bp['boxes'][0].set_facecolor('green')
bp['boxes'][1].set_facecolor('red')
ax.set_ylabel('Log10(Distance + 1)', fontsize=12, fontweight='bold')
ax.set_title('Distance by Operon Label\n(Box Plot)', fontsize=13, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
ax.axhline(y=5, color='orange', linestyle='--', linewidth=2, label='100 kb')

plt.tight_layout()
plt.savefig('graphs/dgeb_distance_analysis.png', dpi=300, bbox_inches='tight')
print("\n💾 Saved to 'graphs/dgeb_distance_analysis.png'")
plt.show()

print("\n" + "=" * 80)
print("🎯 CONCLUSION")
print("=" * 80)
print("""
If DGEB adjacent rows ARE consecutive:
  → Then the dataset includes some pairs that are Mbp apart genomically
  → These might be negative examples (Label=0) or annotation errors
  → The huge "intergenic regions" are real features of the dataset
  
If they're NOT consecutive:
  → Then create_pairs() is creating wrong pairs
  → The whole pairs_df is incorrect
  
Either way, extracting "intergenic regions" between genes that are
Mbp apart makes no biological sense for operon prediction!
""")

