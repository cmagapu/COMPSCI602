"""
Let's systematically debug the methodology to see if the mega-distance pairs
are due to DGEB data quality issues or bugs in our processing pipeline.
"""

print("=" * 80)
print("DEBUGGING METHODOLOGY: DGEB vs OUR PIPELINE")
print("=" * 80)

# Step 1: Check if the mega-distance pairs exist in the ORIGINAL DGEB df
print("\n1. CHECKING ORIGINAL DGEB DATASET")
print("-" * 80)

# Find the worst offenders in pairs_df
worst_pairs = pairs_df.nlargest(5, 'intergenic_distance')[['refseq_A', 'refseq_B', 'intergenic_distance', 'Label']]
print("Top 5 worst distance pairs in pairs_df:")
for i, (_, row) in enumerate(worst_pairs.iterrows()):
    print(f"\n{i+1}. {row['refseq_A']} → {row['refseq_B']}")
    print(f"   Distance: {row['intergenic_distance']:,} bp ({row['intergenic_distance']/1e6:.2f} Mbp)")
    print(f"   Label: {row['Label']} ({'same operon' if row['Label'] == 1 else 'different operons'})")

# Now check if these pairs are adjacent in the ORIGINAL DGEB df
print("\n" + "="*80)
print("2. CHECKING IF WORST PAIRS ARE ADJACENT IN ORIGINAL DGEB")
print("="*80)

for i, (_, row) in enumerate(worst_pairs.iterrows()):
    gene_a = row['refseq_A']
    gene_b = row['refseq_B']
    
    print(f"\n{i+1}. Checking {gene_a} → {gene_b} ({row['intergenic_distance']:,} bp)")
    
    # Search in original df
    mask_a = df['Entry'].str.contains(gene_a, na=False)
    mask_b = df['Entry'].str.contains(gene_b, na=False)
    
    if mask_a.any() and mask_b.any():
        idx_a = df[mask_a].index[0]
        idx_b = df[mask_b].index[0]
        
        print(f"   Gene A at df row: {idx_a}")
        print(f"   Gene B at df row: {idx_b}")
        print(f"   Rows apart: {abs(idx_b - idx_a)}")
        
        if abs(idx_b - idx_a) == 1:
            print("   ✓ ADJACENT in DGEB - This is a DGEB data quality issue!")
            # Check the label in DGEB
            label_in_dgeb = df.iloc[min(idx_a, idx_b)]['Label']
            print(f"   DGEB Label: {label_in_dgeb} (should be 0 for such large distance)")
            if label_in_dgeb == 1 and row['intergenic_distance'] > 1000000:
                print("   🚨 CONFIRMED: DGEB has Label=1 for 1+ Mbp distance - this is wrong!")
        else:
            print("   ✗ NOT adjacent in DGEB - This is a pipeline bug!")
    else:
        print("   ❌ Genes not found in original DGEB")

# Step 3: Check the create_pairs() function logic
print("\n" + "="*80)
print("3. CHECKING CREATE_PAIRS() LOGIC")
print("="*80)

print("\nLet's trace through create_pairs() for a few examples...")

# Look at the first few rows of original df
print("\nFirst 10 rows of original DGEB df:")
for i in range(min(10, len(df))):
    entry = df.iloc[i]['Entry']
    label = df.iloc[i]['Label']
    print(f"  Row {i}: {entry} (Label: {label})")

print("\nFirst 5 pairs created by create_pairs():")
for i in range(min(5, len(pairs_df))):
    row = pairs_df.iloc[i]
    print(f"  Pair {i}: {row['refseq_A']} → {row['refseq_B']} (Label: {row['Label']})")

# Step 4: Check if there's a systematic issue
print("\n" + "="*80)
print("4. SYSTEMATIC ANALYSIS")
print("="*80)

# Count how many mega-distance pairs are adjacent in DGEB vs not
mega_pairs = pairs_df[pairs_df['intergenic_distance'] > 1000000]
print(f"\nTotal pairs with >1 Mbp distance: {len(mega_pairs)}")

adjacent_in_dgeb = 0
not_adjacent = 0
not_found = 0

for _, row in mega_pairs.iterrows():
    gene_a = row['refseq_A']
    gene_b = row['refseq_B']
    
    mask_a = df['Entry'].str.contains(gene_a, na=False)
    mask_b = df['Entry'].str.contains(gene_b, na=False)
    
    if mask_a.any() and mask_b.any():
        idx_a = df[mask_a].index[0]
        idx_b = df[mask_b].index[0]
        
        if abs(idx_b - idx_a) == 1:
            adjacent_in_dgeb += 1
        else:
            not_adjacent += 1
    else:
        not_found += 1

print(f"\nBreakdown of mega-distance pairs:")
print(f"  Adjacent in DGEB: {adjacent_in_dgeb} ({adjacent_in_dgeb/len(mega_pairs)*100:.1f}%)")
print(f"  Not adjacent: {not_adjacent} ({not_adjacent/len(mega_pairs)*100:.1f}%)")
print(f"  Not found: {not_found} ({not_found/len(mega_pairs)*100:.1f}%)")

if adjacent_in_dgeb > 0:
    print(f"\n🚨 {adjacent_in_dgeb} mega-distance pairs ARE adjacent in DGEB!")
    print("   → This suggests DGEB data quality issues")
    
    # Check labels of these adjacent mega-pairs
    mega_adjacent_labels = []
    for _, row in mega_pairs.iterrows():
        gene_a = row['refseq_A']
        gene_b = row['refseq_B']
        
        mask_a = df['Entry'].str.contains(gene_a, na=False)
        mask_b = df['Entry'].str.contains(gene_b, na=False)
        
        if mask_a.any() and mask_b.any():
            idx_a = df[mask_a].index[0]
            idx_b = df[mask_b].index[0]
            
            if abs(idx_b - idx_a) == 1:
                label = df.iloc[min(idx_a, idx_b)]['Label']
                mega_adjacent_labels.append(label)
    
    if mega_adjacent_labels:
        print(f"   Labels of adjacent mega-pairs: {mega_adjacent_labels}")
        label_1_count = sum(1 for l in mega_adjacent_labels if l == 1)
        print(f"   {label_1_count}/{len(mega_adjacent_labels)} have Label=1 (same operon)")
        if label_1_count > 0:
            print("   🚨 CONFIRMED: DGEB has Label=1 for mega-distance pairs!")

if not_adjacent > 0:
    print(f"\n⚠️  {not_adjacent} mega-distance pairs are NOT adjacent in DGEB")
    print("   → This suggests a bug in our create_pairs() or mapping logic")

# Step 5: Check the protein_id extraction
print("\n" + "="*80)
print("5. CHECKING PROTEIN_ID EXTRACTION")
print("="*80)

# Check if protein_id extraction is working correctly
print("\nSample protein_id extractions:")
for i in range(min(5, len(df))):
    entry = df.iloc[i]['Entry']
    protein_id = extract_protein_id(entry)
    print(f"  {entry} → {protein_id}")

# Check if any protein_ids are None
none_count = df['protein_id'].isna().sum()
print(f"\nEntries with None protein_id: {none_count}")

if none_count > 0:
    print("Sample entries with None protein_id:")
    none_entries = df[df['protein_id'].isna()]['Entry'].head()
    for entry in none_entries:
        print(f"  {entry}")

print("\n" + "="*80)
print("🎯 DIAGNOSIS")
print("="*80)

if adjacent_in_dgeb > 0:
    print("""
DIAGNOSIS: DGEB DATA QUALITY ISSUES
✓ Mega-distance pairs ARE adjacent in DGEB
✓ Some have Label=1 (same operon) despite 1+ Mbp distance
→ This is biologically impossible
→ DGEB has annotation errors or data quality issues
→ Our methodology is correct, but DGEB is wrong
""")
elif not_adjacent > 0:
    print("""
DIAGNOSIS: OUR METHODOLOGY BUG
✗ Mega-distance pairs are NOT adjacent in DGEB
→ Our create_pairs() or mapping logic has a bug
→ We're creating wrong pairs somehow
→ Need to fix our pipeline
""")
else:
    print("""
DIAGNOSIS: UNCLEAR
→ Need more investigation
→ Check if mega-pairs exist in original DGEB at all
""")
