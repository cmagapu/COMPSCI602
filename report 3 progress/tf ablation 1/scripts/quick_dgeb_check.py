"""
Quick check: Are the suspicious mega-distance pairs actually 
adjacent rows in the DGEB dataset?
"""

# Check if the 4.1 Mbp pair is adjacent in DGEB
print("🔍 Searching for NP_414895.1 and NP_418697.1 in original DGEB df...\n")

# Find gene A
mask_a = df['Entry'].str.contains('NP_414895.1', na=False)
if mask_a.any():
    idx_a = df[mask_a].index[0]
    print(f"Gene A (NP_414895.1) at df row: {idx_a}")
    print(f"  Entry: {df.iloc[idx_a]['Entry']}")
    print(f"  Label: {df.iloc[idx_a]['Label']}")

# Find gene B
mask_b = df['Entry'].str.contains('NP_418697.1', na=False)
if mask_b.any():
    idx_b = df[mask_b].index[0]
    print(f"\nGene B (NP_418697.1) at df row: {idx_b}")
    print(f"  Entry: {df.iloc[idx_b]['Entry']}")
    if idx_b < len(df) - 1:
        print(f"  Label: {df.iloc[idx_b]['Label']} (label for pair with NEXT gene)")

print(f"\nRows apart in df: {abs(idx_b - idx_a)}")

if abs(idx_b - idx_a) == 1:
    print("\n✓ YES! These ARE adjacent rows in DGEB!")
    print("→ This means DGEB includes genomically non-consecutive pairs")
    print("→ The 4.1 Mbp 'intergenic distance' is a real feature in the dataset")
    print(f"→ Label = {df.iloc[min(idx_a, idx_b)]['Label']} (0=different operons, 1=same operon)")
else:
    print("\n✗ NO! These are NOT adjacent in DGEB")
    print("→ Something else is creating these pairs")

# Check labels of mega-distance pairs
print("\n" + "="*70)
print("Checking labels of pairs with huge distances...")

mega_pairs = pairs_df[pairs_df['intergenic_distance'] > 1000000]
print(f"\nPairs with >1 Mbp distance: {len(mega_pairs)}")
if len(mega_pairs) > 0:
    print(f"  Label distribution:")
    print(f"    Same operon (1): {(mega_pairs['Label'] == 1).sum()}")
    print(f"    Different (0): {(mega_pairs['Label'] == 0).sum()}")

