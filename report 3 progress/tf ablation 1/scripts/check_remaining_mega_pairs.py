"""
Check if the remaining 11 mega-distance pairs are adjacent in DGEB
or if they're artifacts that slipped through the fix.
"""

print("=" * 80)
print("INVESTIGATING REMAINING 11 MEGA-DISTANCE PAIRS")
print("=" * 80)

# Get the remaining mega pairs
mega_fixed = pairs_df_fixed[pairs_df_fixed['intergenic_distance_fixed'] > 1000000].copy()

print(f"\nFound {len(mega_fixed)} pairs with >1 Mbp distance")
print("\nChecking if these are adjacent in DGEB...")

adjacent_count = 0
not_adjacent_count = 0
not_found = 0

for i, (idx, row) in enumerate(mega_fixed.iterrows()):
    print(f"\n{i+1}. {row['protein_id_A']} → {row['protein_id_B']}")
    print(f"   Distance: {row['intergenic_distance_fixed']:,.0f} bp")
    print(f"   Label: {row['Label']}")
    
    # Find in original DGEB by sequence
    seq_a = row['Sequence_A']
    seq_b = row['Sequence_B']
    
    mask_a = df['Sequence'] == seq_a
    mask_b = df['Sequence'] == seq_b
    
    if mask_a.any() and mask_b.any():
        idx_a = df[mask_a].index[0]
        idx_b = df[mask_b].index[0]
        
        print(f"   DGEB row A: {idx_a}")
        print(f"   DGEB row B: {idx_b}")
        print(f"   Rows apart: {abs(idx_b - idx_a)}")
        
        if abs(idx_b - idx_a) == 1:
            print(f"   ✓ ADJACENT in DGEB - This is a DGEB issue!")
            adjacent_count += 1
        else:
            print(f"   ✗ NOT adjacent ({abs(idx_b - idx_a)} rows apart) - Still a mapping bug!")
            not_adjacent_count += 1
    else:
        print(f"   ❌ Not found in DGEB")
        not_found += 1

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"\nRemaining mega-distance pairs: {len(mega_fixed)}")
print(f"  Adjacent in DGEB: {adjacent_count} ({adjacent_count/len(mega_fixed)*100:.0f}%)")
print(f"  Not adjacent: {not_adjacent_count} ({not_adjacent_count/len(mega_fixed)*100:.0f}%)")
print(f"  Not found: {not_found}")

if adjacent_count == len(mega_fixed):
    print("\n✅ ALL remaining mega-pairs ARE adjacent in DGEB!")
    print("   → These are DGEB data quality issues, not our bug!")
    print("   → DGEB includes some pairs that are Mbp apart (probably as hard negatives)")
elif not_adjacent_count > 0:
    print(f"\n⚠️  {not_adjacent_count} mega-pairs are NOT adjacent in DGEB")
    print("   → The duplicate sequence bug is still partially present")
    print("   → Need better mapping strategy")

# Check labels
print("\n" + "=" * 80)
print("LABEL DISTRIBUTION OF MEGA-PAIRS")
print("=" * 80)

label_dist = mega_fixed['Label'].value_counts()
print(f"\nSame operon (Label=1): {label_dist.get(1, 0)}")
print(f"Different operons (Label=0): {label_dist.get(0, 0)}")

if label_dist.get(1, 0) > 0:
    print(f"\n🚨 {label_dist.get(1, 0)} mega-pairs have Label=1 (same operon)")
    print("   → This is biologically impossible!")
    print("   → Definite DGEB data quality issue")

