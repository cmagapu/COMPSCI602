"""
Debug: Figure out what ID format DGEB uses and what GenBank has
"""

print("=" * 80)
print("DEBUGGING ID FORMAT MISMATCH")
print("=" * 80)

# Step 1: What did we extract from DGEB?
print("\n1. IDs EXTRACTED FROM DGEB")
print("-" * 80)

print("\nFirst 20 extracted IDs:")
for i in range(min(20, len(df))):
    entry = df.iloc[i]['Entry']
    extracted = df.iloc[i]['extracted_id']
    print(f"  {entry}")
    print(f"    → {extracted}")

print(f"\nUnique ID patterns in extracted IDs:")
extracted_ids = df['extracted_id'].dropna()
print(f"  Total extracted: {len(extracted_ids)}")

# Group by pattern
aac_count = sum(1 for id in extracted_ids if str(id).startswith('AAC'))
np_count = sum(1 for id in extracted_ids if str(id).startswith('NP_'))
yp_count = sum(1 for id in extracted_ids if str(id).startswith('YP_'))
b_count = sum(1 for id in extracted_ids if str(id).startswith('b'))
other_count = len(extracted_ids) - aac_count - np_count - yp_count - b_count

print(f"  AAC##### format: {aac_count}")
print(f"  NP_##### format: {np_count}")
print(f"  YP_##### format: {yp_count}")
print(f"  b#### format: {b_count}")
print(f"  Other: {other_count}")

# Step 2: What's actually in GenBank?
print("\n2. IDS IN GENBANK PROTEIN_ID_MAP")
print("-" * 80)

print(f"\nTotal IDs in protein_id_map: {len(protein_id_map)}")
print("\nSample of protein_id_map keys:")
for i, key in enumerate(list(protein_id_map.keys())[:30]):
    print(f"  {key}")
    if i >= 29:
        break

# Check formats in protein_id_map
map_aac = sum(1 for k in protein_id_map.keys() if k.startswith('AAC'))
map_np = sum(1 for k in protein_id_map.keys() if k.startswith('NP_'))
map_yp = sum(1 for k in protein_id_map.keys() if k.startswith('YP_'))
map_b = sum(1 for k in protein_id_map.keys() if k.startswith('b'))

print(f"\nID formats in protein_id_map:")
print(f"  AAC##### format: {map_aac}")
print(f"  NP_##### format: {map_np}")
print(f"  YP_##### format: {map_yp}")
print(f"  b#### format: {map_b}")

# Step 3: Try to match a few manually
print("\n3. MANUAL MATCHING TEST")
print("-" * 80)

print("\nTrying to match first 10 DGEB IDs:")
for i in range(min(10, len(df))):
    dgeb_id = df.iloc[i]['extracted_id']
    entry = df.iloc[i]['Entry']
    
    print(f"\n  DGEB: {entry}")
    print(f"    Extracted ID: {dgeb_id}")
    
    # Try direct match
    if dgeb_id in protein_id_map:
        print(f"    ✓ Direct match found!")
        continue
    
    # Try with version number
    found = False
    for suffix in ['', '.1', '.2', '.3']:
        test_id = str(dgeb_id) + suffix
        if test_id in protein_id_map:
            print(f"    ✓ Match with suffix: {test_id}")
            found = True
            break
    
    if not found:
        # Try to find similar IDs in protein_id_map
        similar = [k for k in protein_id_map.keys() if str(dgeb_id) in str(k)]
        if similar:
            print(f"    ~ Similar IDs in map: {similar[:3]}")
        else:
            print(f"    ✗ No match found!")

# Step 4: Check if we need to extract differently
print("\n4. ALTERNATIVE EXTRACTION PATTERNS")
print("-" * 80)

print("\nTrying different extraction patterns on first entry:")
sample_entry = df.iloc[0]['Entry']
print(f"Sample: {sample_entry}")

import re

# Pattern 1: After _prot_ until _
match1 = re.search(r'_prot_([^_]+)_', sample_entry)
print(f"  Pattern 1 '_prot_([^_]+)_': {match1.group(1) if match1 else 'No match'}")

# Pattern 2: After _prot_ until end or _
match2 = re.search(r'_prot_([A-Z0-9.]+)', sample_entry)
print(f"  Pattern 2 '_prot_([A-Z0-9.]+)': {match2.group(1) if match2 else 'No match'}")

# Pattern 3: Just the AAC part with version
match3 = re.search(r'(AAC[0-9]+\.[0-9]+)', sample_entry)
print(f"  Pattern 3 '(AAC[0-9]+\.[0-9]+)': {match3.group(1) if match3 else 'No match'}")

# Pattern 4: NP_ or YP_ with version
match4 = re.search(r'([NY]P_[0-9]+\.[0-9]+)', sample_entry)
print(f"  Pattern 4 '([NY]P_[0-9]+\.[0-9]+)': {match4.group(1) if match4 else 'No match'}")

print("\n" + "=" * 80)
print("🔍 DIAGNOSIS")
print("=" * 80)
print("""
The issue is likely one of:
1. DGEB IDs are in a format that doesn't exist in GenBank
2. Our extraction pattern is wrong
3. GenBank doesn't have the right cross-references
4. Need to map through sequences after all (but carefully!)
""")
