"""
Proper AAC ID mapping from DGEB to GenBank.
Extract AAC IDs correctly and map to genomic positions.
"""

import gzip
import re
import pandas as pd
import numpy as np
from Bio import SeqIO

print("=" * 80)
print("PROPER AAC ID MAPPING - FINAL FIX")
print("=" * 80)

# Step 1: Parse GenBank and create comprehensive ID mappings
print("\n1. Creating comprehensive protein ID mappings from GenBank...")

protein_id_map = {}  # Maps any ID to genomic info

with gzip.open("data/GCF_000005845.2_ASM584v2_genomic.gbff.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                # Extract genomic info
                loc = feature.location
                genomic_info = {
                    "start": int(loc.start),
                    "end": int(loc.end),
                    "strand": int(loc.strand),
                    "locus_tag": None,
                    "protein_id": None,
                    "gene": None
                }
                
                # Collect all possible IDs
                ids_to_store = []
                
                # protein_id (NP_ format)
                if "protein_id" in feature.qualifiers:
                    pid = feature.qualifiers["protein_id"][0]
                    ids_to_store.append(pid)
                    genomic_info["protein_id"] = pid
                
                # old_protein_id (might have AAC format)
                if "old_protein_id" in feature.qualifiers:
                    for old_id in feature.qualifiers["old_protein_id"]:
                        ids_to_store.append(old_id)
                
                # locus_tag (b#### format)
                if "locus_tag" in feature.qualifiers:
                    locus = feature.qualifiers["locus_tag"][0]
                    ids_to_store.append(locus)
                    genomic_info["locus_tag"] = locus
                
                # gene name
                if "gene" in feature.qualifiers:
                    gene = feature.qualifiers["gene"][0]
                    genomic_info["gene"] = gene
                
                # db_xref might have additional IDs
                if "db_xref" in feature.qualifiers:
                    for xref in feature.qualifiers["db_xref"]:
                        ids_to_store.append(xref)
                
                # Store under all IDs
                for id_val in ids_to_store:
                    if id_val:
                        protein_id_map[id_val] = genomic_info.copy()

print(f"   Created mappings for {len(protein_id_map)} protein IDs")

# Check what ID formats we have
id_formats = {}
for id_val in list(protein_id_map.keys())[:100]:
    if id_val.startswith('NP_'):
        id_formats['NP_'] = id_formats.get('NP_', 0) + 1
    elif id_val.startswith('YP_'):
        id_formats['YP_'] = id_formats.get('YP_', 0) + 1
    elif id_val.startswith('AAC'):
        id_formats['AAC'] = id_formats.get('AAC', 0) + 1
    elif id_val.startswith('b'):
        id_formats['b####'] = id_formats.get('b####', 0) + 1
    elif 'GeneID' in id_val:
        id_formats['GeneID'] = id_formats.get('GeneID', 0) + 1

print(f"   ID formats found (sample): {id_formats}")

# Step 2: Extract AAC IDs from DGEB with proper regex
print("\n2. Extracting IDs from DGEB Entry field...")

def extract_id_from_entry(entry):
    """
    Extract protein ID from DGEB Entry field.
    Format: U00096.3_prot_AAC73112.1_1
            or U00096.3_prot_b0001_1
    """
    if '_prot_' not in entry:
        return None
    
    # Extract the part after _prot_ and before the last _
    match = re.search(r'_prot_([^_]+)_', entry)
    if match:
        return match.group(1)
    
    # Fallback: just take what's between _prot_ and end
    parts = entry.split('_prot_')
    if len(parts) > 1:
        id_part = parts[1].split('_')[0]
        return id_part
    
    return None

df['extracted_id'] = df['Entry'].apply(extract_id_from_entry)

print(f"   Extracted IDs for {df['extracted_id'].notna().sum()}/{len(df)} entries")
print(f"   Missing: {df['extracted_id'].isna().sum()}")

print("\n   Sample extractions:")
for i in range(min(10, len(df))):
    print(f"     {df.iloc[i]['Entry']} → {df.iloc[i]['extracted_id']}")

# Step 3: Map DGEB IDs to genomic positions
print("\n3. Mapping DGEB IDs to genomic positions...")

mapped_count = 0
for idx, row in df.iterrows():
    dgeb_id = row['extracted_id']
    if pd.notna(dgeb_id):
        # Try direct lookup
        if dgeb_id in protein_id_map:
            mapped_count += 1
            continue
        
        # Try adding .1, .2, etc. for version numbers
        for version in ['.1', '.2', '.3']:
            test_id = dgeb_id + version if '.' not in dgeb_id else dgeb_id
            if test_id in protein_id_map:
                mapped_count += 1
                df.loc[idx, 'extracted_id'] = test_id
                break

print(f"   Successfully mapped: {mapped_count}/{df['extracted_id'].notna().sum()} IDs")
print(f"   Mapping success rate: {mapped_count/df['extracted_id'].notna().sum()*100:.1f}%")

# Step 4: Create pairs with proper mapping
print("\n4. Creating pairs with proper ID mapping...")

pairs_final = []
for i in range(len(df) - 1):
    row_a = df.iloc[i]
    row_b = df.iloc[i + 1]
    
    id_a = row_a['extracted_id']
    id_b = row_b['extracted_id']
    
    # Look up genomic info
    genomic_a = protein_id_map.get(id_a)
    genomic_b = protein_id_map.get(id_b)
    
    pair = {
        "Entry_A": row_a["Entry"],
        "Entry_B": row_b["Entry"],
        "dgeb_id_A": id_a,
        "dgeb_id_B": id_b,
        "Sequence_A": row_a["Sequence"],
        "Sequence_B": row_b["Sequence"],
        "Label": row_a["Label"]
    }
    
    if genomic_a:
        pair['protein_id_A'] = genomic_a['protein_id']
        pair['start_A'] = genomic_a['start']
        pair['end_A'] = genomic_a['end']
        pair['strand_A'] = genomic_a['strand']
        pair['locus_tag_A'] = genomic_a['locus_tag']
    else:
        pair['protein_id_A'] = None
        pair['start_A'] = None
        pair['end_A'] = None
        pair['strand_A'] = None
        pair['locus_tag_A'] = None
    
    if genomic_b:
        pair['protein_id_B'] = genomic_b['protein_id']
        pair['start_B'] = genomic_b['start']
        pair['end_B'] = genomic_b['end']
        pair['strand_B'] = genomic_b['strand']
        pair['locus_tag_B'] = genomic_b['locus_tag']
    else:
        pair['protein_id_B'] = None
        pair['start_B'] = None
        pair['end_B'] = None
        pair['strand_B'] = None
        pair['locus_tag_B'] = None
    
    pairs_final.append(pair)

pairs_df_final = pd.DataFrame(pairs_final)

print(f"   Created {len(pairs_df_final)} pairs")

# Calculate distances
def calc_distance(row):
    if pd.notna(row['start_A']) and pd.notna(row['start_B']):
        if row['strand_A'] == row['strand_B']:
            if row['strand_A'] == 1:
                return abs(row['start_B'] - row['end_A'])
            else:
                return abs(row['start_A'] - row['end_B'])
    return None

pairs_df_final['intergenic_distance'] = pairs_df_final.apply(calc_distance, axis=1)
pairs_df_final['same_strand'] = (pairs_df_final['strand_A'] == pairs_df_final['strand_B'])

# Step 5: Compare with buggy version
print("\n" + "=" * 80)
print("5. COMPARISON: PROPER MAPPING vs BUGGY SEQUENCE MAPPING")
print("=" * 80)

# Drop pairs without mappings
initial = len(pairs_df_final)
pairs_df_final = pairs_df_final.dropna(subset=['intergenic_distance']).reset_index(drop=True)
final = len(pairs_df_final)

print(f"\nPairs with valid mappings: {final}/{initial} ({final/initial*100:.1f}%)")

print(f"\nBUGGY (sequence mapping):")
print(f"  N: {len(pairs_df)}")
print(f"  Median distance: {pairs_df['intergenic_distance'].median():,.0f} bp")
print(f"  Max distance: {pairs_df['intergenic_distance'].max():,.0f} bp")
print(f"  >1 Mbp: {(pairs_df['intergenic_distance'] > 1000000).sum()} pairs")

print(f"\nPROPER (AAC ID mapping):")
print(f"  N: {len(pairs_df_final)}")
print(f"  Median distance: {pairs_df_final['intergenic_distance'].median():,.0f} bp")
print(f"  Max distance: {pairs_df_final['intergenic_distance'].max():,.0f} bp")  
print(f"  >1 Mbp: {(pairs_df_final['intergenic_distance'] > 1000000).sum()} pairs")

# Check mega-distance pairs
mega_final = pairs_df_final[pairs_df_final['intergenic_distance'] > 1000000]

print(f"\n" + "=" * 80)
print("6. MEGA-DISTANCE PAIRS (>1 Mbp)")
print("=" * 80)

if len(mega_final) > 0:
    print(f"\n⚠️  Still have {len(mega_final)} mega-distance pairs")
    print("\nTop 5:")
    for i, (_, row) in enumerate(mega_final.nlargest(5, 'intergenic_distance').iterrows()):
        print(f"  {i+1}. {row['protein_id_A']} → {row['protein_id_B']}")
        print(f"     Distance: {row['intergenic_distance']:,.0f} bp")
        print(f"     Label: {row['Label']}")
        print(f"     DGEB IDs: {row['dgeb_id_A']} → {row['dgeb_id_B']}")
    
    print("\n   These might be legitimate DGEB pairs (hard negatives)")
else:
    print("\n✅ SUCCESS! No mega-distance pairs!")
    print("   The duplicate sequence bug is FIXED!")

print("\n" + "=" * 80)
print("✅ PROPER MAPPING COMPLETE")
print("=" * 80)
print(f"\nUse 'pairs_df_final' for analysis!")
print(f"This version properly maps DGEB IDs to genomic positions.")

