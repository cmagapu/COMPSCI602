"""
Paste this into a new cell in your notebook to investigate the 
"intergenic region" artifacts and understand why they're clearly wrong.
"""

# Investigate the suspicious "intergenic regions"
print("=" * 80)
print("INVESTIGATING INTERGENIC REGION ARTIFACTS")
print("=" * 80)

# Get E. coli genome size
genome_size = len(genome_record.seq)
print(f"\n📏 E. coli K-12 MG1655 genome size: {genome_size:,} bp ({genome_size/1e6:.2f} Mbp)")

# Check the most suspicious pairs
suspicious_pairs = [
    ("NP_414556.1", "NP_416895.2", 2499032),  # 2.5 Mbp!
    ("NP_414895.1", "NP_418697.1", 4118680),  # 4.1 Mbp! (from FIMO output)
    ("YP_009518741.1", "YP_009518770.2", 1122119),  # 1.1 Mbp
]

print("\n" + "=" * 80)
print("CHECKING SUSPICIOUS 'INTERGENIC REGIONS'")
print("=" * 80)

for gene_a, gene_b, reported_len in suspicious_pairs:
    print(f"\n{'='*70}")
    print(f"Pair: {gene_a} → {gene_b}")
    print(f"Reported 'intergenic' length: {reported_len:,} bp ({reported_len/1e6:.2f} Mbp)")
    print(f"{'='*70}")
    
    if gene_a in gene_annotations and gene_b in gene_annotations:
        ann_a = gene_annotations[gene_a]
        ann_b = gene_annotations[gene_b]
        
        print(f"\n📍 Gene A ({gene_a}):")
        print(f"   Position: {ann_a['start']:,} - {ann_a['end']:,}")
        print(f"   Strand: {'+' if ann_a['strand'] == 1 else '-'}")
        print(f"   Length: {ann_a['end'] - ann_a['start']:,} bp")
        
        print(f"\n📍 Gene B ({gene_b}):")
        print(f"   Position: {ann_b['start']:,} - {ann_b['end']:,}")
        print(f"   Strand: {'+' if ann_b['strand'] == 1 else '-'}")
        print(f"   Length: {ann_b['end'] - ann_b['start']:,} bp")
        
        # Calculate actual distance
        if ann_a['strand'] == 1:
            distance = ann_b['start'] - ann_a['end']
        else:
            distance = ann_a['start'] - ann_b['end']
        
        print(f"\n📏 Calculated intergenic distance: {abs(distance):,} bp ({abs(distance)/1e6:.2f} Mbp)")
        print(f"   → This is {abs(distance)/genome_size*100:.1f}% of the ENTIRE E. coli genome!")
        
        # Count how many genes are in between
        genes_between = []
        for pid, ann in gene_annotations.items():
            if ann['strand'] == ann_a['strand']:  # Same strand
                if ann_a['strand'] == 1:
                    # Forward strand: check if gene is between A and B
                    if ann_a['end'] < ann['start'] < ann_b['start']:
                        genes_between.append(pid)
                else:
                    # Reverse strand: check if gene is between B and A  
                    if ann_b['end'] < ann['start'] < ann_a['start']:
                        genes_between.append(pid)
        
        print(f"\n🧬 Genes on same strand between A and B: {len(genes_between)}")
        if len(genes_between) > 0:
            print(f"   First few: {', '.join(genes_between[:5])}")
        
        print(f"\n🔴 WHY THIS IS WRONG:")
        if abs(distance) > 10000:
            print(f"   • Real intergenic regions are typically 50-500 bp in bacteria")
            print(f"   • This {abs(distance)/1e3:.0f} kb 'region' contains {len(genes_between)} other genes!")
            print(f"   • Genes A and B are NOT consecutive in the genome")
            print(f"   • This is an ARTIFACT of how pairs were created")
    else:
        print("   (Genes not found in annotations)")

# Now check what REAL intergenic regions look like
print("\n" + "=" * 80)
print("WHAT REAL INTERGENIC REGIONS LOOK LIKE")
print("=" * 80)

# Count genes by strand
forward_genes = [g for g, ann in gene_annotations.items() if ann['strand'] == 1]
reverse_genes = [g for g, ann in gene_annotations.items() if ann['strand'] == -1]

print(f"\nTotal genes: {len(gene_annotations):,}")
print(f"  Forward strand (+): {len(forward_genes):,}")
print(f"  Reverse strand (-): {len(reverse_genes):,}")

# Calculate real intergenic distances for consecutive genes on the GENOME
print("\nCalculating real intergenic distances for genomically consecutive genes...")

# Sort genes by position
genes_by_position = sorted(gene_annotations.items(), key=lambda x: x[1]['start'])

real_intergenic_distances = []
for i in range(len(genes_by_position) - 1):
    gene_a, ann_a = genes_by_position[i]
    gene_b, ann_b = genes_by_position[i + 1]
    
    # Only for same strand, non-overlapping
    if ann_a['strand'] == ann_b['strand']:
        distance = ann_b['start'] - ann_a['end']
        if distance > 0:  # Non-overlapping
            real_intergenic_distances.append(distance)

real_intergenic_distances = np.array(real_intergenic_distances)

print(f"\nReal intergenic distances (same strand, consecutive in genome):")
print(f"  N = {len(real_intergenic_distances):,}")
print(f"  Mean: {real_intergenic_distances.mean():.1f} bp")
print(f"  Median: {np.median(real_intergenic_distances):.1f} bp")
print(f"  Min: {real_intergenic_distances.min()} bp")
print(f"  Max: {real_intergenic_distances.max():,} bp")
print(f"  95th percentile: {np.percentile(real_intergenic_distances, 95):.1f} bp")
print(f"  99th percentile: {np.percentile(real_intergenic_distances, 99):.1f} bp")

# Count how many are > 10kb
long_regions = (real_intergenic_distances > 10000).sum()
print(f"\n  Regions > 10 kb: {long_regions} ({long_regions/len(real_intergenic_distances)*100:.1f}%)")

very_long = (real_intergenic_distances > 100000).sum()
print(f"  Regions > 100 kb: {very_long} ({very_long/len(real_intergenic_distances)*100:.2f}%)")

mega_long = (real_intergenic_distances > 1000000).sum()
print(f"  Regions > 1 Mb: {mega_long} ({mega_long/len(real_intergenic_distances)*100:.3f}%)")

# Visualize
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# 1. Distribution of REAL intergenic distances
ax = axes[0]
ax.hist(real_intergenic_distances[real_intergenic_distances < 1000], bins=50, 
        color='green', alpha=0.7, edgecolor='black')
ax.set_xlabel('Intergenic Distance (bp)', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title('Real Intergenic Distances\n(< 1 kb, consecutive in genome)', 
             fontsize=13, fontweight='bold')
ax.axvline(x=real_intergenic_distances.mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Mean: {real_intergenic_distances.mean():.0f} bp')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# 2. Log-scale to see the full distribution
ax = axes[1]
ax.hist(np.log10(real_intergenic_distances + 1), bins=50, 
        color='blue', alpha=0.7, edgecolor='black')
ax.set_xlabel('Log10(Intergenic Distance + 1)', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title('Real Intergenic Distances\n(log scale)', fontsize=13, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# 3. Compare extracted vs real
ax = axes[2]

# Get lengths from intergenic_regions.fasta
with open("data/intergenic_regions.fasta", "r") as f:
    extracted_lengths = []
    for line in f:
        if line.startswith(">"):
            length = int(line.strip().split("len=")[1])
            extracted_lengths.append(length)

extracted_lengths = np.array(extracted_lengths)

# Plot both distributions on log scale
ax.hist(np.log10(real_intergenic_distances + 1), bins=50, alpha=0.5, 
        label='Real (genomically consecutive)', color='green', edgecolor='black')
ax.hist(np.log10(extracted_lengths + 1), bins=50, alpha=0.5, 
        label='Extracted (pairs_df)', color='red', edgecolor='black')
ax.set_xlabel('Log10(Intergenic Distance + 1)', fontsize=12, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax.set_title('Real vs Extracted\n(Why extracted is wrong)', fontsize=13, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

# Add annotation showing the problem
ax.axvline(x=np.log10(100000), color='orange', linestyle='--', linewidth=2,
          label='100 kb threshold')
ax.text(np.log10(100000), ax.get_ylim()[1]*0.9, '← Biologically\nunrealistic', 
       ha='right', fontsize=10, color='red', fontweight='bold')

plt.tight_layout()
plt.savefig('graphs/intergenic_region_artifacts.png', dpi=300, bbox_inches='tight')
print("\n💾 Saved to 'graphs/intergenic_region_artifacts.png'")
plt.show()

# Summary
print("\n" + "=" * 80)
print("🎯 CONCLUSION: WHY 4.1 Mbp IS CLEARLY WRONG")
print("=" * 80)
print(f"""
1. E. coli genome is only {genome_size/1e6:.2f} Mbp total
   → A 4.1 Mbp "intergenic region" is 88% of the ENTIRE GENOME!

2. Real intergenic regions in bacteria:
   → Median: {np.median(real_intergenic_distances):.0f} bp
   → 95th percentile: {np.percentile(real_intergenic_distances, 95):.0f} bp
   → >1 kb is already unusual

3. These "intergenic regions" contain HUNDREDS of other genes!
   → They're not actually intergenic at all
   → The genes in the pair are NOT consecutive in the genome

4. Root cause: The create_pairs() function
   → Creates pairs from consecutive rows in the DataFrame
   → But DataFrame is ordered by protein ID, not genomic position
   → So "consecutive" proteins may be Mbp apart in the genome!

5. This is why TF binding sites appear:
   → More "intergenic" sequence = more chances for random TF motif matches
   → The 435 hits in 4.1 Mbp = 0.00011 hits/bp (very sparse!)
   → These are likely false positives in non-intergenic regions

BIOLOGICAL REALITY:
• Bacterial operons have 0-200 bp intergenic distances
• Anything > 500 bp is likely NOT an operon
• 4.1 Mbp is not an "intergenic region", it's the entire genome!
""")

