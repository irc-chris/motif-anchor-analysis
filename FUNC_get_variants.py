import pandas as pd

# Read original file (no header)
df = pd.read_csv('/mnt/altnas/work/ishawnia/empvsag/motif-anchor-analysis/show-ragini4/GM_anchors_with_motifs.tsv', sep='\t', header=None)

# Create new dataframe with required columns
output = pd.DataFrame({
    'variant_id': df[3] + '_' + df[7].astype(str) + '_' + df[8] + '_' + df[9],
    'CHROM': df[3],
    'POS': df[7],
    'REF': df[8],
    'ALT': df[9]
})

# Save to new TSV file
output.to_csv('/mnt/altnas/work/ishawnia/empvsag/motif-anchor-analysis/show-ragini4/GM_variants.tsv', sep='\t', index=False)