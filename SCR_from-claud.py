import pandas as pd
import sys

# ---------------------------- 
# Configuration
# ---------------------------- 
folder = sys.argv[1]
chip_diff_threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
bad_ag_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
good_ag_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
FULL = "outside" # "full"

HEPG2_FILE = f'/mnt/altnas/work/ishawnia/empvsag/{FULL}-CTCF-HG2-15.tsv'
GM_FILE = f'/mnt/altnas/work/ishawnia/empvsag/{FULL}-CTCF-GM12878-15.tsv'

print(f"Processing with chip_diff_threshold={chip_diff_threshold}, bad_ag_threshold={bad_ag_threshold}, good_ag_threshold={good_ag_threshold}")

# ---------------------------- 
# Helper functions
# ---------------------------- 
def filter_high_chip_diff(df, threshold):
    """Keep only rows where abs(DIFF_LOG2_data1) > threshold."""
    return df[abs(df['DIFF_LOG2_data1']) > threshold].copy()

def create_good_bad_datasets(df, genome_name):
    """
    Split dataframe into good and bad predictions.
    
    Good: abs(DIFF_LOG2_data1 - DIFF_LOG2_data2) <= 1
    Bad: abs(DIFF_LOG2_data2) < 0.1
    """
    df = df.copy()
    df['GENOME'] = genome_name
    
    # Good predictions: within ±1 of empirical
    good = df[abs(df['DIFF_LOG2_data1']) - abs(df['DIFF_LOG2_data2']) <= good_ag_threshold].copy()
    good['prediction_quality'] = 'good'
    
    # Bad predictions: predicted close to 0
    bad = df[abs(df['DIFF_LOG2_data2']) < bad_ag_threshold].copy()
    bad['prediction_quality'] = 'bad'
    
    return good, bad

# ---------------------------- 
# Load and process data
# ---------------------------- 
# Load dataframes
df_hepg2 = pd.read_csv(HEPG2_FILE, sep='\t')
df_gm = pd.read_csv(GM_FILE, sep='\t')

print(f"Loaded HepG2: {df_hepg2.shape[0]} rows")
print(f"Loaded GM12878: {df_gm.shape[0]} rows")

# Filter for high chip-seq differences
df_hepg2_filtered = filter_high_chip_diff(df_hepg2, chip_diff_threshold)
df_gm_filtered = filter_high_chip_diff(df_gm, chip_diff_threshold)

print(f"After filtering |DIFF_LOG2_data1| > {chip_diff_threshold}:")
print(f"  HepG2: {df_hepg2_filtered.shape[0]} rows")
print(f"  GM12878: {df_gm_filtered.shape[0]} rows")
df_hepg2_filtered.to_csv(f'{folder}/1-hepg2_filtered.tsv', sep='\t', index=False)
df_gm_filtered.to_csv(f'{folder}/1-gm12878_filtered.tsv', sep='\t', index=False)

# Create good and bad datasets
hepg2_good, hepg2_bad = create_good_bad_datasets(df_hepg2_filtered, 'HepG2')
gm_good, gm_bad = create_good_bad_datasets(df_gm_filtered, 'GM12878')

print(f"\nDataset sizes:")
print(f"  HepG2 good: {hepg2_good.shape[0]}, bad: {hepg2_bad.shape[0]}")
print(f"  GM12878 good: {gm_good.shape[0]}, bad: {gm_bad.shape[0]}")

# ---------------------------- 
# Combine and save
# ---------------------------- 
# Combine all datasets
combined_anchors = pd.concat([hepg2_good, hepg2_bad, gm_good, gm_bad], ignore_index=True)

# Add btwn_lines column (1 for good, 0 for bad)
combined_anchors['btwn_lines'] = (combined_anchors['prediction_quality'] == 'good').astype(int)

# Save combined file
combined_anchors.to_csv(
    f'{folder}/1-motif_anchor_analysis_selected_anchors.tsv', 
    sep='\t', 
    index=False
)

# Extract unique anchor positions for each genome
just_anchor_variants = combined_anchors[['CHR', 'POS1', 'POS2', 'GENOME']]

hepg2_variants = (just_anchor_variants[just_anchor_variants['GENOME'] == 'HepG2']
                  .drop(columns=['GENOME'])
                  .drop_duplicates())
hepg2_variants.to_csv(
    f'{folder}/1-hepg2_selected_anchors_variants.bed', 
    sep='\t', 
    index=False, 
    header=['#chrom', 'start', 'end']
)

gm_variants = (just_anchor_variants[just_anchor_variants['GENOME'] == 'GM12878']
               .drop(columns=['GENOME'])
               .drop_duplicates())
gm_variants.to_csv(
    f'{folder}/1-gm12878_selected_anchors_variants.bed', 
    sep='\t', 
    index=False, 
    header=['#chrom', 'start', 'end']
)

print(f"\nSaved files to {folder}/")
print(f"  Combined anchors: {combined_anchors.shape[0]} rows")
print(f"  HepG2 unique anchors: {hepg2_variants.shape[0]}")
print(f"  GM12878 unique anchors: {gm_variants.shape[0]}")