import pandas as pd
import sys

# ---------------------------- 
# Configuration
# ---------------------------- 
if len(sys.argv) < 3:
    print("Usage: python script.py <folder> <FULL> [chip_diff_threshold] [bad_ag_threshold] [good_ag_threshold]")
    sys.exit(1)

folder = sys.argv[1]
FULL = sys.argv[2]  # Required: "full" or other value
chip_diff_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
bad_ag_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
good_ag_threshold = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
exp = sys.argv[6] if len(sys.argv) > 6 else None

print(f"Arguments: folder={folder}, FULL={FULL}, chip_diff={chip_diff_threshold}, "
      f"bad_ag={bad_ag_threshold}, good_ag={good_ag_threshold}")

HEPG2_FILE = f'/mnt/altnas/work/ishawnia/empvsag/{FULL}-CTCF-HG2-15.tsv'
GM_FILE = f'/mnt/altnas/work/ishawnia/empvsag/{FULL}-CTCF-GM12878-15.tsv'

print(f"Processing with chip_diff_threshold={chip_diff_threshold}, bad_ag_threshold={bad_ag_threshold}, good_ag_threshold={good_ag_threshold}")

# ---------------------------- 
# Helper functions
# ---------------------------- 
def filter_high_chip_diff(df, threshold):
    """Keep only rows where abs(DIFF_LOG2_data1) > threshold."""
    return df[abs(df['DIFF_LOG2_data1']) > threshold].copy()

def create_good_bad_datasets(df, genome_name, good_ag_threshold=1, bad_ag_threshold=0.1):
    """
    Add binary columns for each group (one-hot encoding) and assign single prediction_quality label.
    A point can satisfy multiple group criteria, but gets ONE final label based on priority.
    
    Groups (can overlap):
    group_good: abs(DIFF_LOG2_data1 - DIFF_LOG2_data2) <= good_ag_threshold
    group_bad: abs(DIFF_LOG2_data2) < bad_ag_threshold
    group_cloud: both > 0.5 or both < -0.5
    
    Priority for prediction_quality label: good > cloud > bad > other
    """
    df = df.copy()
    df['GENOME'] = genome_name
    
    # Create binary columns for each group (one-hot encoding - can overlap)
    df['group_good'] = (abs(df['DIFF_LOG2_data1'] - df['DIFF_LOG2_data2']) <= good_ag_threshold).astype(int)
    df['group_bad'] = (abs(df['DIFF_LOG2_data2']) < bad_ag_threshold).astype(int)
    df['group_cloud'] = (
        ((df['DIFF_LOG2_data1'] > 0.5) & (df['DIFF_LOG2_data2'] > 0.5)) | 
        ((df['DIFF_LOG2_data1'] < -0.5) & (df['DIFF_LOG2_data2'] < -0.5))
    ).astype(int)
    
    # Assign single prediction_quality label with priority: good > cloud > bad > other
    df['prediction_quality'] = 'other'
    df.loc[df['group_bad'] == 1, 'prediction_quality'] = 'bad'
    df.loc[df['group_cloud'] == 1, 'prediction_quality'] = 'cloud'
    df.loc[df['group_good'] == 1, 'prediction_quality'] = 'good'
    
    # Split into three separate dataframes based on final label (mutually exclusive)
    good = df[df['prediction_quality'] == 'good'].copy()
    bad = df[df['prediction_quality'] == 'bad'].copy()
    cloud = df[df['prediction_quality'] == 'cloud'].copy()
    
    return df


def plot_simple_scatter(df, folder):
    """
    Simple scatter plot of empirical vs predicted, colored by prediction quality.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with DIFF_LOG2_data1, DIFF_LOG2_data2, and prediction_quality columns
    folder : str
        Folder path to save the figure
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Get unique prediction qualities and assign colors
    qualities = df['prediction_quality'].unique()
    colors = plt.cm.Set1(range(len(qualities)))
    color_map = dict(zip(qualities, colors))
    
    # Plot each quality group
    for quality in qualities:
        subset = df[df['prediction_quality'] == quality]
        ax.scatter(
            subset['DIFF_LOG2_data1'],
            subset['DIFF_LOG2_data2'],
            c=[color_map[quality]],
            label=quality,
            alpha=0.6,
            s=20
        )
    
    # Add y=x reference line
    ax.axline((0, 0), slope=1, color='gray', linestyle='--', alpha=0.5, label='y=x')
    
    ax.set_xlabel("Empirical Log2 Fold Change")
    ax.set_ylabel("Predicted Log2 Fold Change")
    ax.set_title(f"Empirical vs Predicted (n={len(df)})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{folder}/simple_quality_check.png', dpi=300)
    plt.close()
    
    print(f"Saved plot to {folder}/simple_quality_check.png")


# ---------------------------- 
# Load and process data
# ---------------------------- 
# Load dataframes
df_hepg2 = pd.read_csv(HEPG2_FILE, sep='\t')
df_gm = pd.read_csv(GM_FILE, sep='\t')
if exp:
    df_gm = df_gm[df_gm['transcription_factor_data1'].str.endswith(exp)].copy()

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

# Label datasets
df_hepg2_labeled = create_good_bad_datasets(df_hepg2_filtered, 'HepG2')
df_gm_labeled = create_good_bad_datasets(df_gm_filtered, 'GM12878')

print(f"\nDataset sizes:")
print(f"  HepG2 - good: {(df_hepg2_labeled['prediction_quality']=='good').sum()}, "
      f"bad: {(df_hepg2_labeled['prediction_quality']=='bad').sum()}, "
      f"cloud: {(df_hepg2_labeled['prediction_quality']=='cloud').sum()}, "
      f"other: {(df_hepg2_labeled['prediction_quality']=='other').sum()}")
print(f"  GM12878 - good: {(df_gm_labeled['prediction_quality']=='good').sum()}, "
      f"bad: {(df_gm_labeled['prediction_quality']=='bad').sum()}, "
      f"cloud: {(df_gm_labeled['prediction_quality']=='cloud').sum()}, "
      f"other: {(df_gm_labeled['prediction_quality']=='other').sum()}")

# ---------------------------- 
# Combine and save
# ---------------------------- 
# Combine all datasets
# combined_anchors = pd.concat([hepg2_good, hepg2_bad, gm_good, gm_bad], ignore_index=True)
combined_anchors = pd.concat([df_hepg2_labeled, df_gm_labeled])
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
plot_simple_scatter(combined_anchors, folder)