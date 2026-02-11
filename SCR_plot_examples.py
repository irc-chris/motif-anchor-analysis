import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np
import seaborn as sns


def plot_confusion_heatmap(ax, x_param, adf, unique=False):
    """
    Plot a confusion matrix heatmap.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to plot on
    x_param : str
        Column name for x-axis (columns of confusion matrix)
    adf : pandas.DataFrame
        DataFrame containing the data
    unique : bool
        Whether this is unique anchors only
    """
    # Check if dataframe is empty
    if adf.empty:
        ax.text(0.5, 0.5, 'No data available', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title("No Data")
        return
    
    # Create confusion matrix
    confusion_data = pd.crosstab(
        adf["strongest_effect"], 
        adf[x_param]
    )
    
    # Check if confusion matrix is empty
    if confusion_data.empty or confusion_data.size == 0:
        ax.text(0.5, 0.5, 'No data available', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title("No Data")
        return
    
    # Create heatmap
    sns.heatmap(
        confusion_data, 
        annot=True, 
        fmt='d', 
        cmap='YlOrRd', 
        ax=ax,
        cbar_kws={'label': 'Count'}
    )
    
    # Set titles with sample size
    n = len(adf)
    anchor_type = "Unique Anchors" if unique else "All Anchor-Motif Overlaps"
    ax.set_title(f"{anchor_type}\nSNP Effect vs Prediction Accuracy (n={n})")
    ax.set_xlabel("Predicted Within ±1 of Empirical?")
    ax.set_ylabel("Strongest SNP Effect on CTCF Motif")


def plot_combined_boxplot(ax, data, snp_effect_groups, effect_labels, title_prefix, positions_emp, positions_pred):
    """
    Plot combined empirical and predicted boxplots by SNP effect.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to plot on
    data : pandas.DataFrame
        DataFrame containing the data
    snp_effect_groups : list
        List of SNP effect values to group by
    effect_labels : list
        Labels for each SNP effect group
    title_prefix : str
        Prefix for the plot title
    positions_emp : list
        X-axis positions for empirical boxplots
    positions_pred : list
        X-axis positions for predicted boxplots
    """
    # Prepare data for boxplots
    boxplot_data_emp = [
        data.loc[data["strongest_effect"] == g, "DIFF_LOG2_data1"] 
        for g in snp_effect_groups
    ]
    boxplot_data_pred = [
        data.loc[data["strongest_effect"] == g, "DIFF_LOG2_data2"] 
        for g in snp_effect_groups
    ]
    
    # Plot empirical boxplots
    bp1 = ax.boxplot(boxplot_data_emp, positions=positions_emp, widths=0.6,
                    patch_artist=True, 
                    boxprops=dict(facecolor='lightblue', alpha=0.7),
                    medianprops=dict(color='darkblue', linewidth=2))
    
    # Plot predicted boxplots
    bp2 = ax.boxplot(boxplot_data_pred, positions=positions_pred, widths=0.6,
                    patch_artist=True,
                    boxprops=dict(facecolor='lightcoral', alpha=0.7),
                    medianprops=dict(color='darkred', linewidth=2))
    
    # Get n's for each group
    n_per_group = [len(d) for d in boxplot_data_emp]

    # Set x-axis labels
    ax.set_xticks([1.5, 4.5, 7.5])
    ax.set_xticklabels(effect_labels)
    ax.set_xlabel("Strongest SNP Effect on CTCF Motif")
    ax.set_ylabel("Log2 Fold Change")
    ax.set_title(f"{title_prefix}: Empirical vs Predicted Strength by SNP Effect\n(n={', '.join(map(str, n_per_group))})")
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', alpha=0.7, label='Empirical'),
        Patch(facecolor='lightcoral', alpha=0.7, label='Predicted')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

# ---------------------------- 
# Configuration
# ---------------------------- 
folder = sys.argv[1]
GENOMES = ['GM', 'HG2']
SNP_EFFECT_LABELS = ['Effect <1', 'Effect 1-5', 'Effect ≥5']

# ---------------------------- 
# Load data
# ---------------------------- 
combined_df_og = pd.read_csv(f"{folder}/1-motif_anchor_analysis_selected_anchors.tsv", sep="\t")
all_genome_dfs = []

# ---------------------------- 
# Process each genome
# ---------------------------- 
for genome in GENOMES:
    # Filter by genome
    combined_df = combined_df_og[combined_df_og['GENOME'].str.startswith(genome[0])]
    print(f"\nProcessing genome: {genome}")
    print(f"  Base anchors: {combined_df.shape[0]}")
    
    # Load SNP motif data
    combined_anchors_df = pd.read_csv(
        f"{folder}/3-{genome}_snp_motif_details.tsv", 
        sep="\t", 
        header=0
    )
    
    # Merge anchor and SNP data
    combined_with_snp_df = pd.merge(
        combined_df, 
        combined_anchors_df,
        on=['CHR', 'POS1', 'POS2'], 
        how='left'
    )
    
    # Preserve integer types
    for col in combined_anchors_df.columns:
        if col not in ['CHR', 'POS1', 'POS2'] and col in combined_with_snp_df.columns:
            if combined_anchors_df[col].dtype in ['int64', 'int32', 'int16', 'int8']:
                combined_with_snp_df[col] = combined_with_snp_df[col].fillna(0).astype(int)
    
    combined_with_snp_df = combined_with_snp_df.sort_values(by=['variant_id'])
    
    # Create color mapping for variants
    variant_ids = np.sort(combined_with_snp_df['variant_id'].unique())
    cmap = plt.cm.get_cmap('turbo', len(variant_ids))
    variant_to_color = {vid: cmap(i) for i, vid in enumerate(variant_ids)}
    colors = combined_with_snp_df['variant_id'].map(variant_to_color)
    
    print(f"  Unique variants: {len(variant_ids)}")
    print(f"  Total data points: {combined_with_snp_df.shape[0]}")
    
    # ---------------------------- 
    # Create plots (2 rows x 3 columns)
    # ---------------------------- 
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1 (0,0): Empirical vs Predicted Log2FC
    for group_name, group_df in combined_with_snp_df.groupby('GENOME_x'):
        axes[0,0].scatter(
            group_df['DIFF_LOG2_data1'],
            group_df['DIFF_LOG2_data2'],
            c=colors[group_df.index],
            label=group_name,
            alpha=0.7,
            s=30
        )
    n_plot1 = len(combined_with_snp_df)
    axes[0,0].axline((0, 0), slope=1, color='gray', linestyle='--', alpha=0.5, label='y=x')
    axes[0,0].set_xlabel("Empirical Log2 Fold Change")
    axes[0,0].set_ylabel("Predicted Log2 Fold Change")
    axes[0,0].set_title(f"Empirical vs Predicted Anchor Strength\n(n={n_plot1})")
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].legend()


    # Plot 2 (0,1) and Plot 3 (0,2): Combined boxplots
    good_pred = combined_with_snp_df[combined_with_snp_df['btwn_lines'] == 1]
    bad_pred = combined_with_snp_df[combined_with_snp_df['btwn_lines'] == 0]
    unique_groups = [0, 1, 2]
    positions_emp = [1, 4, 7]
    positions_pred = [2, 5, 8]


    plot_combined_boxplot(axes[0,1], bad_pred, unique_groups, SNP_EFFECT_LABELS, 
                        "Bad Predictions", positions_emp, positions_pred)
    plot_combined_boxplot(axes[0,2], good_pred, unique_groups, SNP_EFFECT_LABELS, 
                        "Good Predictions", positions_emp, positions_pred)


    # Plot 4 (1,0): Scatter of SNP effect vs predicted values
    n_plot3 = len(combined_with_snp_df)
    axes[1,0].scatter(
        combined_with_snp_df['strongest_effect'],
        combined_with_snp_df['DIFF_LOG2_data2'],
        c=colors,
        alpha=0.6,
        s=30
    )
    axes[1,0].set_xlim(-0.5, 2.5)
    axes[1,0].set_xticks([0, 1, 2])
    axes[1,0].set_xticklabels(['<1', '1-5', '≥5'])
    axes[1,0].set_xlabel("Strongest SNP Effect on CTCF Motif")
    axes[1,0].set_ylabel("Predicted Log2 Fold Change")
    axes[1,0].set_title(f"SNP Effect vs Predicted Strength\n(n={n_plot3})")
    axes[1,0].grid(True, alpha=0.3, axis='y')
    
    # Plot 5 (1,1): Confusion matrix - all overlaps
    plot_confusion_heatmap(axes[1,1], "btwn_lines", combined_with_snp_df, unique=False)
    
    # Plot 6 (1,2): Confusion matrix - unique anchors
    combined_with_snp_unique_anchors = combined_with_snp_df.drop_duplicates(subset=['CHR', 'POS1', 'POS2'])
    plot_confusion_heatmap(axes[1,2], "btwn_lines", combined_with_snp_unique_anchors, unique=True)
    
    plt.tight_layout()
    plt.savefig(f"{folder}/6-{genome}_motif_anchor_snp_effects_summary.png", dpi=300)
    plt.close()
    
    all_genome_dfs.append(combined_with_snp_df)

# ---------------------------- 
# Combined genome plot
# ---------------------------- 
fig, ax = plt.subplots(figsize=(8, 6))
total_df = pd.concat(all_genome_dfs, ignore_index=True)
plot_confusion_heatmap(ax, "btwn_lines", total_df, unique=False)

plt.tight_layout()
plt.savefig(f"{folder}/6-combined_btwn_lines_vs_strongest_effect.png", dpi=300)
plt.close()

print(f"\nAll plots saved to {folder}/")