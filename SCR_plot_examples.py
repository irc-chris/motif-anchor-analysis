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
    """
    # Create a confusion matrix / contingency table
    confusion_data = pd.crosstab(
        adf["strongest_effect"], 
        adf[x_param]
    )
    
    # Create heatmap
    sns.heatmap(
        confusion_data, 
        annot=True, 
        fmt='d', 
        cmap='YlOrRd', 
        ax=ax,
        cbar_kws={'label': 'Count'}
    )
    
    if unique:
        ax.set_title(f"Unique Anchors: SNP Effect vs {x_param}")
    else:
        ax.set_title("SNP Effect vs Overlapping Anchor Prediction")
    ax.set_xlabel(f"Is Overlapping Anchor Predicted {x_param}?")
    ax.set_ylabel("Strongest SNP Effect on Motif")


folder = sys.argv[1]  # folder to read / write files
# Concatenate all dataframes
combined_df_og = pd.read_csv(f"{folder}/1-motif_anchor_analysis_selected_anchors.tsv", sep="\t")
dfs = []

for genome in ['GM', 'HG2']:

    combined_df = combined_df_og[combined_df_og['GENOME'].str.startswith(genome[0])]

    print(f"Processing genome: {genome}", genome[0])

    combined_anchors_df = pd.read_csv(f"{folder}/3-{genome}_snp_motif_details.tsv", sep="\t", header=0)

    combined_with_snp_df = pd.merge(combined_df, combined_anchors_df,
                        on=['CHR', 'POS1', 'POS2'], how='left')

    # Convert float columns back to int if they were originally int and have no decimals
    for col in combined_anchors_df.columns:
        if col not in ['CHR', 'POS1', 'POS2'] and col in combined_with_snp_df.columns:
            if combined_anchors_df[col].dtype in ['int64', 'int32', 'int16', 'int8']:
                combined_with_snp_df[col] = combined_with_snp_df[col].fillna(0).astype(int)

    combined_with_snp_df = combined_with_snp_df.sort_values(by=['variant_id'])

    variant_ids = combined_with_snp_df['variant_id'].unique()
    variant_ids = np.sort(variant_ids)

    cmap = plt.cm.get_cmap('turbo', len(variant_ids))  # perceptually good
    variant_to_color = {
        vid: cmap(i)
        for i, vid in enumerate(variant_ids)
    }


    # Create the figure and two subplots
    fig, axes = plt.subplots(1, 5, figsize=(20, 6))

    # Define colors for each group
    colors = combined_with_snp_df['variant_id'].map(variant_to_color)

    # Plotting on the left subplot
    for group_name, group_df in combined_with_snp_df.groupby('GENOME_x'):
        axes[0].scatter(
            group_df['DIFF_LOG2_data1'],
            group_df['DIFF_LOG2_data2'],
            c=colors[group_df.index],
            label=group_name,
            alpha=0.7
        )

    print('Number of rows for genome:', combined_df.shape[0])
    print("Unique variant IDs:", len(variant_ids))
    print(combined_with_snp_df.shape)

    axes[0].set_xlabel("DIFF_LOG2_Emp")
    axes[0].set_ylabel("DIFF_LOG2_Pred")
    axes[0].set_title("Random Sample from 'along-x' and 'btwn_lines'\n from ChipSeq vs AG")
    axes[0].grid(True)

    unique_groups = [0,1,2]
    axes[1].boxplot(
        [combined_with_snp_df.loc[combined_with_snp_df["strongest_effect"] == g, "DIFF_LOG2_data2"] for g in unique_groups],
        labels=['Effect <1', 'Effect <5', 'Effect >=5']
    )
    axes[1].set_ylabel("Length")

    axes[1].set_title("Predicted Log2 Diff by SNP Effect")
    axes[1].set_xlabel("Strongest SNP Effect on Motif")
    axes[1].set_ylabel("DIFF_LOG2_Pred")


    axes[2].set_title("Predicted Log2 Diff vs SNP Effect")
    axes[2].scatter(
        combined_with_snp_df['strongest_effect'],
        combined_with_snp_df['DIFF_LOG2_data2'],
        c=colors,
        alpha=0.7
    )
    axes[2].set_xlim(-1,3)
    axes[2].set_xlabel("Strongest SNP Effect on Motif")
    axes[2].set_ylabel("DIFF_LOG2_Pred")

    plot_confusion_heatmap(axes[3], "btwn_lines", combined_with_snp_df)
    combined_with_snp_unique_anchors = combined_with_snp_df.drop_duplicates(subset=['CHR', 'POS1', 'POS2'])
    plot_confusion_heatmap(axes[4], "btwn_lines", combined_with_snp_unique_anchors, unique=True)

    plt.tight_layout()
    plt.savefig(f"{folder}/{genome}_motif_anchor_snp_effects_summary.png")
    dfs.append(combined_with_snp_df)

fig, ax = plt.subplots(figsize=(8, 6))
total = pd.concat(dfs, ignore_index=True)
plot_confusion_heatmap(ax, "btwn_lines", total)
ax.set_xlabel("Is Anchor Predicted btwn_lines?")
ax.set_ylabel("Strongest SNP Effect on Overlapping Motif")
ax.set_title("Good Prediction vs Strongest SNP Effect on Motif")
plt.tight_layout()
plt.savefig(f"{folder}/btwn_lines_vs_strongest_effect.png")