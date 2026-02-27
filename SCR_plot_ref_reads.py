import sys
import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

# Read the dataframe
df = pd.read_csv(sys.argv[1] + "/6-combined_anchor_snp_data.tsv", sep='\t')
folder = sys.argv[1]
thresh = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0

labels = {"1": "Reference has more reads", "0": "Alternative has more reads", "2": "Multiple SNPs / Phasings"}

# Get unique genomes
genomes = df['GENOME_y'].unique()

# Create figure
fig2, axes2 = plt.subplots(len(genomes), 5, figsize=(20, 9))

if len(genomes) == 1:
    axes2 = axes2.reshape(1, -1)

for idx, genome in enumerate(genomes):
    genome_df = df[df['GENOME_y'] == genome]
    
    # Define the 5 different category combinations to plot
    plot_configs = [
        {'categories': [0, 1, 2], 'title_suffix': 'All Categories'},
        {'categories': [0, 1], 'title_suffix': 'REF vs ALT Only'},
        {'categories': [0], 'title_suffix': 'ALT Only'},
        {'categories': [1], 'title_suffix': 'REF Only'},
        {'categories': [2], 'title_suffix': 'Weird SNP or Phased)'}
    ]
    
    cmap = plt.get_cmap('Set1')
    
    for plot_idx, config in enumerate(plot_configs):
        categories = config['categories']
        
        # Filter data for this plot's categories
        plot_df = genome_df[genome_df['ref_more_reads'].isin(categories)]
        
        # Shuffle the dataframe to mix up plotting order
        plot_df_shuffled = plot_df.sample(frac=1, random_state=42)
        
        # Calculate metrics
        along_x = abs(plot_df['DIFF_LOG2_data2']) < 0.1
        btwn_lines_mask = abs(plot_df['DIFF_LOG2_data1'] - plot_df['DIFF_LOG2_data2']) <= thresh
        cloud_consistency = (plot_df['DIFF_LOG2_data1'] > 0.5) & (plot_df['DIFF_LOG2_data2'] > 0.5) | (plot_df['DIFF_LOG2_data1'] < -0.5) & (plot_df['DIFF_LOG2_data2'] < -0.5)

        n_total = len(plot_df)

        # Stats for within ±thresh
        n_good_thresh = btwn_lines_mask.sum()
        perc_good_thresh = n_good_thresh / n_total if n_total > 0 else 0
        n_good_cloud_thresh = cloud_consistency.sum()
        perc_good_cloud_thresh = n_good_cloud_thresh / n_total if n_total > 0 else 0

        # Plot on figure (within ±thresh)
        ax2 = axes2[idx, plot_idx]
        for val in sorted(categories):
            subset = plot_df_shuffled[plot_df_shuffled['ref_more_reads'] == val]
            ax2.scatter(
                subset['DIFF_LOG2_data1'],
                subset['DIFF_LOG2_data2'],
                color=cmap(val),
                alpha=0.3,
                label=labels.get(str(val)),
                s=15
            )
        ax2.legend(title="ref_more_reads", loc='best')
        ax2.axline((0, 1), slope=1, color='gray', linestyle='--', alpha=0.5)
        ax2.axline((0, -1), slope=1, color='gray', linestyle='--', alpha=0.5)
        ax2.axline((0, 0.5), slope=0, color='gray', linestyle='--', alpha=0.5)
        ax2.axline((0, -0.5), slope=0, color='gray', linestyle='--', alpha=0.5)
        ax2.set_xlabel("Empirical Log2 Fold Change")
        ax2.set_ylabel("Predicted Log2 Fold Change")
        ax2.set_title(f"{genome}: {config['title_suffix']}\n(n={n_total}), within ±{thresh}: {n_good_thresh} [{perc_good_thresh:.1%}]\nClouds of Consistency: {n_good_cloud_thresh} [{perc_good_cloud_thresh:.1%}]")
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(-6,6)
        ax2.set_ylim(-7,7)

        ax2.add_patch(Rectangle((0.5, 0.5), 5.5, 5.5,
            linewidth=1, 
            edgecolor='black', 
            facecolor='none',
            linestyle='--',
            alpha=0.7))


fig2.tight_layout()
fig2.savefig(f'{folder}/8-ref-comparison_within{thresh}.png', dpi=300)