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

# Create two separate figures
fig1, axes1 = plt.subplots(len(genomes), 5, figsize=(20, 9))
fig2, axes2 = plt.subplots(len(genomes), 5, figsize=(20, 9))

if len(genomes) == 1:
    axes1 = axes1.reshape(1, -1)
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
        
        # Calculate both metrics
        within1_mask = abs(plot_df['DIFF_LOG2_data1'] - plot_df['DIFF_LOG2_data2']) <= 1
        along_x = abs(plot_df['DIFF_LOG2_data2']) < 0.1
        btwn_lines_mask = abs(plot_df['DIFF_LOG2_data1'] - plot_df['DIFF_LOG2_data2']) <= thresh
        cloud_consistency = (plot_df['DIFF_LOG2_data1'] > 0.5) & (plot_df['DIFF_LOG2_data2'] > 0.5) | (plot_df['DIFF_LOG2_data1'] < -0.5) & (plot_df['DIFF_LOG2_data2'] < -0.5)
        
        n_total = len(plot_df)
        
        # Stats for within ±1
        n_good_1 = within1_mask.sum()
        perc_good_1 = n_good_1 / n_total if n_total > 0 else 0
        n_along_x = along_x.sum()
        perc_along_x = n_along_x / n_total if n_total > 0 else 0
        n_good_cloud = cloud_consistency.sum()
        perc_good_cloud = n_good_cloud / n_total if n_total > 0 else 0
        
        # Stats for within ±thresh
        n_good_thresh = btwn_lines_mask.sum()
        perc_good_thresh = n_good_thresh / n_total if n_total > 0 else 0
        n_along_x_thresh = along_x.sum()  # Same along_x applies
        perc_along_x_thresh = n_along_x_thresh / n_total if n_total > 0 else 0
        n_good_cloud_thresh = cloud_consistency.sum()
        perc_good_cloud_thresh = n_good_cloud_thresh / n_total if n_total > 0 else 0

        # Plot on first figure (within ±1)
        ax1 = axes1[idx, plot_idx]
        for val in sorted(categories):
            subset = plot_df_shuffled[plot_df_shuffled['ref_more_reads'] == val]
            ax1.scatter(
                subset['DIFF_LOG2_data1'],
                subset['DIFF_LOG2_data2'],
                color=cmap(val),
                alpha=0.3,
                label=labels.get(str(val)),
                s=15
            )
        ax1.legend(loc='best')
        ax1.axline((0, 1), slope=1, color='gray', linestyle='--', alpha=0.5)
        ax1.axline((0, -1), slope=1, color='gray', linestyle='--', alpha=0.5)
        ax1.axline((0, 0.5), slope=0, color='gray', linestyle='--', alpha=0.5)
        ax1.axline((0, -0.5), slope=0, color='gray', linestyle='--', alpha=0.5)
        ax1.set_xlabel("Empirical Log2 Fold Change")
        ax1.set_ylabel("Predicted Log2 Fold Change")
        ax1.set_title(f"{genome}: {config['title_suffix']} (n={n_total}),\n Gradient of Goodness: {n_good_1} [{perc_good_1:.1%}]\n Belt of Badness: {n_along_x} [{perc_along_x:.1%}]\n Clouds of Consistency: {n_good_cloud} [{perc_good_cloud:.1%}]")
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(-6,6)
        ax1.set_ylim(-7,7)
        
        # Plot on second figure (within ±thresh)
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

        # For axes1
        ax1.add_patch(Rectangle((0.5, 0.5), 5.5, 5.5, 
            linewidth=1, 
            edgecolor='black', 
            facecolor='none',
            linestyle='--',
            alpha=0.7))

        ax1.add_patch(Rectangle((-0.5, -0.5), -5.5, -5.5, 
            linewidth=1, 
            edgecolor='black', 
            facecolor='none',
            linestyle='--',
            alpha=0.7))

        # For axes2
        ax2.add_patch(Rectangle((0.5, 0.5), 5.5, 5.5, 
            linewidth=1, 
            edgecolor='black', 
            facecolor='none',
            linestyle='--',
            alpha=0.7))


# Save both figures
fig1.tight_layout()
fig1.savefig(f'{folder}/8-ref-comparison_within1.png', dpi=300)

fig2.tight_layout()
fig2.savefig(f'{folder}/8-ref-comparison_within{thresh}.png', dpi=300)