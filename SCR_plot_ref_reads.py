import sys
import pandas as pd

import matplotlib.pyplot as plt

# Read the dataframe
df = pd.read_csv(sys.argv[1], sep='\t')
folder = sys.argv[1].split('/')[0]

labels = {"1": "Reference has more reads", "0": "Alternative has more reads", "2": "Multiple SNPs / Phasings"}

# Get unique genomes
genomes = df['GENOME_x'].unique()

# Create figure with 2 rows (genomes) x 5 columns (plot types)
fig, axes = plt.subplots(len(genomes), 5, figsize=(25, 8))

# If only one genome, make sure axes is 2D
if len(genomes) == 1:
    axes = axes.reshape(1, -1)

for idx, genome in enumerate(genomes):
    genome_df = df[df['GENOME_x'] == genome]
    
    # Define the 5 different category combinations to plot
    plot_configs = [
        {'categories': [0, 1, 2], 'title_suffix': 'All Categories'},
        {'categories': [0, 1], 'title_suffix': 'Categories 0 & 1 Only'},
        {'categories': [0], 'title_suffix': 'Category 0 Only'},
        {'categories': [1], 'title_suffix': 'Category 1 Only'},
        {'categories': [2], 'title_suffix': 'Category 2 Only (Weird SNPs)'}
    ]
    
    cmap = plt.get_cmap('Set1')
    
    for plot_idx, config in enumerate(plot_configs):
        ax = axes[idx, plot_idx]  # Get the correct subplot
        ax.set_xlim(-6, 6)
        ax.set_ylim(-6, 6)
        
        categories = config['categories']
        
        # Filter data for this plot's categories
        plot_df = genome_df[genome_df['ref_more_reads'].isin(categories)]
        
        # Shuffle the dataframe to mix up plotting order
        plot_df_shuffled = plot_df.sample(frac=1, random_state=42)
        
        # Plot each category from the shuffled data
        for val in sorted(categories):
            subset = plot_df_shuffled[plot_df_shuffled['ref_more_reads'] == val]
            
            ax.scatter(
                subset['DIFF_LOG2_data1'],
                subset['DIFF_LOG2_data2'],
                color=cmap(val),
                alpha=0.5,
                label=labels.get(str(val)),
                s=20
            )
        
        ax.legend(title="ref_more_reads", loc='best')
        ax.axline((0, 0), slope=1, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel("Empirical Log2 Fold Change")
        ax.set_ylabel("Predicted Log2 Fold Change")
        perc_bad = len(plot_df[plot_df['prediction_quality'] == 'bad'])/len(plot_df) if len(plot_df) > 0 else 0
        perc_good = len(plot_df[plot_df['prediction_quality'] == 'good'])/len(plot_df) if len(plot_df) > 0 else 0
        ax.set_title(f"{genome}: {config['title_suffix']}\n(n={len(plot_df)}, good={perc_good:.2%}, bad={perc_bad:.2%})")
        ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{folder}/8-ref_more_reads.png", dpi=300)
plt.close()