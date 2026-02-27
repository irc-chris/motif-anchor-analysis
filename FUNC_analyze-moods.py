import sys
import pandas as pd

import matplotlib.pyplot as plt

# Read the CSV file from command line argument
df = pd.read_csv(sys.argv[1], sep='\t')
genome = "HG2" # "GM"
# Create the plots
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Histogram of delta_score
axes[0].hist(df['delta_score'], bins=30, edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Delta Score')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Distribution of Delta Score')
axes[0].grid(axis='y', alpha=0.3)

# Plot 2: Box plot of delta_score by strongest_effect
df.boxplot(column='delta_score', by='strongest_effect', ax=axes[1])
axes[1].set_xlabel('Strongest Effect')
axes[1].set_ylabel('Delta Score')
axes[1].set_title('Delta Score by Strongest Effect n={}')
plt.sca(axes[1])
plt.xticks(rotation=45, ha='right')

plt.tight_layout()

# Save the plot
plt.savefig(f"{genome}-delta_score_distribution.png", dpi=300, bbox_inches='tight')
print(f"Plot saved as '{genome}-new-delta_score_distribution.png'")