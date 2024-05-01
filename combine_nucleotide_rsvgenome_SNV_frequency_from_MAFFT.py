import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Define the start and end positions of each gene for RSV
genes = {
    'NS1': (99, 518),
    'NS2': (626, 1000),
    'N': (1140, 2315),
    'P': (2348, 3073),
    'M': (3263, 4033),
    'SH': (4303, 4500),
    'G': (4690, 5589),
    'F': (5666, 7390),
    'M2': (7618, 8443),
    'L': (8509, 15009)
}

# Define colors for each gene for visibility
colors = {
    'NS1': 'orange',
    'NS2': 'red',
    'N': 'blue',
    'P': 'green',
    'M': 'purple',
    'SH': 'pink',
    'G': 'lightblue',
    'F': 'lightgreen',
    'M2': 'yellow',
    'L': 'magenta'
}

# Load the dataframes from the CSV files
df_counts = pd.read_csv("/home/tamengkel/compvir/mafft/all_rsv_nucleotide_counts_240401.csv")

print(df_counts)

# Calculate SNV ratio
df_counts["SNV"] = (df_counts["Total"] - df_counts[["A", "T", "G", "C",]].max(axis=1)) / df_counts["Total"]
df_counts['SNV_log'] = df_counts['SNV'].replace(0, 0.01)  # Replace 0 with a small number like 1e-6

# Create a new column for color
df_counts['Color'] = 'red'  # default color

# Assign colors based on position
df_counts.loc[df_counts['Position'] % 3 == 1, 'Color'] = 'red'
df_counts.loc[df_counts['Position'] % 3 == 2, 'Color'] = 'blue'
df_counts.loc[df_counts['Position'] % 3 == 0, 'Color'] = 'darkgoldenrod'

# Filter for the positions start_position to end_position
start_position = 0
end_position = 13500
step = 500

# Figure with two subplots, one on top of the other
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={'height_ratios': [1, 10]})

# Plot A (top plot)
for gene, (start, end) in genes.items():
    ax1.add_patch(Rectangle((start, 0), end-start, 1, color=colors[gene], label=gene))
    ax1.text((start+end)/2, 0.5, gene, ha='center', va='center', fontsize=10)

# Formatting for HCV genome plot
ax1.set_xlim(start_position, end_position)
ax1.set_ylim(0, 1)
ax1.set_yticklabels([])
ax1.set_xticks(range(start_position, end_position, step))
#ax1.set_xlabel('Nucleotides')
ax1.set_title('RSV Genome Organization')

# Set the title and labels
ax2.set_title("Shannon Entropy")

# Calculate the tick positions and labels within the defined range
tick_labels = list(range(start_position, end_position + 1, step))

# Set the custom ticks and labels on the x-axis
ax2.set_xticks(tick_labels)
ax2.set_xticklabels(tick_labels, rotation=45)

# Plotting
# Create a figure with two subplots, one on top of the other
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={'height_ratios': [1, 10]})

# Plot A (top plot)
for gene, (start, end) in genes.items():
    ax1.add_patch(Rectangle((start, 0), end-start, 1, color=colors[gene], label=gene))
    ax1.text((start+end)/2, 0.5, gene, ha='center', va='center', fontsize=10)

for color in ['red', 'blue', 'darkgoldenrod']:
    subset = df_counts[df_counts['Color'] == color]
    ax2.scatter(subset['Position'], subset['SNV_log'], alpha=0.5, s=3, color=color)

# Formatting for HCV genome plot
ax1.set_xlim(start_position, end_position)
ax1.set_ylim(0, 1)
ax1.set_yticklabels([])
ax1.set_xticks(range(start_position, end_position, step))
#ax1.set_xlabel('Nucleotides')
ax1.set_title('RSV Genome Organization')

# Calculate the tick positions and labels within the defined range
tick_labels = list(range(start_position, end_position + 1, step))

# Set the custom ticks and labels on the x-axis
ax2.set_xticks(tick_labels)
ax2.set_xticklabels(tick_labels, rotation=45)

# Add the entropy values to the plot on a secondary y-axis using the filtered_df DataFrame
#ax2 = ax.twinx()
ax2.set_yscale("log")
ax2.set_ylim(0.01, 1)
ax2.set_yticks([0.01, 0.05, 0.1, 0.5, 1])
ax2.set_yticklabels(['1%', "5%", '10%', "50%", '100%'])
#ax2.scatter(df_counts['Position'], df_counts['SNV_log'], alpha=0.5, s=3)  # s controls the size of dots, alpha controls transparency
ax2.set_xlabel('Position')
ax2.set_ylabel("SNV frequency")
#ax2.set_ylabel('Shannon Entropy')
#ax2.title('SNV Ratio')
#ax2.legend(loc='upper right')
#df_counts.to_csv("counts_nucleotide_entropy.csv", index=False)

plt.show()


