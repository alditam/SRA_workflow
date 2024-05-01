import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Define the start and end positions of each gene
genes = {
    'C': (0, 572),
    'E1': (573, 1148),
    'E2': (1149, 2237),
    'p7': (2238, 2426),
    'NS2': (2427, 3077),
    'NS3': (3078, 4969),
    'NS4A': (4969, 5132),
    'NS4B': (5133, 5915),
    'NS5A': (5916, 7260),
    'NS5B': (7260, 9032)
}

# Define colors for each gene for visibility
colors = {
    'UTR': 'grey',
    'C': 'yellow',
    'E1': 'orange',
    'E2': 'gold',
    'p7': 'beige',
    'NS2': 'yellowgreen',
    'NS3': 'lightgreen',
    'NS4A': 'lightsteelblue',
    'NS4B': 'cornflowerblue',
    'NS5A': 'salmon',
    'NS5B': 'plum'
}

# Load the dataframes from the CSV files
df_counts = pd.read_csv("/home/tamengkel/compvir/mafft/all_rsv_nucleotide_counts_240401.csv")

# Define a function to calculate Shannon entropy for a given distribution of nucleotides
def shannon_entropy(distribution):
    total_counts = np.sum(distribution)
    if total_counts == 0:
        return 0  # Return 0 entropy for a window with no counts
    probabilities = distribution / total_counts
    probabilities2 = np.where(probabilities == 0, 1e-10, probabilities)
    entropy = -np.sum(probabilities * np.log2(probabilities2))
    return entropy

# Apply shannon_entropy function to each row
df_counts["Entropy"] = df_counts.apply(lambda row: shannon_entropy(row[['A', 'T', 'G', 'C']]), axis=1)
df_counts = df_counts[['Position', 'A', 'T', 'G', 'C', "Total", "Entropy"]]

# Ensure "Entropy" values are in float
df_counts["Entropy"] = df_counts["Entropy"].astype(float)

# Define the window size
window_size = 27

# Apply rolling window to the "Entropy" column
df_counts["Smoothed_Entropy"] = df_counts["Entropy"].rolling(window=window_size, center=True, min_periods=1).mean()

# Filter for the positions start_position to end_position
start_position = 0
end_position = 9000
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
ax1.set_title('Hepatitis C Virus Genome Organization')

# Set the title and labels
ax2.set_title("Shannon Entropy")

# Calculate the tick positions and labels within the defined range
tick_labels = list(range(start_position, end_position + 1, step))

# Set the custom ticks and labels on the x-axis
ax2.set_xticks(tick_labels)
ax2.set_xticklabels(tick_labels, rotation=45)

# Add the entropy values to the plot on a secondary y-axis using the filtered_df DataFrame
#ax2 = ax.twinx()
ax2.set_xlabel('Nucleotides')
ax2.set_ylim(0,2)
ax2.plot(df_counts['Position'], df_counts['Smoothed_Entropy'], color='green')
ax2.set_ylabel('Shannon Entropy')
ax2.legend(loc='upper right')

# Plotting
# Create a figure with two subplots, one on top of the other
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
ax1.set_title('Hepatitis C Virus Genome Organization')

# Set the title and labels
ax2.set_title("Shannon Entropy")

# Calculate the tick positions and labels within the defined range
tick_labels = list(range(start_position, end_position + 1, step))

# Set the custom ticks and labels on the x-axis
ax2.set_xticks(tick_labels)
ax2.set_xticklabels(tick_labels, rotation=45)

# Add the entropy values to the plot on a secondary y-axis using the filtered_df DataFrame
#ax2 = ax.twinx()
ax2.set_xlabel('Nucleotides')
ax2.set_ylim(0,1)
ax2.plot(df_counts['Position'], df_counts['Smoothed_Entropy'], color='green')
#ax2.set_ylabel('Shannon Entropy')
ax2.legend(loc='upper right')
#df_counts.to_csv("counts_nucleotide_entropy.csv", index=False)

plt.show()