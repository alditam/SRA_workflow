import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

genes = {
    'NS1': (33, 172),
    'NS2': (208, 333),
    'N': (380, 771),
    'P': (782, 1024),
    'M': (1087, 1344),
    'SH': (1434, 1500),
    'G': (1563, 1863),
    'F': (1888, 2463),
    'M2': (2539, 2814),
    'L': (2836, 5003)
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
df_counts = pd.read_csv("/home/tamengkel/compvir/mafft/all_rsv_amino_acid_counts_240401.csv")  # Update path as necessary

# Shannon entropy calculation remains the same as in your original script
def shannon_entropy(distribution):
    total_counts = np.sum(distribution)
    if total_counts == 0:
        return 0  # Return 0 entropy for a window with no counts
    probabilities = distribution / total_counts
    probabilities2 = np.where(probabilities == 0, 1e-10, probabilities)
    entropy = -np.sum(probabilities * np.log2(probabilities2))
    return entropy

# Apply shannon_entropy function to each row
df_counts["Entropy"] = df_counts.apply(lambda row: shannon_entropy(row[['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']]), axis=1)

# Ensure "Entropy" values are in float
df_counts["Entropy"] = df_counts["Entropy"].astype(float)

# Define the window size
window_size = 9

# Apply rolling window to the "Entropy" column
df_counts["Smoothed_Entropy"] = df_counts["Entropy"].rolling(window=window_size, center=True, min_periods=1).mean()

# Make sure index of data frame in correct order
df_counts = df_counts[['Position', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', "Total", "Entropy", "Smoothed_Entropy"]]

# Filter for the positions start_position to end_position for RSV
start_position = 0
end_position = 4500  # Assuming end position based on RSV genome size
step = 150  # Adjust step as necessary for RSV

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 6), sharex=True, gridspec_kw={'height_ratios': [1, 10]})

# Plot RSV genome (top plot)
for gene, (start, end) in genes.items():
    ax1.add_patch(Rectangle((start, 0), end-start, 1, color=colors[gene], label=gene))
    ax1.text((start+end)/2, 0.5, gene, ha='center', va='center', fontsize=10)

# Formatting for RSV genome plot
ax1.set_xlim(start_position, end_position)
ax1.set_ylim(0, 1)
ax1.set_yticklabels([])
ax1.set_xticks(range(start_position, end_position, step))
ax1.set_title('Respiratory Syncytial Virus Genome Organization')

# entropy plotting RSV
tick_labels = list(range(start_position, end_position + 1, step))
ax2.set_xlabel("Position")
ax2.set_xticks(tick_labels)
ax2.set_xticklabels(tick_labels, rotation=45)
ax2.set_ylim(0, 4.32)  # Adjust max entropy value as needed
ax2.plot(df_counts['Position'], df_counts['Smoothed_Entropy'], color='green')
ax2.set_ylabel('Shannon Entropy')
ax2.legend(loc='upper right')

df_counts.to_csv("/home/tamengkel/compvir/mafft/240429_rsv_amino_entropy_values.csv", sep="\t", index=False)
#plt.show()