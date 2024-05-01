from Bio import SeqIO
import pandas as pd

refseq = ["RSV-A_NC_038235.1", "KT992094.1", "KU707921.1", "OM326756.1", "MW039343.1",
          "KF713490.1", "KF713492.1", "MT994243.1", "KX348546.1", "OR466361.1",
          "MZ516004.1", "MT506703.1", "OR466341.1", "MK733767.1", "MK733766.1", "KJ817798.1", "MK733769.1", "MK733768.1",
          "KJ817801.1", "KJ817800.1", "RSV-B_NC_001781.1", "OR466363.1", "MN365425.1", "MN365303.1", "MZ515670.1"]

# Read the FASTA file and store the sequences
alignment = SeqIO.parse("/home/tamengkel/compvir/mafft/rsv_sra_combined_mafft_msa.fasta", "fasta")
sequences = [str(record.seq) for record in alignment if record.id not in refseq]

# Get the length of the alignment
alignment_length = len(sequences[0])

# Get the total number of sequences
total_sequences = len(sequences)

# Initialize dictionaries to store nucleotide counts
nucleotide_counts = {"A": [], "T": [], "G": [], "C": []}

# Iterate over each position in the alignment
for i in range(alignment_length):
    # Initialize count variables for each nucleotide
    counts = {"A": 0, "T": 0, "G": 0, "C": 0}

    # Count the occurrences of each nucleotide at position i
    for sequence in sequences:
        nucleotide = sequence[i].upper()
        if nucleotide in counts:
            counts[nucleotide] += 1

    # Store the counts in the respective dictionary
    nucleotide_counts["A"].append(counts["A"])
    nucleotide_counts["T"].append(counts["T"])
    nucleotide_counts["G"].append(counts["G"])
    nucleotide_counts["C"].append(counts["C"])

# Create a dataframe from the nucleotide counts
df = pd.DataFrame(nucleotide_counts)

# Calculate the total amount of nucleotides at each position
df["Total"] = df[["A", "T", "G", "C"]].sum(axis=1)

# Add the position column
df["Position"] = range(1, alignment_length + 1)

# Reorder the column to include the new Total column
df = df[["Position", "A", "T", "G", "C", "Total"]]

# Save the dataframe as an Excel file
df.to_csv("/home/tamengkel/compvir/mafft/all_rsv_nucleotide_counts_240401.csv", index=False)

# Convert counts to percentages
'''df_percentage = df.copy()
for nucleotide in ["A", "T", "G", "C"]:
    df_percentage[nucleotide] = df_percentage[nucleotide] / df_percentage[["A", "T", "G", "C"]].sum(axis=1) * 100'''

