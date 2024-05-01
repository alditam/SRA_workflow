from Bio import SeqIO
import pandas as pd

refseq = ["HCV_GT1_NC_004102.1", "HCV_GT2_NC_009823.1", "HCV_GT2b_MG406988", "HCV_GT3_NC_009824.1", "HCV_GT4_NC_009825.1",
          "HCV_GT5_NC_009826.1", "HCV_GT6_NC_009827.1", "HCV_GT7_NC_030791.1", "HCV_GT1_H77_NC_038882.1", "E1E2_2r",
          "E1E2_GT1b-Con1", "E1E2_H2a-3", "E1E2_Jc1", "E1E2_H2b-4", "E1E2_H2b-5", "E1E2_2b-J8", "E1E2_2k", "E1E2_GT4a_ED43",
          "E1E2_GT5a_SA13", "E1E2_GT1b-J4", "E1E2_GT3a_S52", "E1E2_H77"]

# Read the FASTA file and store the sequences
alignment = SeqIO.parse("/home/tamengkel/compvir/mafft/240305_combined_hcv_mafft_msa.fasta", "fasta")
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
df.to_csv("/home/tamengkel/compvir/mafft/240305_all_hcv_nucleotide_counts.csv", index=False)

# Convert counts to percentages
'''df_percentage = df.copy()
for nucleotide in ["A", "T", "G", "C"]:
    df_percentage[nucleotide] = df_percentage[nucleotide] / df_percentage[["A", "T", "G", "C"]].sum(axis=1) * 100'''

