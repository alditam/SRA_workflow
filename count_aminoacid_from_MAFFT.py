from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd

refseq = ["HCV_GT1_NC_004102.1", "HCV_GT2_NC_009823.1", "HCV_GT2b_MG406988", "HCV_GT3_NC_009824.1", "HCV_GT4_NC_009825.1",
          "HCV_GT5_NC_009826.1", "HCV_GT6_NC_009827.1", "HCV_GT7_NC_030791.1", "HCV_GT1_H77_NC_038882.1", "E1E2_2r",
          "E1E2_GT1b-Con1", "E1E2_H2a-3", "E1E2_Jc1", "E1E2_H2b-4", "E1E2_H2b-5", "E1E2_2b-J8", "E1E2_2k", "E1E2_GT4a_ED43",
          "E1E2_GT5a_SA13", "E1E2_GT1b-J4", "E1E2_GT3a_S52", "E1E2_H77"]

# Read the FASTA file and store the sequences
alignment = SeqIO.parse("/home/tamengkel/compvir/mafft/240305_hcv_all_amino_combined_hcv_mafft_msa.fst", "fasta")
sequences = [str(record.seq) for record in alignment if record.id not in refseq]

# Get the length of the alignment
alignment_length = len(sequences[0])

# Initialize a dictionary to store amino acid counts
amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # List of amino acids
amino_acid_counts = {aa: [0] * alignment_length for aa in amino_acids}

# Iterate over each position in the alignment
for i in range(alignment_length):
    # Initialize count variables for each amino acid
    counts = {aa: 0 for aa in amino_acids}

    # Count the occurrences of each amino acid at position i
    for sequence in sequences:
        amino_acid = sequence[i].upper()
        if amino_acid in counts:
            counts[amino_acid] += 1

    # Store the counts in the respective dictionary
    for aa in amino_acids:
        amino_acid_counts[aa][i] = counts[aa]

# Create a dataframe from the amino acid counts
df_aa = pd.DataFrame(amino_acid_counts)

# Calculate the total amount of nucleotides at each position
df_aa["Total"] = df_aa[['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']].sum(axis=1)

# Add the position column
df_aa["Position"] = range(1, alignment_length + 1)

# Reorder the columns
df_aa = df_aa[["Position", 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', "Total"]]

# Save the dataframe as a CSV file
df_aa.to_csv("/home/tamengkel/compvir/mafft/all_hcv_amino_acid_counts_240311.csv", index=False)

