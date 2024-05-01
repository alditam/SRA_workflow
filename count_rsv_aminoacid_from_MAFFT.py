from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd

refseq = ["RSV-A_NC_038235.1", "KT992094.1", "KU707921.1", "OM326756.1", "MW039343.1",
          "KF713490.1", "KF713492.1", "MT994243.1", "KX348546.1", "OR466361.1",
          "MZ516004.1", "MT506703.1", "OR466341.1", "MK733767.1", "MK733766.1", "KJ817798.1", "MK733769.1", "MK733768.1",
          "KJ817801.1", "KJ817800.1", "RSV-B_NC_001781.1", "OR466363.1", "MN365425.1", "MN365303.1", "MZ515670.1"]

# Read the FASTA file and store the sequences
alignment = SeqIO.parse("/home/tamengkel/compvir/mafft/rsv_amino_sra_combined_mafft_msa.fst", "fasta")
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
df_aa.to_csv("/home/tamengkel/compvir/mafft/all_rsv_amino_acid_counts_240401.csv", index=False)

