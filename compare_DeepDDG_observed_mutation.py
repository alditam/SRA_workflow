import pandas as pd

# Load the CSV files into pandas DataFrames
df_ddg = pd.read_csv('/home/tamengkel/compvir/deepddg_result_sync_position/sync_08_E2_GT2b_J8.csv')
df_mafft = pd.read_csv("/home/tamengkel/compvir/mafft/all_hcv_amino_acid_counts_240311.csv")

# Filter both DataFrames to only include positions from 380 to 770
df_ddg = df_ddg[(df_ddg['ResID'] >= 380) & (df_ddg['ResID'] <= 770)]
df_mafft = df_mafft[(df_mafft['Position'] >= 380) & (df_mafft['Position'] <= 770)]

# Initialize the DataFrame for Table C
df_c = pd.DataFrame(columns=[
    'Position',
    'WT',
    'Mutation',
    'Stabilizing Mutations',
    'Destabilizing mutations',
    'Positive ddG',
    'Negative ddG'
])

# Function to find the mutation with the highest positive ddG among observed mutations
def highest_positive_ddg_mutation(position, wild_type, observed_mutations):
    highest_ddg = 0
    highest_mutation = None
    for mutation in observed_mutations:
        if mutation == wild_type:
            continue
        mutation_ddg = df_ddg[(df_ddg['ResID'] == position) & (df_ddg['Mut'] == mutation)]['ddG'].max()
        if mutation_ddg > highest_ddg:
            highest_ddg = mutation_ddg
            highest_mutation = mutation
    return highest_mutation, highest_ddg

# Iterate through df_mafft for counting mutations
for _, mafft_row in df_mafft.iterrows():
    position = mafft_row['Position']
    wild_type = df_ddg[df_ddg['ResID'] == position]['WT'].iloc[0] if not df_ddg[df_ddg['ResID'] == position].empty else None

    # Skip if no wild type is identified for the position
    if wild_type is None:
        continue

    # Identify observed mutations from df_mafft
    observed_mutations = {aa for aa in 'ACDEFGHIKLMNPQRSTVWY' if mafft_row[aa] > 0}

    # Find the mutation with the highest positive ddG among observed mutations
    mutation, highest_ddg = highest_positive_ddg_mutation(position, wild_type, observed_mutations)

    # Count stabilizing and destabilizing mutations from the observed mutations only
    stabilizing_ddgs = df_ddg[(df_ddg['ResID'] == position) & (df_ddg['ddG'] > 0) & (df_ddg['Mut'].isin(observed_mutations))]
    destabilizing_ddgs = df_ddg[(df_ddg['ResID'] == position) & (df_ddg['ddG'] < 0) & (df_ddg['Mut'].isin(observed_mutations))]
    total_stabilizing = len(stabilizing_ddgs)
    total_destabilizing = len(destabilizing_ddgs)

    # Append the results to df_c
    df_c = df_c.append({
        'Position': position,
        'WT': wild_type,
        'Mutation': mutation if highest_ddg > 0 else 0,
        'Stabilizing Mutations': total_stabilizing,
        'Destabilizing mutations': total_destabilizing,
        'Positive ddG': highest_ddg if highest_ddg > 0 else 0,
        'Negative ddG': destabilizing_ddgs['ddG'].min() if not destabilizing_ddgs.empty else 0
    }, ignore_index=True)

# Save the resulting DataFrame to a TSV file
df_c.to_csv('/home/tamengkel/compvir/ergebnis/deepddg/deepddg_mafft_observed_mutations_GT2b_J8_240501.tsv', sep='\t', index=False)
