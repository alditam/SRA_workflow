import subprocess
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import timeit
import os
import numpy as np
import glob
from Bio import SeqIO
import pandas as pd


def parse_blastn_output(files):
    query_sequences = {}
    for file in files:
        try:
            with open(file, 'r') as f:
                for line in f:
                    columns = line.strip().split('\t')
                    query_id = columns[0]
                    sframe = int(columns[6])
                    qseq = columns[11]


                    if sframe == -1:
                        qseq = str(Seq(qseq).reverse_complement())

                    if query_id in query_sequences:
                        query_sequences[query_id].append(qseq)
                    else:
                        query_sequences[query_id] = [qseq]
                        
        except FileNotFoundError:
            print(f"Warning: File not found: {file}. Skipping.")
    return query_sequences

def write_fasta_file(file, sequences):
    seq_records = []
    for query_id, qseqs in sequences.items():
        for i, qseq in enumerate(qseqs):
            seq_records.append(SeqRecord(Seq(qseq), id=f"{query_id}_{i+1}", description=''))
    SeqIO.write(seq_records, file, "fasta")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_files_pattern output_fasta.fa")
        sys.exit(1)

    input_files_pattern = sys.argv[1]
    output_file = sys.argv[2]

    input_files = glob.glob(input_files_pattern)
    query_sequences = parse_blastn_output(input_files)
    write_fasta_file(output_file, query_sequences)

if __name__ == "__main__":
    main()



#/home/tamengkel/compvir_prak/transient_fasta

# read all files using globs


'''files = glob.glob("/home/tamengkel/compvir_prak/blastn/filtered/*.csv")

for file in files:
    # read as dataframe within a loop
    df1 = pd.read_csv(file, sep="\t",
                      names=["qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "evalue", "pident", "score", "qdiff", "qperc"])
    df2 = df1[["qseqid", "qlen", "qstart", "qend"]]
    print(df2)'''

