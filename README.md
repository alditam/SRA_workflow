# SRA_workflow
Python scripts explanation in pseudo-code format.

# files
references_hcv_blastn.fasA = HCV reference sequence for basic local alignment search with BLASTn
references_reporters_codon.fasA = HCV reference sequence with adjusted codon used for multiple sequence alignment with MAFFT
RSV-AB-CDSs-final.fasA = RSV reference sequence used for multiple sequence alignment with MAFFT
sra_snakemake.py   = snakemake based workflow to process SRA files
sra_workflow.py    = python based workflow to process SRA files

# sra_snakemake.py

Initialize necessary libraries and set global constants
Import essential libraries (subprocess, timeit, os, numpy, pandas, BioPython, threading, concurrent.futures)
Define CPU count and output directory constants

Read list of sample identifiers (SRR IDs) from a file
Read SRR IDs from "SraAccList.txt"

Define function to locate existing transcript files
Function get_existing_transcript(wildcards):
    Check if transcript files exist for a given SRR ID
    Return path to the transcript file if it exists

Snakemake rules for the workflow
Rule all:
    Define expected outputs for all samples based on their SRR IDs

Rule prefetch:
    Download SRA files for each SRR ID using the "prefetch" command

Rule fastq_dump:
    Convert SRA files to FASTQ format using the "fastq-dump" command
    Split reads into separate files for paired-end data

Rule fastp:
    Process FASTQ files with "fastp" to trim and filter reads
    Specify minimum quality score and number of threads

Checkpoint spades:
    Assemble reads into transcripts using "spades.py" in RNA mode
    Specify input paired-end files and output directory

Rule blastn:
    Perform BLAST search of assembled transcripts against a reference HCV database
    Filter results by identity and evalue, specify output format

Rule extract_sequence:
    Extract specific sequences from BLAST output using a custom script
    Produce a final FASTA file containing sequences of interest

Implementation details
Utilize threading and concurrency for efficient processing
Log all steps to monitor progress and troubleshoot issues


























































