import subprocess
import timeit
import seqExtraction_usingIdenticalHeader
import os
import numpy as np
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import threading
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor

# globale variable für threads/anzahl cpus = 32
CPU = "32"

# output directory
OUTDIR =  "/mnt/twincore/compvironas/user/aldi/hcv_sra_pilot/"

# Read SRR IDs from a file
with open("SraAccList_1330-1809_240126.txt", "r") as f:
    ID_list = f.read().strip().split("\n")

srr_id = ID_list

def get_existing_transcript(wildcards):
    checkpoint_output = checkpoints.spades.get(srr_id=wildcards.srr_id).output.spades_outdir
    if os.path.isfile(f"{OUTDIR}spades/{wildcards.srr_id}_spades/transcripts.fasta"):
        found_files = glob_wildcards(os.path.join(checkpoint_output, f"{OUTDIR}spades/{srr_id}_spades/transcripts.fasta"))
        results = f"{OUTDIR}spades/{wildcards.srr_id}_spades/transcripts.fasta"
    return results



rule all:
    input:
        spades_outdir = expand(f"{OUTDIR}spades/{{srr_id}}_spades/", srr_id=ID_list)

rule prefetch:
    output:
        sra = f"{OUTDIR}prefetch/sra/{{srr_id}}.sra"
    shell:
        "prefetch {wildcards.srr_id}"

rule fastq_dump:
    input:
        sra = f"{OUTDIR}prefetch/sra/{{srr_id}}.sra"
    output:
        fastq1 = (f"{OUTDIR}fastq/{{srr_id}}_1.fastq"),
        fastq2 = (f"{OUTDIR}fastq/{{srr_id}}_2.fastq")
    shell:
        "fastq-dump --split-spot --skip-technical --split-files {input.sra} --outdir {OUTDIR}fastq/"

rule fastp:
    input:
        read1 = f"{OUTDIR}fastq/{{srr_id}}_1.fastq",
        read2 = f"{OUTDIR}fastq/{{srr_id}}_2.fastq"
    output:
        read1 = f"{OUTDIR}fastp/{{srr_id}}_1.fastq",
        read2 = f"{OUTDIR}fastp/{{srr_id}}_2.fastq"
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} -q 20 -w 16"

checkpoint spades:
    input:
        read1 = f"{OUTDIR}fastp/{{srr_id}}_1.fastq",
        read2 = f"{OUTDIR}fastp/{{srr_id}}_2.fastq"
    output:
        spades_outdir = directory(f"{OUTDIR}spades/{{srr_id}}_spades/")
    shell:
        "spades.py --rna -1 {input.read1} -2 {input.read2} -o {output.spades_outdir} --threads 16"

rule blastn:
    input:
        spades_transcript = get_existing_transcript,
        hcv_db = f"{OUTDIR}refseq/HCV_genomes_RefSeq.fasta"
    output:
        blastn_outfile = f"{OUTDIR}blastn/{{srr_id}}_blastn"
    shell:
        "blastn -query {input.spades_transcript} -db {input.hcv_db} -out {output.blastn_outfile} -outfmt '6 qseqid qlen qstart qend sseqid slen sframe evalue score pident bitscore qseq' -perc_identity 80 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1"


rule extract_sequence:
    input:
        blastn_outfile = f"{OUTDIR}blastn/{{srr_id}}_blastn"
    output:
        extract_outfile = f"{OUTDIR}blastn_seq_extract/{{srr_id}}_extract_blastn.fasta"
    shell:
        "python3 extract_sequence_blastn.py {wildcards.srr_id} {input.blastn_outfile} {output.extract_outfile}"

