import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Global variables
CPU = "32"
OUTDIR = "/mnt/twincore/compvironas/user/aldi/hcv_sra_pilot/"
HCV_DB = f"{OUTDIR}refseq/HCV_genomes_RefSeq.fasta"

# Read SRR IDs from a file
with open("SraAccList_1330-1809_240126.txt", "r") as f:
    srr_ids = f.read().strip().split("\n")

def prefetch(srr_id):
    prefetch_cmd = f"prefetch {srr_id} -O {OUTDIR}prefetch/sra/"
    subprocess.run(prefetch_cmd, shell=True)

def fastq_dump(srr_id):
    sra_file = f"{OUTDIR}prefetch/sra/{srr_id}.sra"
    subprocess.run(f"fastq-dump --split-spot --skip-technical --split-files {sra_file} --outdir {OUTDIR}fastq/", shell=True)

def fastp(srr_id):
    read1 = f"{OUTDIR}fastq/{srr_id}_1.fastq"
    read2 = f"{OUTDIR}fastq/{srr_id}_2.fastq"
    subprocess.run(f"fastp -i {read1} -I {read2} -o {OUTDIR}fastp/{srr_id}_1.fastq -O {OUTDIR}fastp/{srr_id}_2.fastq -q 20 -w 16", shell=True)

def spades(srr_id):
    read1 = f"{OUTDIR}fastp/{srr_id}_1.fastq"
    read2 = f"{OUTDIR}fastp/{srr_id}_2.fastq"
    subprocess.run(f"spades.py --rna -1 {read1} -2 {read2} -o {OUTDIR}spades/{srr_id}_spades/ --threads 16", shell=True)

def get_existing_transcript(srr_id):
    transcript_path = f"{OUTDIR}spades/{srr_id}_spades/transcripts.fasta"
    if os.path.isfile(transcript_path):
        return transcript_path
    return None

def blastn(srr_id):
    transcript_path = get_existing_transcript(srr_id)
    if transcript_path:
        blastn_outfile = f"{OUTDIR}blastn/{srr_id}_blastn"
        blastn_cmd = f"blastn -query {transcript_path} -db {HCV_DB} -out {blastn_outfile} -outfmt '6 qseqid qlen qstart qend sseqid slen sframe evalue score pident bitscore qseq' -perc_identity 80 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1"
        subprocess.run(blastn_cmd, shell=True)
        return blastn_outfile
    return None

def extract_sequence(srr_id):
    blastn_outfile = blastn(srr_id)
    if blastn_outfile:
        extract_outfile = f"{OUTDIR}blastn_seq_extract/{srr_id}_extract_blastn.fasta"
        subprocess.run(f"python3 extract_sequence_blastn.py {srr_id} {blastn_outfile} {extract_outfile}", shell=True)

def main():
    for srr_id in srr_ids:
        print(f"Processing {srr_id}...")
        prefetch(srr_id)
        fastq_dump(srr_id)
        fastp(srr_id)
        spades(srr_id)
        extract_sequence(srr_id)

if __name__ == "__main__":
    main()
