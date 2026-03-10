#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 10:39:35 2025

@author: nickduggett
"""


from Bio import SeqIO
import os
import subprocess

def split_fasta(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    contig_files = []

    for record in SeqIO.parse(input_file, "fasta"):
        output_path = os.path.join(output_dir, f"{record.id}.fasta")
        with open(output_path, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        contig_files.append(output_path)
        print(f"Written: {output_path}")

    return contig_files


def create_mash_sketch_from_directory(fasta_dir, sketch_output, kmer_size):
    # Get all .fasta files in the directory
    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    fasta_paths = [os.path.join(fasta_dir, f) for f in fasta_files]

    if not fasta_paths:
        raise FileNotFoundError("No .fasta files found in the directory.")

    # Run mash sketch
    cmd = ["mash", "sketch", "-k", str(kmer_size),"-S","22","-s","1000","-p","4","-o", sketch_output] + fasta_paths
    subprocess.run(cmd, check=True)
    print(f"Mash sketch created: {sketch_output}.msh")

import sys

#HACK: to allow calling from seqfinder.py and sketching of databases on the fly
if len(sys.argv)>1:
    kmer_size=str(16)
    input_fasta = sys.argv[1]
    output_directory = sys.argv[2]
    sketch_name = os.path.join(os.path.dirname(input_fasta),os.path.splitext(os.path.basename(input_fasta))[0]+"_k"+kmer_size)
else:  
    # Example usage
    kmer_size=str(16)
    input_fasta = "/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_6.0.0/references/AMR/AMR_with_plasmids_20250709.fna"
    output_directory=os.path.join(os.path.dirname(input_fasta),"split_genes")
    sketch_name = os.path.join(os.path.dirname(input_fasta),os.path.splitext(os.path.basename(input_fasta))[0]+"_k"+kmer_size)


contigs = split_fasta(input_fasta, output_directory)
create_mash_sketch_from_directory(output_directory, sketch_name, kmer_size)
