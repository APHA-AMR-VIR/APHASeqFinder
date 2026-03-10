#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Wed Aug  6 10:37:38 2025

@author: grahamhill

"""
import sys
import os
from Bio import SeqIO

if len(sys.argv)>1:
    input_db = sys.argv[1]
    output_db = sys.argv[2]
else:
    input_db = "/home/grahamhill/APHASeqFinder/wgsim_of_db/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6.fna"
    output_db = "/home/grahamhill/APHASeqFinder/wgsim_of_db/AMRDatabase_20230526_and_EnteroPLasmids_20230526_short-tetA6_cdhit80_clustered.fna"
#cdhit
#currently set to only cluster things of same length, which is less than ideal
#Have to figure out a method for gapped alignment that works with SeqFinder v abricate
if not os.path.exists(output_db):
    os.system("cd-hit-est -i "+input_db+" -o "+output_db+" -T 4 -uS 0 -aL 1")

clusters = {}
core = []
with open(output_db+".clstr","r") as file:
    for line in file:
        if line.startswith(">Cluster"):
            cluster = line.strip()
            this_list = [] 
        else:
            this_list.append(line.split(",")[1][1:13])
            clusters.update({cluster:this_list})
        if "*" in line:
            core.append(line.split(",")[1][1:13])
            
db = SeqIO.parse(input_db,'fasta')

db = list(db)

dbdict = {}
name_cluster = {}
for gene in db:
    dbdict.update({">"+gene.id[0:11]:gene})
    name = gene.id[12:16]
    if name not in name_cluster:
        name_cluster.update({name:[gene]})
    else:
        this_list = name_cluster[name]
        this_list.append(gene)
        name_cluster.update({name:this_list})

clustered_db = SeqIO.parse(output_db,'fasta')
clustered_db = list(clustered_db)

#TODO: system for dealing with nucl and peptide mutations
#TODO: make a dictionary of SNPs for each clustered gene
#TODO: write a function to score the SNPs of mapped genes and pick the best as the true hit
#TODO: write a function to go through the output files and re-call based on above
"""
from Bio import Align
reference = "ATA"
example =    "AAA"
aligner = Align.PairwiseAligner(scoring='blastn')
aligner.mode = 'global'
print("Aligning sequences with "+str(aligner.algorithm))
alignments = aligner.align(example,reference)
for alignment in alignments:
    print(alignment)
"""
def identify_snps(ref_seq, gene_seq):
    if len(ref_seq) != len(gene_seq):
        #TODO: full alignment
        print("Not same length, won't work")
        return None
    #can't handle degenerate bases
    '''codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }'''
    syn_snps = []
    nonsyn_snps = []
    
    k = 0
    for i in range(0, len(ref_seq), 3):
        #codon tracker
        k = k+1
        ref_codon = ref_seq[i:i+3]
        gene_codon = gene_seq[i:i+3]
        
        # Skip incomplete codons at end
        if len(ref_codon) < 3:
            for j in range(0, len(ref_codon)):
                #base tracker
                l = ((k-1)*3)+j+1
                if not ref_codon[j] == gene_codon[j]:
                    syn_snps.append("nt"+str(ref_codon[j])+str(l)+str(gene_codon[j]))
            return syn_snps, nonsyn_snps  

        if ref_codon != gene_codon:
            ref_aa = SeqIO.SeqRecord(ref_codon).translate(table=11)
            gene_aa = SeqIO.SeqRecord(gene_codon).translate(table=11)
            #ref_aa = codontable[ref_codon]
            #gene_aa = codontable[gene_codon]
            #if ref_aa == gene_aa:
            if ref_aa.seq == gene_aa.seq:
                for j in range(0, 3):
                    #base tracker
                    l = ((k-1)*3)+j+1
                    if not ref_codon[j] == gene_codon[j]:
                        syn_snps.append("nt"+str(ref_codon[j])+str(l)+str(gene_codon[j]))
            else:
                #nonsyn_snps.append(ref_aa+str(k)+gene_aa)
                nonsyn_snps.append(str(ref_aa.seq)+str(k)+str(gene_aa.seq))

    return syn_snps, nonsyn_snps
"""
# Example input
reference = "AAA"
example =    "TTT"

syn_snps, nonsyn_snps = identify_snps(reference, example)
print("Synonymous SNPs:", syn_snps)
print("Non-synonymous SNPs:", nonsyn_snps)
"""
core_clustered = {}
for core_gene in core:
    for clust in clusters:
        if core_gene in clusters[clust]:
            core_clustered[core_gene] = clusters[clust]
 
"""
degenerates = set("RYSWKMBDHVN")
for gene in dbdict.values():
    for c in degenerates:
        if c in gene.seq:
            print(c)
            print(gene.seq.translate())
            if "X" in gene.seq.translate():
                print("!!!!!!!!!!!!!!!!!!")
            print(gene.id)
"""

table = []
for core_gene in core_clustered:
    gene_list = core_clustered[core_gene]
    if len(gene_list) > 1:
        ref_gene = dbdict.get(core_gene)
        #print(ref_gene.id)
        for gene in gene_list:
            if not gene == core_gene:
                this_gene = dbdict.get(gene)
                syn_snps, nonsyn_snps = identify_snps(ref_gene.seq, this_gene.seq)
                table.append([ref_gene.id,this_gene.id,",".join(syn_snps),",".join(nonsyn_snps)])
          
import pandas as pd
table = pd.DataFrame(table,columns=['ref_gene','grp_gene','syn_snps','nonsyn_snps'])
#TODO: fake and cheap alignment scoring by subtracting snps from found snps and then calc the best match by BLOSUM90, nucl tiebreak


def number_to_code(n):
    charset = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    code = ''
    for i in range(4):
        n, remainder = divmod(n, 26)
        code = charset[remainder] + code
    return code

i = 0
for gene in clustered_db:
    gene.id = gene.id[:12] + number_to_code(i) + gene.id[12:]
    gene.description = "" 
    i = i+1

SeqIO.write(clustered_db, output_db, "fasta")