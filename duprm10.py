# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
# This script filter out sequence fragments from both mates fastq files if 
# that are repeated for both mates at the same time, keeping only one copy
#
# Usage:
# python2.6 $HOME/Desktop/software/duprm10.py $HOME/Desktop/data/try/AFT-16-01273-13_1.fastq $HOME/Desktop/data/try/AFT-16-01273-13_2.fastq $HOME/Desktop/data/try/A1.fastq $HOME/Desktop/data/try/A2.fastq $HOME/Desktop/data/try/A3.fastq $HOME/Desktop/data/try/A4.fastq
# arg1 and arg2 are the input mates fastq files
# arg3 and arg4 are the filter in mates fastq files
# arg3 and arg4 are the filter out mates fastq files

from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
import itertools
import os,sys
import time

print "Removing exact duplicates"
#start = time.time()
args=sys.argv
Unique_seqs1=set()
Unique_ids1=set()
Unique_seqs2=set()
Unique_ids2=set()
outfile1 = open(args[3],"w")
outfile2 = open(args[4],"w")
outfileR1 = open(args[5],"w")
outfileR2 = open(args[6],"w")
fastq_iter1 = SeqIO.parse(open(args[1]),"fastq")
fastq_iter2 = SeqIO.parse(open(args[2]),"fastq")
for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
    if str(rec1.seq) not in Unique_seqs1:
       Unique_seqs1.add(str(rec1.seq))
       Unique_ids1.add(str(rec1.id))
    if str(rec2.seq) not in Unique_seqs2:
       Unique_seqs2.add(str(rec2.seq))
       Unique_ids2.add(str(rec2.id))
fastq_iter1 = SeqIO.parse(open(args[1]),"fastq")
fastq_iter2 = SeqIO.parse(open(args[2]),"fastq") 
for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
    if (str(rec1.id) in Unique_ids1) or (str(rec2.id) in Unique_ids2):
        SeqIO.write([rec1],outfile1,"fastq")
        SeqIO.write([rec2],outfile2,"fastq")
    else :
        SeqIO.write([rec1],outfileR1,"fastq")
        SeqIO.write([rec2],outfileR2,"fastq")
outfile1.close()
outfile2.close()
outfileR1.close()
outfileR2.close()
#end = time.time()
#print end-start

