#python2.6 $HOME/Desktop/software/duprm10.py $HOME/Desktop/data/try/AFT-16-01273-13_1.fastq $HOME/Desktop/data/try/AFT-16-01273-13_2.fastq $HOME/Desktop/data/try/A1.fastq $HOME/Desktop/data/try/A2.fastq $HOME/Desktop/data/try/A3.fastq $HOME/Desktop/data/try/A4.fastq
from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord

import sys

print "Removing exact duplicates"
#start = time.time()
args=sys.argv
Unique_seqs=[]
outfile = open(args[2],"w")
outfileR = open(args[3],"w")
fastq_iter = SeqIO.parse(open(args[1]),"fastq")

for rec in fastq_iter:
    if str(rec.seq) not in Unique_seqs:
        Unique_seqs.append(str(rec.seq))
        SeqIO.write([rec],outfile,"fastq")
    else :
        SeqIO.write([rec],outfileR,"fastq")
#        ##

#fastq_iter = SeqIO.parse(open(args[1]),"fastq")
#for rec in fastq_iter:
#    if str(rec.id) in Unique_ids:
#        SeqIO.write([rec],outfile,"fastq")


outfile.close()
outfileR.close()

