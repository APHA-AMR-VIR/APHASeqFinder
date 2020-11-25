# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
# This script top the quality mapping value to the givin arg3 for the
# sam file provided as arg1.A new sam file with name arg2 is generated. 


import sys

args=sys.argv

if len(args)<2:
    fname ="/media/First3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results/fastqs_2/AMR_fasta_sequence_small/Mapped/CVL100/CVL100.sam" #"/media/First3TB/Work/WorkVLA/Data/Virology/Denise/Mapped/RV1124/RV1124.sorted.sam"
    fnameOut ="/media/First3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results/fastqs_2/AMR_fasta_sequence_small/Mapped/CVL100/CVL100.F4C.sam"    
    mapqth = 60
else:    
    fname =sys.argv[1]
    fnameOut =sys.argv[2]
    mapqth=int(sys.argv[3])



print "Processing sam file: "+fname
fileIn = open(fname, 'r')
fileOut = open(fnameOut,'w')
for line in fileIn:
    if not line[0]=="@":
        line = line.strip()
        line = line.split("\t")
        if int(line[4])<mapqth:
            line[4]=str(mapqth) 
        line = "\t".join(line)+"\n"
    fileOut.write(line)
fileIn.close()
fileOut.close()
