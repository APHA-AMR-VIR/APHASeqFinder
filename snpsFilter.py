# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
# This script select snps from a vcf format file (arg4) that verifies
# a test for quality based on the DP4 measument being the reference number of 
# fragments forward (reverse) smaller than arg2 percentage of the alternative
# number of forward (reverse) fragments. The minium coverage (DP parameter) must be
# greater than arg3 and the quality of the consensus call has to be greater than arg1. 
# The output is a csv format file with the snps that verify the conditions metioned.


import csv, sys

args=sys.argv
thqual=int(sys.argv[1]) #200
thprop=float(sys.argv[2]) #0.05
thmin=int(sys.argv[3]) #2
fnameI=sys.argv[4]
fnameO=sys.argv[5]
fOut=open(fnameO,"w")
writer = csv.writer(fOut, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
with open(fnameI) as infile:
	for line in infile:
		if line[0]!="#" and "INDEL" not in line:
			line=line.split()
			refGenome=line[0]
			pos=int(line[1])
			ref=line[3]
			alt=line[4]
			qual=float(line[5])
			det=line[7].split(";")
			cov=int([s for s in det if "DP=" in s][0].split("DP=")[1])
			gcov=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
			if len(alt)>1:
				print line
			if len(alt)==1 and qual>=thqual and gcov[0]<=thprop*gcov[2] and gcov[1]<=thprop*gcov[3] and (gcov[2]>thmin or gcov[3]>thmin):
				writer.writerow([refGenome,pos,ref,alt,qual,cov]+gcov)
fOut.close()
#C:\Python27\python snpsFilter.Py 200 0.05 2 try.snp try_SN.csv
