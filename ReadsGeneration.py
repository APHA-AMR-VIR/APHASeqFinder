#python C:\MRSA\ReadsGeneration.py C:\MRSA\fnas2\
from os import listdir
import os,sys,os.path

def makereads(dir,fname,readlen,nreads,gap):
	import os
	import random
	import Bio
	from Bio.Seq import Seq
	fastfile1=open(os.path.join(dir,fname[0:-4]+"_1.fastq"),"w")
	fastfile2=open(os.path.join(dir,fname[0:-4]+"_2.fastq"),"w")
	genome=file(os.path.join(dir,fname)).readlines()
	genome.pop(0)
	genome=''.join(genome)
	genome=genome.replace("\n","")
	gLen = len(genome)
	print fname+" is "+str(gLen)+" bases"
	max=0
	for i in range(1,nreads):
		rnum=random.randrange(0,gLen-2*readlen-gap)
		fastfile1.write("@"+str(i)+"_"+str(rnum)+"\n")
		fastfile2.write("@"+str(i)+"_"+str(rnum)+"\n")
		fastfile1.write(genome[rnum:rnum+readlen]+"\n")
		read2=Seq(genome[rnum+readlen+gap:rnum+readlen+gap+readlen])
		fastfile2.write(str(read2.reverse_complement())+"\n")		
		fastfile1.write("+"+"\n")
		fastfile2.write("+"+"\n")	
		fastfile1.write("I"*readlen+"\n")
		fastfile2.write("I"*readlen+"\n")		
		if max<rnum:
			max=rnum
	fastfile1.close()
	fastfile2.close()
	print "Max= "+str(max+2*readlen+gap)
	print "END"

rlen=150
nreads=2000000
gap=0
args= sys.argv
pathO=sys.argv[1]
fnames= os.listdir(pathO)
fnames = [x for x in fnames if "fna" in x]
for fname in fnames:
  print fname
  makereads(pathO,fname,rlen,nreads,gap)

  
		