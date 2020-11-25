# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
# This script collects information on the vcf file which is recorded in the document:
# args6_Mapped_args2_Info.txt. Also a light version of the vcf file (in csv format)
# is generated for downstream analysis. 
#


import os, csv,  sys, os.path, math,numpy
from Bio import SeqIO


def mean(data):
    return float(sum(data))/len(data)

def std(data):
    n = len(data)
    c = mean(data)
    ss = sum((x-c)**2 for x in data)	
    return math.sqrt(ss/(n-1))

def countReadsFastq(fname):
    num_lines = sum(1 for line in open(fname))
    return num_lines/4
 
def countUniqueReadsInFastq(fname):
    outfile=open(fname[:-5]+"unique.fastq","w")
    #done_ids=set()
    done_seqs=[]
    fastq_iter=SeqIO.parse(open(fname),"fastq")
    cont=0
    contT=0
    for rec in fastq_iter:
        contT=contT+1
        #if str(rec.id) not in done_ids:
        if str(rec.seq) not in done_seqs: 
           done_seqs=done_seqs+[str(rec.seq)]
           SeqIO.write([rec],outfile,"fastq")
           cont=cont+1
    outfile.close()
    cmd = "gzip "+ fname[:-5]+"unique.fastq"
    os.system(cmd)           
    return contT,cont
    
def medianReadLength(fname):
    fastq_iter = SeqIO.parse(open(fname),"fastq")
    lens=[]
    for rec in fastq_iter:
        lens.append(len(str(rec.seq)))
    return numpy.median(lens)

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."
def readfqFileSeveralContigs(fname):
    fileIn = open(fname, 'rb')
    lines= fileIn.readlines()
    pluses = [i for i in range(0,len(lines)) if lines[i]=="+\n"]
    ini=0
    ids=[]
    seqs=[]
    for plus in pluses:
        ids.append(lines[ini][1:-1])
        seqs.append("".join([lines[i][:-1] for i in range(ini+1,plus)]))
        ini = plus+(plus-ini)
    print ids
    return ids,seqs

def readRefSeveralContigs(fname):
    fileIn = open(fname, 'rb')
    lines= fileIn.readlines()
    ids=[]
    seqs=[]
    seq=""
    for line in lines:
        if line[0]==">":
            if len(ids)>0:
                seqs.append(seq)
                seq=""
            ids.append(line[1:-1].strip())
        else:
            seq=seq+line[:-1]
    seqs.append(seq)
    return "".join(seqs)

def pileupStats(fname,fqfname,refname,rn):
    print "Procesing vcf file: "+fname
    ref = readRefSeveralContigs(refname)
    reflen = len(ref)
    print "reflen "+str(reflen)
    
    fqlen=len("".join(readfqFileSeveralContigs(fqfname)[1]))
    print "fqlen "+str(fqlen)
    
    fileIn = open(fname, 'r')
    vcflen=0
    for line in fileIn:
        if line[0]!="#" and "INDEL" not in line:
            vcflen=vcflen+1
    fileIn.close()    
    
    header =["refName","pos","ref","alt_vcf","cov","qual","gcovRF","gcovRR","gcovAF","gcovAR"]
    stats =[header]+vcflen*[["NA",0,"NA","NA",0,0,0,0,0,0]]    
    
    fileIn = open(fname, 'r')
    indels =0
    cont=1
    for line in fileIn:
        if line[0]!="#" and "INDEL" not in line:
            if cont%100000==0:
                print str(cont)+" of "+str(vcflen)+" positions done" 
            line=line.split()
            pos=int(line[1])
            ref_vcf=line[3]
            if line[4]==".":
                alt_vcf=line[3]
            else:
                alt_vcf=line[4]
            qual=float(line[5])
            det=line[7].split(";")
            cov=int([s for s in det if "DP=" in s][0].split("DP=")[1])
            if "DP4=" in line[7]:
                gcov=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
            else:
                gcov=[0,0,0,0]
            #else:
            #    gcov=[0,0,0,0]
            stats[cont] = [line[0],pos,ref_vcf,alt_vcf,cov,qual]+gcov
            cont=cont+1
        else:
            if line[0]!="#" and "INDEL" in line:
                indels=indels+1
    writeCSV(fname.replace('.pileup.vcf','_alignment_stats.csv'),stats)
    



##############################################
##############################################
########## Alignment stats
#arg[1] isolate new name $nf
#arg[2] run name $7
#arg[3] reference genome $2 reference genome
#arg[4] stats file name $4 stats file name
#arg[5] isolate original name $of
#arg[6] isolate directory $isoldir
#arg[7] results path $5
#Rscript $1/myStats.r $nf $3 $2 $4 $of $isoldir $5

rn=3
args=sys.argv
#f1 = countReadsFastq(os.path.join(args[6],args[1]+"_R1.fastq"))
#f1med=medianReadLength(os.path.join(args[6],args[1]+"_R1.fastq"))

#f2 = countReadsFastq(os.path.join(args[6],args[1]+"_R2.fastq"))
#f2med=medianReadLength(os.path.join(args[6],args[1]+"_R2.fastq"))

#fT=f1+f2


#ff1 = countReadsFastq(os.path.join(args[6],args[1]+".nodup_R1.fastq"))
#ff1med=medianReadLength(os.path.join(args[6],args[1]+".nodup_R1.fastq"))

#ff2 = countReadsFastq(os.path.join(args[6],args[1]+".nodup_R2.fastq"))
#ff2med=medianReadLength(os.path.join(args[6],args[1]+".nodup_R2.fastq"))

#ffT=ff1+ff2

#fff1 = countReadsFastq(os.path.join(args[6],args[1]+".TrimIn_R1.fastq"))
#fff1med=medianReadLength(os.path.join(args[6],args[1]+".TrimIn_R2.fastq"))

#fff2 = countReadsFastq(os.path.join(args[6],args[1]+".TrimIn_R1.fastq"))
#fff2med=medianReadLength(os.path.join(args[6],args[1]+".TrimIn_R2.fastq"))

#fffT=fff1+fff2


#fmap= countReadsFastq(os.path.join(args[6],args[1]+".sorted.fastq"))
#fmapmed= medianReadLength(os.path.join(args[6],args[1]+".sorted.fastq"))

#permap=round(100*(float(fmap)/fffT),rn)

alignStats = pileupStats(os.path.join(args[6],args[1]+".pileup.vcf"),args[1]+".fq",args[3],rn)

#results = [args[1],args[2],f1,f1med,f2,f2med,fT, ff1,ff1med,ff2,ff2med, ffT, fff1,fff1med,fff2,fff2med, fffT,fmap,permap,fmapmed]
#results = [args[1],args[2],ff1,ff1med,ff2,ff2med, ffT, fff1,fff1med,fff2,fff2med, fffT,fmap,permap,fmapmed]
#results = [args[1],args[2],fff1,fff1med,fff2,fff2med,fffT,fmap,permap,fmapmed]
#results = [args[1],args[2],fmap,fmapmed]
#results = [str(x) for x in results]
#fileIn = open(os.path.join(args[6],"Mapped_"+args[2]+"_Info.txt"), 'a')
#fileIn.write("\t".join(map(str,results))+"\n")
#fileIn.close()

#fileIn = open(os.path.join(args[7],args[4]), 'a')
#fileIn.write("\t".join(map(str,results))+"\n")
#fileIn.close()


#start = time.clock()	
#print doSam2("/home/javier/Desktop/try/AF2122-NCBI.sam")	
#end = time.clock()
#print end - start
#res=pileupStats("/home/javier/Desktop/try/Africanum-NCBI.pileup.vcf","/home/javier/Desktop/try/AF2122.fna","/home/javier/Desktop/try/Africanum-NCBI.fq",3)
#print res
