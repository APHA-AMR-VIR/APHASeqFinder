import os,sys,re,csv,itertools
from operator import *
from itertools import *

args=sys.argv

if len(args)<2:
    fqFile="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results_20150612/20s/20-0B-Contigs/20-7/20-7.fq"
    fnaFile="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMR_pipe/references/20-0B-Contigs.fna"
    csvFile="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results_20150612/20s/20-0B-Contigs/20-7/20-7_alignment_stats.csv"
else:    
    fqFile=sys.argv[1]
    fnaFile=sys.argv[2]
    
def readTable(fname,ch):
    infile=open(fname,"rb")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return dataOut
    
def readfqFileSeveralContigsOld(fname):
    fileIn = open(fname, 'rb')
    lines= fileIn.readlines()
    idsPos=[]
    plusPos=[]
    cont=0
    for line in lines:
        if line[0]=="@":
            idsPos.append(cont)
        if line=="+\n":
            plusPos.append(cont)
        cont=cont+1
    ids=[]
    seqs=[]
    for idPos,pluPos in zip(idsPos,plusPos):
        ids.append(lines[idPos][1:-1])
        seq=lines[idPos+1:pluPos]
        seq="".join([x[:-1] for x in seq])
        seqs.append(seq)
    return ids,seqs
    
def readfqFileSeveralSequences(fname):
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
    return ids,seqs

def readfnaFileSeveralSequences(fname):
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
            seq=seq+line[:-1].strip()
    seqs.append(seq)
    return ids,seqs

def readfnaFile(fname):
    fileIn = open(fname, 'rb')
    dataOut = fileIn.readlines()
    dataOut=[x.rstrip() for x in dataOut[1:]]
    dataOut="".join(dataOut)
    fileIn.close()
    return dataOut

def findDeletionPos(seq,ch):
    gaps = [m.start() for m in re.finditer(ch.upper(), seq.upper())]
    intervals = []
    for k, g in groupby(enumerate(gaps), lambda (i,x):i-x):
        interval = map(itemgetter(1), g)
        intervals= intervals + [[interval[0],interval[-1]]]
    return intervals

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."


ids,seqs=readfqFileSeveralSequences(fqFile)

refids,refseqs = readfnaFileSeveralSequences(fnaFile)

cvs = readTable(csvFile,",")
prefN=cvs[0].index("refName")
ppos=cvs[0].index("pos")
pref=cvs[0].index("ref")
palt=cvs[0].index("alt_vcf")
pQual=cvs[0].index("qual")

table=[]
#ids = [x[0] for x in cvs[1:]]
tabledel=[["Gene","Star","End","Length","Sequence"]]

print len(ids)
print len(list(set(ids)))

for idref,seqref in zip(refids,refseqs):
    if idref in ids:
        seq = seqs[ids.index(idref)]
        Ns=seq.upper().count("N")
        if Ns>0:
            intervals=findDeletionPos(seq,"N")
            for interval in intervals:
                tabledel.append([idref,interval[0],interval[1],interval[1]-interval[0]+1,seqref[interval[0]:interval[1]]])
    else:
        seq = ""
        Ns = 0
    per = round(float((len(seq)-Ns)/float(len(seqref))),2)
    table.append([idref,len(seqref),len(seq),Ns,per])
   
table.sort(key=lambda x: x[0]) #-x[4])
table=[["id","real-len","mapped-len","Ns","% mapped-Ns/real"]]+table

writeCSV(fqFile[:-3]+"_ARGEnesTable.csv",table)
writeCSV(fqFile[:-3]+"_DeletedSequences.csv",tabledel)




































