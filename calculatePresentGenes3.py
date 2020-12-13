# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
# This script calculate the percentage of the genes present from a database 
# (arg2 in fasta format) in a sample (arg1 and arg3, fq and stats file respectively). It will 
# also provide the mean coverage for each gene and the number of good snps and their
# corresponding annotation. 

import sys,re,csv
from operator import itemgetter
from itertools import groupby
import numpy

args=sys.argv

if len(args)<2:
    fqFile="/media/Second3TB/Work/WorkVLA/Data/MRSA/Results/WGS_Results/APHA/GeneList_20151014/MRSA01/MRSA01.fq"
    fnaFile="/media/Second3TB/Work/WorkVLA/Data/MRSA/references/GeneList_20151014.fas"
    csvFile="/media/Second3TB/Work/WorkVLA/Data/MRSA/Results/WGS_Results/APHA/GeneList_20151014/MRSA01/MRSA01_alignment_stats.csv"
    thqual = 150
    thSnpCovProp=0.2
    thMinGC=2
    gSAnnot="yes"
else:    
    fqFile=sys.argv[1]
    fnaFile=sys.argv[2]
    csvFile=sys.argv[3]
    thqual = int(sys.argv[4])
    thSnpCovProp = float(sys.argv[5])
    thMinGC= int(sys.argv[6])
    gSAnnot=sys.argv[7]
    
def readTable(fname,ch,dig=0):
    infile=open(fname,"r")
    data = csv.reader(infile, delimiter=ch)
    if dig!=0:
        dataOut = [[strToint(x) for x in row] for row in data]
    else:
        dataOut = [row for row in data]
    infile.close()
    return(dataOut)
    
def strToint(s):
    if "." in s:
        return(float(s))
    if s.isdigit():
        return(int(s))
    return(s)
    
def readfqFileSeveralContigsOld(fname):
    fileIn = open(fname, 'r')
    lines= fileIn.readlines()
    lines=[x.strip() for x in lines]
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
        ids.append(lines[idPos][1:])
        seq=lines[idPos+1:pluPos]
        seq="".join([x for x in seq])
        seqs.append(seq)
    return(ids,seqs)
    
def readfqFileSeveralSequences(fname):
    fileIn = open(fname, 'r')
    lines= fileIn.readlines()
    lines=[x.strip() for x in lines]
    pluses = [i for i in range(0,len(lines)) if lines[i]=="+"]
    ini=0
    ids=[]
    seqs=[]
    for plus in pluses:
        ids=ids+[lines[ini][1:]]
        seqs=seqs+["".join([lines[i] for i in range(ini+1,plus)])]
        ini = plus+(plus-ini)
    return(ids,seqs)

def readfnaFileSeveralSequences(fname):
    fileIn = open(fname, 'r')
    lines= fileIn.readlines()
    lines=[x.strip() for x in lines]
    ids=[]
    seqs=[]
    seq=""
    for line in lines:
        #print("**************")
        if line[0]==">":
            if len(ids)>0:
                seqs.append(seq)
                seq=""
            ids.append(line[1:])
        else:
            seq=seq+line
    seqs.append(seq)
    return(ids,seqs)

def findDeletionPos(seq,ch):
    gaps = [m.start() for m in re.finditer(ch.upper(), seq.upper())]
    intervals = []
    for k,g in groupby(enumerate(gaps), lambda x: x[1]-x[0]):
        interval = list(map(itemgetter(1), g))
        intervals= intervals + [[interval[0],interval[-1]]]
    return(intervals)

def writeCSV(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")

def checkSnps(qual,gcovRF,gcovRR,gcovAF,gcovAR):
    if qual>=thqual and gcovRF<=thSnpCovProp*gcovAF and gcovRR<=thSnpCovProp*gcovAR and (gcovAF>=thMinGC or gcovAR>=thMinGC):
        return(True)
    else:
        return(False)

def isACGT(c):
    return c in ["a","A","c","C","g","G","t","T"]
    
def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'} 
    letters = list(s)
    new=[]
    for l in letters:
        if isACGT(l):
            new=new+[basecomplement[l]]
        else:
            new=new+[l]
    return(''.join(new))

def revcom(s):
    return(complement(s[::-1]))

def translate_dna(sequence,pos=-1):

    codontable = {
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    #proteinsequence = ''
    #start = sequence.find('ATG')
    #sequencestart = sequence[int(start):]
    #stop = sequencestart.find('TAA')
    #cds = str(sequencestart[:int(stop)+3])

    #for n in range(0,len(cds),3):
    #    if cds[n:n+3] in codontable:
    #        proteinsequence += codontable[cds[n:n+3]]
    #        print proteinsequence
    #    sequence = ''
    codonInfo="NA"
    proteinsequence = ''
    for n in range(0,len(sequence),3):
        if sequence[n:n+3] in codontable:
            proteinsequence += codontable[sequence[n:n+3]]
            if pos>-1 and pos in range(n,n+3):
                codonInfo=[sequence[n:pos]+"("+sequence[pos]+")"+sequence[pos+1:n+3],codontable[sequence[n:n+3]]]
        #else:
        #    print "not found "+sequence[n:n+3]+" at position "+str(n)
    return(proteinsequence,codonInfo)

    
def gsAnnot(x,seqref):
    refSeq=seqref
    geneSeq=seqref[:(x[ppos]-1)]+x[palt]+seqref[x[ppos]:]
    
    
    refPro,refCodon=translate_dna(refSeq,x[ppos]-1)
    genePro,geneCodon=translate_dna(geneSeq,x[ppos]-1)
    if refPro==genePro:
        annot="syn"
    else:
        annot="non"    
    return(x+[refCodon,geneCodon,annot])

################################################################################################

ids,seqs=readfqFileSeveralSequences(fqFile)

refids,refseqs = readfnaFileSeveralSequences(fnaFile)

cvs_table=readTable(csvFile,",",1)

prefName=cvs_table[0].index("refName")
ppos=cvs_table[0].index("pos")
pref=cvs_table[0].index("ref")
palt=cvs_table[0].index("alt_vcf")
pqual=cvs_table[0].index("qual")
pcov=cvs_table[0].index("cov")
pgcovRF=cvs_table[0].index("gcovRF")
pgcovRR=cvs_table[0].index("gcovRR")
pgcovAF=cvs_table[0].index("gcovAF")
pgcovAR=cvs_table[0].index("gcovAR")

table=[]
tabledel=[["Gene","Star","End","Length","Sequence","SequenceToBlast"]]
allgoodsnps=[]
for idref,seqref in zip(refids,refseqs):
    lenRef = len(seqref)
    
    if idref in ids:
 
        geneSeq=seqs[ids.index(idref)]
        
        cov=[x for x in geneSeq if x!="n"]
        
        Ns=[x for x in geneSeq if x=="n"]
        lenNs=len(Ns)+(lenRef-len(geneSeq))
        other=[x for x in geneSeq if (not isACGT(x) and x!="n")]
        
        mapp=[x for x in cvs_table if x[prefName]==idref]
        covVector = [x[pcov] for x in mapp]
        covAver=round(numpy.mean(covVector),2)
        snps=[x for x in mapp if x[pref]!=x[palt] and x[palt].upper() in ["A","C","G","T"]]
        if gSAnnot.upper()=="YES":
            goodsnps = [x for x in snps if checkSnps(x[pqual],x[pgcovRF],x[pgcovRR],x[pgcovAF],x[pgcovAR])]
            allgoodsnps=allgoodsnps+goodsnps
            goodsnpsAnnot=[gsAnnot(x,seqref) for x in goodsnps]
            gsnpssummary=";".join([str(x[ppos])+"-"+x[pref]+"-"+x[palt]+"-"+str(x[-3:]) for x in goodsnpsAnnot])
            syn=gsnpssummary.count("syn")
            non=gsnpssummary.count("non")
        
        if len(Ns)>0:
            intervals=findDeletionPos(geneSeq,"N")
            for interval in intervals:
                if interval[0]<20:
                    seqStart=0
                else:
                    seqStart=interval[0]-20
                tabledel.append([idref,interval[0]+1,interval[1]+1,interval[1]-interval[0]+1,seqref[interval[0]:interval[1]+1],seqref[seqStart:interval[1]+1]])
    else:
        cov=""
        lenNs=0
        other=""
        snps=""
        covAver=0
        if gSAnnot.upper()=="YES":
            goodsnps=""
            gsnpssummary=""
            syn=0
            non=0
    if gSAnnot.upper()=="YES":
        table.append([idref,lenRef,len(cov),covAver,lenNs,round(100*float(len(cov))/lenRef,2),len(snps),len(other),len(goodsnps),syn,non,gsnpssummary])
    else:
        table.append([idref,lenRef,len(cov),covAver,lenNs,round(100*float(len(cov))/lenRef,2),len(snps),len(other)])
table.sort(key=lambda x: x[0]) #-x[4])
if gSAnnot.upper()=="YES":
    table=[["id","real-len","mapped-len","meanCov","Gaps","% mapped-gaps/real","snps","other","goodSnps","#syn","#non","annot"]]+table
else:
    table=[["id","real-len","mapped-len","meanCov","Gaps","% mapped-gaps/real","snps","other"]]+table

out_file=fqFile[:-3]+"_CompareTo_"+fnaFile.split("/")[-1].split(".")[0]+".csv"
writeCSV(out_file,table)
writeCSV(fqFile[:-3]+"_DeletedSequences.csv",tabledel)
writeCSV(fqFile[:-3]+"_goodsnps.csv",allgoodsnps)
































