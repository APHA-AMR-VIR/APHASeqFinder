import os,sys,re,csv,itertools
from operator import *
from itertools import *

args=sys.argv

if len(args)<2:
    fqFile="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results/fastqs/AMR_fasta_sequence/Mapped/ARD120B/ARD120B.fq"
    fnaFile="/media/First3TB/Work/WorkVLA/Projects/NSOR1041/pipelineFASTARef/software/references/AMR_fasta_sequence.fna"
    csvFile="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results/fastqs/AMR_fasta_sequence/Mapped/ARD120B/ARD120B_alignment_stats.csv"
    thqual = 150
    thSnpCovProp=0.2
    thMinGC=2
else:    
    fqFile=sys.argv[1]
    fnaFile=sys.argv[2]
    csvFile=sys.argv[3]
    thqual = int(sys.argv[4])
    thSnpCovProp = float(sys.argv[5])
    thMinGC= int(sys.argv[6])
    
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

def checkSnps(qual,pgcovRF,gcovRR,pgcovAF,pgcovAR):
    if qual>=thqual and pgcovRF<=thSnpCovProp*pgcovAF and gcovRR<=thSnpCovProp*pgcovAR and (pgcovAF>=thMinGC or pgcovAR>=thMinGC):
        return True
    else:
        return False

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
    return ''.join(new)

def revcom(s):
    return complement(s[::-1])

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
    return proteinsequence,codonInfo

    
def gsAnnot(x,seqref):
    refSeq=seqref
    geneSeq=ref[:(x[ppos]-1)]+x[palt]+ref[x[ppos]:]
    
    
    refPro,refCodon=translate_dna(refSeq,x[ppos]-1)
    genePro,geneCodon=translate_dna(geneSeq,x[ppos]-1)
    if refPro==genePro:
        annot="syn"
    else:
        annot="non"    
    return x+[refCodon,geneCodon,annot]



ids,seqs=readfqFileSeveralSequences(fqFile)

refids,refseqs = readfnaFileSeveralSequences(fnaFile)

cvs = readTable(csvFile,",")
prefName=cvs[0].index("refName")
ppos=cvs[0].index("pos")
pref=cvs[0].index("ref")
palt=cvs[0].index("alt_vcf")
pqual=cvs[0].index("qual")

pgcovRF=stats[0].index("gcovRF")
pgcovRR=stats[0].index("gcovRR")
pgcovAF=stats[0].index("gcovAF")
pgcovAR=stats[0].index("gcovAR")

table=[]
ids = [x[0] for x in cvs[1:]]
ids = list(set(ids))

for idref,seqref in zip(refids,refseqs):
    lenRef = len(seqref)
    if idref in ids:
        geneSeq= seqs[ids.index(idref)]
        cov=[x in geneSeq if x!="n"]
        Ns=[x in geneSeq if x=="n"]+(lenRef-len(geneSeq))
        other=[x in geneSeq if (x not isACGT(x) and x!="n")]
        
        mapp=[x for x in cvs if x[prefName]==idref]
        
        snps=[x for x in mapp if x[pref]!=x[palt] and x[palt].upper() in ["A","C","G","T"]]
        goodsnps = [x for x in snps if checkSnps(x[pqual],x[pgcovRF],x[pgcovRR],x[pgcovAF],x[pgcovAR])]
        goodsnpsAnnot=[gsAnnot(x,seqref) for x in goodsnps]
        gsnpssummary=";".join([str(x[ppos])+"("+str(x[ppos]-min(loc))+")"+"-"+x[pref]+"-"+x[palt]+"-"+str(x[-3:]) for x in goodsnpsAnnot])
        if len(cov)!=len(mapp):
            print "****************************************************************"
            print len(cov)
            print len(mapp)
        
    else:
        seq = ""
        #Ns = ""
        eq = ""
        snps=""
        mapp=""
    per = round(float(len(mapp))/lenRef,2)
    table.append([idref,lenRef,len(cov),len(Ns),round(100*float(len(cov))/lenRef,2),len(snps),len(other),len(goodsnps),gsnpssummary])

table.sort(key=lambda x: x[0]) #-x[4])
table=[["id","real-len","mapped-len","Gaps","% mapped-gaps/real","snps","other","goodSnps","annot"]]+table

writeCSV(fqFile[:-3]+"_ARGEnesTable.csv",table)
































