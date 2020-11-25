import csv,numpy
import os, sys, os.path
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

args=sys.argv
if len(args)<2:
    fastaFile="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/MultiGenes/reference/MGESeqs_20151012.fas"
    coorTable="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/MultiGenes/reference/MGENames_20151012.csv"
    idi="Accession"
    loc="Location"
    newid="Probe Name"

    
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

def translate_dna(sequence):
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
    proteinsequence = ''
    for n in range(0,len(sequence),3):
        if sequence[n:n+3] in codontable:
            proteinsequence += codontable[sequence[n:n+3]]
        else:
            print "not found "+sequence[n:n+3]+" at position "+str(n)
    
    return proteinsequence


def readTable(fname,ch):
    infile=open(fname,"rb")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return(dataOut)
    
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
            ids.append(line[1:].strip()) #.split(" ")[0])
        else:
            seq=seq+line[:-1]
    seqs.append(seq)
    print "Length of ref parts:"
    print map(len,seqs)
    return ids,seqs #"".join(seqs)

idis,seqs=readRefSeveralContigs(fastaFile)

info=readTable(coorTable,",")
pidi=info[0].index(idi)
ploc=info[0].index(loc)
pnewid=info[0].index(newid)

toWriteFNA=""
#toWriteFAA=""

for gene in info[1:]:
    comp=1
    if "(complement)" in gene[ploc]:
        gene[ploc]=gene[ploc].replace("(complement)","")
        comp=-1
    loc=map(int,gene[ploc].split(".."))
    idi=">"+gene[pnewid]+"\n"
    toWriteFNA=toWriteFNA+idi
    fastaSeq=[seqs[i] for i in range(len(idis)) if gene[pidi] in idis[i]]
    if len(fastaSeq)>1:
        print gene[pidi]
        print "hell"
    fastaSeq=fastaSeq[0]    
    seq=fastaSeq[(min(loc)-1):max(loc)]
    if comp==-1:
        print seq
        seq=revcom(seq)
        print seq
    for i in range(0,len(seq),100):
        toWriteFNA=toWriteFNA+seq[i:i+100]+"\n"
    
    #pseq=translate_dna(seq)
    #toWriteFAA=toWriteFAA+idi
    #for i in range(0,len(pseq),100):
    #    toWriteFAA=toWriteFAA+pseq[i:i+100]+"\n"

FNA=open(coorTable[:-3]+"fna","w")
FNA.write(toWriteFNA)
FNA.close()
#FAA=open(coorTable[:-3]+"faa","w")
#FAA.write(toWriteFAA)
#FAA.close()
