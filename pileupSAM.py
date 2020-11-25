import os, sys

args=sys.argv

if len(args)<2:
    fname ="/media/First3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results/fastqs_2/AMR_fasta_sequence_small/Mapped/CVL100/CVL100.sorted.sam" #"/media/First3TB/Work/WorkVLA/Data/Virology/Denise/Mapped/RV1124/RV1124.sorted.sam"
    reffile = "/media/First3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results/fastqs_2/AMR_fasta_sequence_small/Mapped/CVL100/AMR_fasta_sequence_small.fna" #"/media/First3TB/Work/WorkVLA/Data/Virology/Denise/Mapped/RV1124/RV1124.fna"
    mapqth = 0
else:    
    fname =sys.argv[1]
    reffile=sys.argv[2]
    mapqth=int(sys.argv[3])  


revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

def readRef(fname):
    fileIn = open(fname, 'rb')
    dataOut = fileIn.readlines()
    dataOut=[x.rstrip() for x in dataOut]
    dataOut="".join(dataOut[1:])
    fileIn.close()
    return dataOut
    
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
    
def writeListToTextFile(lista,fileName):
    f = open(fileName, 'w')
    f.write("\n".join(lista))
    f.close()
    print "file "+fileName+" saved."
def isNumber(ch):
    if ch in "0123456789":
        return 1
    else:
        return 0
def indChar(strin):
     stin = [isNumber(x) for x in strin]
     return [i for i, x in enumerate(stin) if x == 0]
def cigar(strin):
    inds = indChar(strin)    
    return [strin[:inds[0]]]+[strin[inds[i]+1:inds[i+1]] for i in range(len(inds)-1)],[strin[inds[i]] for i in range(len(inds)-1)]+[strin[-1]]
    
ids,seqs=readfnaFileSeveralSequences(reffile)

pileup = [[[ids[j],str(i+1),seqs[j][i],str(0),"","",""] for i in range(len(seqs[j]))] for j in range(len(ids))]

fileIn = open(fname, 'r')
cont1=0
for line in fileIn:
    if not line[0]=="@":
        cont1=cont1+1
        if "betaL-g0446_CTX-M-58" in line:
            print line
fileIn.close()

cont=0
with open(fname, "r") as sam:
    for line in sam:
        cont=cont+1
        if cont%1000000==0:
            print str(cont)+" of "+str(cont1)+" rows done" 
        if not line[0]=="@":
            line = line.strip()
            line=line.split("\t")
            flag=int(line[1])
            mapq=int(line[4])
            pos=int(line[3])
            refInd = ids.index(line[2])
            if pos>0 and mapq>=mapqth: #line[3]!="0"
                seq = line[9].upper()
                qual =line[10]
                cig = line[5]
                cigaro = cigar(cig)
                char="."
                if flag & 16:
                    #seq=revcompl(seq)
                    char=","
                #for i in range(len(cigaro[0])):
                #    print cigaro
                seqPos=0
                refPos=0
                for i in range(len(cigaro[0])):
                    if cigaro[1][i]=="M":
                        lenMatch=int(cigaro[0][i])
                        for p in range(lenMatch):
                            if pileup[refInd][pos-1+refPos+p][2]==seq[seqPos+p]:  
                                pileup[refInd][pos-1+refPos+p][4]=pileup[refInd][pos-1+refPos+p][4]+char
                            else:
                                pileup[refInd][pos-1+refPos+p][4]=pileup[refInd][pos-1+refPos+p][4]+seq[seqPos+p]
                            pileup[refInd][pos-1+refPos+p][5]=pileup[refInd][pos-1+refPos+p][5]+qual[seqPos+p]
                    if cigaro[1][i]=="I":
                        seqPos=seqPos+int(cigaro[0][i])
                    if cigaro[1][i] in ["S","H"]:
                        seqPos=seqPos+int(cigaro[0][i])
                        #refPos=refPos+int(cigaro[0][i])
                    if cigaro[1][i] in ["M"]:
                        seqPos=seqPos+int(cigaro[0][i])
                        refPos=refPos+int(cigaro[0][i])
                    if cigaro[1][i] in ["D","N","P"]:
                        refPos=refPos+int(cigaro[0][i])
                
sam.close()
pileup = [item for sublist in pileup for item in sublist]
for i in range(len(pileup)): pileup[i][3]=str(len(pileup[i][4]))
pileup=["\t".join(x) for x in pileup]
writeListToTextFile(pileup,fname+"_mapqth_"+str(mapqth)+".pileup")

            