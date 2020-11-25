args=sys.argv

fqFile=args[1] #"/media/Second3TB/Work/WorkVLA/Data/WGS_Analysis/results/WGS/ECO-MLST-MG1655-alleles/Mapped/CVL101/CVL101.fq"
fastaFile=fqFile[:-2]+".fasta"


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
    return ids,seqs

ids,seqs=readfqFileSeveralContigs(fqFile)

fasta=open(fastaFile,"w")
for idi,seq in zip(ids,seqs):
    fasta.write(">"+idi+"\n")
    for i in range(0,len(seq),60):
        fasta.write(seq[i:i+60]+"\n")
fasta.close()
    


            
    


