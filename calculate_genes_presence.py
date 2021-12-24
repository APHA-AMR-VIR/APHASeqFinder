#!/usr/bin/env python

'''
APHASeqfinder
version 4.0.0
submitted to github on 23/12/2021
Javier Nunez, AMR Team, Bacteriology
Animal and Plant Health Agency
This script resolve the absence or presence of the database genes
'''


import sys,re,csv,os
from operator import itemgetter
from itertools import groupby
import numpy

def read_csv(fname,ch,dig=0):
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

def write_csv(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")
    
def read_fq(fname):
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

def read_fasta(fname):
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

def pileup_vcf(vcf_file_name,fq_file_name,ref_name):
    print("Procesing vcf file: "+vcf_file_name)
    ref = "".join(read_fasta(ref_name)[1])
    reflen = len(ref)
    print("reflen "+str(reflen))
    
    fqlen=len("".join(read_fq(fq_file_name)[1]))
    print("fqlen "+str(fqlen))
    
    fileIn = open(vcf_file_name, 'r')
    vcflen=0
    for line in fileIn:
        if line[0]!="#" and "INDEL" not in line:
            vcflen=vcflen+1
    fileIn.close()    
    
    header =["refName","pos","ref","alt_vcf","cov","qual","gcovRF","gcovRR","gcovAF","gcovAR"]
    stats =[header]+vcflen*[["NA",0,"NA","NA",0,0,0,0,0,0]]    
    
    fileIn = open(vcf_file_name, 'r')
    cont=1
    for line in fileIn:
        if line[0]!="#" and "INDEL" not in line:
            if cont%100000==0:
                print(str(cont)+" of "+str(vcflen)+" positions done") 
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
                gcov=[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(",")
                gcov=[int(x) for x in gcov]
            else:
                gcov=[0,0,0,0]
            #else:
            #    gcov=[0,0,0,0]
            stats[cont] = [line[0],pos,ref_vcf,alt_vcf,cov,qual]+gcov
            cont=cont+1
    return(stats)

################################################################################################
#run_cmd(['python',os.path.join(soft_path,'calculatePresentGenes3.py'),sample_folder,ref_copy,mlst_fasta,'150','0.2','2','yes'])
args=sys.argv

if len(args)<2:
    sample_folder="/home/javi/APHASeqFinder_new_version/try_data/with_mlst/101-288"
    ref_fasta="/home/javi/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fasta"
    mlst_fasta="/home/javi/APHASeqFinder_new_version/references/mlst/ECO-MLST-MG1655-alleles.fna"
    thqual = 150
    thSnpCovProp=0.2
    thMinGC=2
else:    
    sample_folder=sys.argv[1]
    ref_fasta=sys.argv[2]
    mlst_fasta=sys.argv[3]
    thqual = int(sys.argv[4])
    thSnpCovProp = float(sys.argv[5])
    thMinGC= int(sys.argv[6])


vcf_file_name=os.path.join(sample_folder,sample_folder.split(os.sep)[-1]+'.mpileup.vcf')
fq_file_name=os.path.join(sample_folder,sample_folder.split(os.sep)[-1]+'.fq')
ids,seqs=read_fq(fq_file_name)

refids,refseqs=read_fasta(ref_fasta)

mlst_ids,mlst_seqs=read_fasta(mlst_fasta)



cvs_table=pileup_vcf(vcf_file_name,fq_file_name,ref_fasta)
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
        
        len_Ns_not_extremes=len([x for x in geneSeq if x=="n"])
        lenNs=len_Ns_not_extremes+(lenRef-len(geneSeq))
        other=[x for x in geneSeq if (not isACGT(x) and x!="n")]
        
        mapp_gene=[x for x in cvs_table if x[prefName]==idref]
        covAver=round(numpy.mean([x[pcov] for x in mapp_gene if x[pcov]>0]),2)
        snps=[x for x in mapp_gene if x[pref]!=x[palt] and x[palt].upper() in ["A","C","G","T"]]
        goodsnps = [x for x in snps if checkSnps(x[pqual],x[pgcovRF],x[pgcovRR],x[pgcovAF],x[pgcovAR])]
        allgoodsnps=allgoodsnps+goodsnps
        goodsnpsAnnot=[gsAnnot(x,seqref) for x in goodsnps]
        gsnpssummary=";".join([str(x[ppos])+"-"+x[pref]+"-"+x[palt]+"-"+str(x[-3:]) for x in goodsnpsAnnot])
        syn=gsnpssummary.count("syn")
        non=gsnpssummary.count("non")
        
        ###### deletions inside the gene 
        if len_Ns_not_extremes>0:
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
        goodsnps=""
        gsnpssummary=""
        syn=0
        non=0
    table.append([idref,lenRef,len(cov),covAver,lenNs,round(100*float(len(cov))/lenRef,2),len(snps),len(other),len(goodsnps),syn,non,gsnpssummary])


#### mlst normalization
mlst_tab=[x[3] for x in table if x[0] in mlst_ids]
print("Normalising by the mlst median coverage")
median_cov_ref_mlst=numpy.median(mlst_tab)
   
for i in range(len(table)):
    norm_value=round(table[i][3]/float(median_cov_ref_mlst),2)
    table[i]=table[i][:4]+[norm_value]+table[i][4:]

mlst_tab=[x for x in table if x[0] in mlst_ids]
table=[x for x in table if x[0] not in mlst_ids]
########################
    
table.sort(key=lambda x: x[0]) #-x[4])
table=[["id","ref_len","mapped_len","mean_depth","norm_depth","non_calls","perc_mapped","snps","other","good_snps","syn","non","annotation"]]+table
out_file=fq_file_name[:-3]+"_CompareTo_"+ref_fasta.split(os.sep)[-1].split(".")[0]+".csv"
write_csv(out_file,table)

mlst_tab=[["id","ref_len","mapped_len","mean_depth","norm_depth","non_calls","perc_mapped","snps","other","good_snps","syn","non","annotation"]]+mlst_tab
out_file_mlst=fq_file_name[:-3]+"_mlst_only.csv"
write_csv(out_file_mlst,mlst_tab)

write_csv(fq_file_name[:-3]+"_DeletedSequences.csv",tabledel)
write_csv(fq_file_name[:-3]+"_goodsnps.csv",allgoodsnps)
































