import os,sys,re,csv,itertools
import os, fnmatch

args=sys.argv

if len(args)<2:
    patho="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMR_Results/SurreyWGS/NCTC11168_geneList_genes"
    word="_geneList_genes.csv" #"NC_012225_geneList" "6-1_versus_NC_012225_geneList"
    col="goodSnps" #"% mapped-gaps/real" #"#non" #"#syn" #"goodSnps" #"% mapped-gaps/real"
else:
    patho=sys.argv[1]
    word=sys.argv[2]
    col=sys.argv[3]
    runName=sys.argv[4]
    

print patho
print word
print col

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
    
def readTable(fname,ch):
    infile=open(fname,"rb")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return dataOut

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."
        

files = find('*'+word+'*', patho)

first =readTable(files[0],",")

to=2

genes = [x[0:to] for x in first[1:]]
strains = first[0][0:to]


for fil in files:
    print fil
    strains.append(fil.split(os.sep)[-1].split("_")[0])
    table = readTable(fil,",")
    tableNames=[x[0:to] for x in table]
    for gene in genes:
        if gene[0:to] in tableNames:
            ind = tableNames.index(gene[0:to])
            toadd=table[ind][table[0].index(col)]
            gene.append(table[ind][table[0].index(col)])
        else:
            print "Not foound "+str(gene)
            gene.append("NA")#

genes = [strains] + genes
if col =="% mapped-gaps/real": col="PercentageCovered"
if col =="% Cov": col="PercCov"
if col =="good snps": col="SNPs"
if col =="#syn": col="SynonimousSNPs"
if col =="#non": col="NonSynonimousSNPs"
writeCSV(os.path.join(patho,runName+"_Comp_"+col+"_Run_"+patho.split(os.sep)[-2]+"_VS_"+patho.split(os.sep)[-1]+".csv"),genes)
        

    
