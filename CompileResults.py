import os,sys,re,csv,itertools
import os, fnmatch

patho="/media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMRGenes/Results"


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

files = find('*_ARGEnesTable.csv', patho)

genes =readTable(files[0],",")
genes = [[x[0],x[1]] for x in genes[1:]]
strains = ["Gene","Length"]

for fil in files:
    strains.append(fil.split(os.sep)[-1].split("_")[0])
    table = readTable(fil,",")
    tableNames=[x[0] for x in table]
    for gene in genes:
        if gene[0] in tableNames:
            ind = tableNames.index(gene[0])
            gene.append(float(table[ind][-1])*100)
        else:
            gene.append("NA")

genes = [strains] + genes
writeCSV(os.path.join(patho,"AMRGenes_Compilation.csv"),genes)
        

    