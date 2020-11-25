import os,sys,re,csv,itertools
import os, fnmatch

args=sys.argv

if len(args)<2:
    patho="/media/Second3TB/Work/WorkVLA/Data/WGS_Analysis/results/PlasmidsDNA/Denovo"
    word=".contigs.fa"
    pathBlast="/media/Second3TB/Work/WorkVLA/Data/WGS_Analysis/software_20150210/blast"
    ref="NC_000913.fna"
else:
    patho=sys.argv[1]
    word=sys.argv[2]
    

print patho
print word

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

files = find('*'+word, patho)
print files
first =readTable(files[0],",")

for fil in files:
    cmd = "blastall"+ " -p blastn -d "+os.path.join(pathBlast,ref)+" -i "+fil + " -o " + fil+".crunch" + " -m 8 -e 0.01"
    print cmd    
    os.system(cmd)
        

    