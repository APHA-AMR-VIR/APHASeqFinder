import os,sys,csv

fnaFile="/media/Second3TB/Work/WorkVLA/Software/NCGM_1984.fna"
glistFile="/media/Second3TB/Work/WorkVLA/Software/WGS_GeneMapper/references/NCGM_1984_geneList.csv"



def readTable(fname,ch):
    infile=open(fname,"rb")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return dataOut




annot = [["Id","start","end","length"]]
ident = ">"+nameNoExt+"\n"
line = inFile.readline()
seq=""
start=1
Id=line[1:-1]

for line in inFile:
    if line[0]==">":
        end=len(seq)
        annot.append([Id,start,end,end-start+1])
        seq = seq + "O"*300
        Id = line[1:-1]
        start = len(seq)+1
    else:
        seq = seq + line[0:-1]
end=len(seq)
annot.append([Id,start,end,end-start+1])
inFile.close()

writeCSV(anoFile,annot)

rows=[seq[x:x+100] for x in range(0,len(seq),100)]
outFile.write(ident+"\n".join(rows))
outFile.close()





#outFile.close()
        
    