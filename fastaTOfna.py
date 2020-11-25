import os,sys,csv

args=sys.argv

if len(args)<2:
    faFile="/media/First3TB/Work/WorkVLA/Projects/NSOR1041/pipelineFASTARef/references/try.fa"

else:    
    faFile=sys.argv[1]

def writeCSV(fname,matrix):
    with open(fname, "wb") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print "file "+fname+" saved."

nameNoExt=".".join(faFile.split(os.sep)[-1].split(".")[0:-1])
fnaFile=".".join(faFile.split(".")[:-1])+".fna" 
anoFile=".".join(faFile.split(".")[:-1])+".ann" 

inFile = open(faFile,"r")
outFile = open(fnaFile,"w")
othFile = open(anoFile,"w")


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
        
    