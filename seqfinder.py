'''
Dec 2020
Javier Nunez, AMR Team, Bacteriology
Animal and Plant Health Agency
WGS_mapper tool
This software map paired end fastq files against a reference providing:
    - QC summary table
    - Consensus and VCF files
    - Deletions and insertions annotation
    - SNP calling
    
To run this from a terminal, cd into the folder that contains this file. Then:

python WGS_mapper.py path/WGS_mapper_arguments_bTB.args runName

example:

python WGS_mapper.py /media/javier/3TB/Software/WGS_BTB/WGS_mapper_arguments_bTB.args 20180601_SB4020__10087
 
'''

import sys,os,fnmatch,csv,os.path
import multiprocessing as mp
from pathlib import Path
import pandas as pd

###################################
# functions used in this script
def find_file(pattern, path):
    result = []
    for f in os.listdir(path):
        if fnmatch.fnmatch(f, pattern):
            result.append(f) #os.path.join(root, name))
    return result

def run_cmd(lis):
    print("********************************")
    print(" ".join(lis))
    print("********************************")
    os.system(" ".join(lis))
    
def trimming(sample_name,sample_folder,r1_file,r2_file):
    r1_file_trimmed_paired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r1_file_trimmed_unpaired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_unpaired.fastq.gz"))
    r2_file_trimmed_paired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r2_file_trimmed_unpaired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_unpaired.fastq.gz"))
    fastq_trimmed=os.path.join(sample_folder,sample_name+'.trimmed.fastq.gz')
    run_cmd(['java -jar',os.path.join(soft_path,'third_party_software','trimmomatic-0.39.jar'),'PE','-threads 1',r1_file,r2_file,r1_file_trimmed_paired,r1_file_trimmed_unpaired,r2_file_trimmed_paired,r2_file_trimmed_unpaired,'SLIDINGWINDOW:10:20 MINLEN:80'])
    run_cmd(['cat',r1_file_trimmed_paired,r1_file_trimmed_unpaired,r2_file_trimmed_paired,r2_file_trimmed_unpaired,'>',fastq_trimmed])
    return(fastq_trimmed)    

def sam_flag_changed(fname,mapqth=60):
    print("Changing flags to sam file: "+fname)
    fnameOut=fname.replace(".F4.sam",".60.sam")
    fileIn = open(fname, 'r')
    fileOut = open(fnameOut,'w')
    for line in fileIn:
        if not line[0]=="@":
            line = line.strip()
            line = line.split("\t")
            if int(line[4])<mapqth:
                line[4]=str(mapqth) 
            line = "\t".join(line)+"\n"
        fileOut.write(line)
    fileIn.close()
    fileOut.close()
    
    
def snps_filter(thqual,thprop,thmin,fnameI,fnameO):
    fOut=open(fnameO,"w")
    writer = csv.writer(fOut, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
    with open(fnameI) as infile:
        for line in infile:
            if line[0]!="#" and "INDEL" not in line:
                line=line.split()
                refGenome=line[0]
                pos=int(line[1])
                ref=line[3]
                alt=line[4]
                qual=float(line[5])
                det=line[7].split(";")
                cov=int([s for s in det if "DP=" in s][0].split("DP=")[1])
                gcov=[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(",")
                gcov=[int(x) for x in gcov]
                if len(alt)>1:
                    print(line)
                if len(alt)==1 and qual>=thqual and gcov[0]<=thprop*gcov[2] and gcov[1]<=thprop*gcov[3] and (gcov[2]>thmin or gcov[3]>thmin):
                    writer.writerow([refGenome,pos,ref,alt,qual,cov]+gcov)
    fOut.close()

def writeCSV(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")
def readfqFileSeveralContigs(fname):
    print('readfqFileSeveralSequences')
    print(fname)
    fileIn = open(fname, 'r')
    lines= fileIn.readlines()
    lines=[x.strip() for x in lines]
    pluses = [i for i in range(0,len(lines)) if lines[i]=="+"]
    ini=0
    ids=[]
    seqs=[]
    for plus in pluses:
        ids=ids+[lines[ini][1:-1]]
        seqs=seqs+["".join([lines[i] for i in range(ini+1,plus)])]
        ini = plus+(plus-ini)
    return(ids,seqs)

def readRefSeveralContigs(fname):
    fileIn = open(fname, 'r')
    lines= fileIn.readlines()
    lines=[x.strip() for x in lines]
    ids=[]
    seqs=[]
    seq=""
    for line in lines:
        if line[0]==">":
            if len(ids)>0:
                seqs.append(seq)
                seq=""
            ids.append(line[1:])
        else:
            seq=seq+line
    seqs.append(seq)
    return "".join(seqs)

def pileupStats(fname,fqfname,refname,rn):
    print("Procesing vcf file: "+fname)
    ref = readRefSeveralContigs(refname)
    reflen = len(ref)
    print("reflen "+str(reflen))
    
    fqlen=len("".join(readfqFileSeveralContigs(fqfname)[1]))
    print("fqlen "+str(fqlen))
    
    fileIn = open(fname, 'r')
    vcflen=0
    for line in fileIn:
        if line[0]!="#" and "INDEL" not in line:
            vcflen=vcflen+1
    fileIn.close()    
    
    header =["refName","pos","ref","alt_vcf","cov","qual","gcovRF","gcovRR","gcovAF","gcovAR"]
    stats =[header]+vcflen*[["NA",0,"NA","NA",0,0,0,0,0,0]]    
    
    fileIn = open(fname, 'r')
    indels =0
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
        else:
            if line[0]!="#" and "INDEL" in line:
                indels=indels+1
    out_file=fname.replace('.mpileup.vcf','_alignment_stats.csv')
    writeCSV(out_file,stats)
    return(out_file)
    
def one_sample(r1_file):
    r2_file=r1_file.replace("R1","R2")
    sample_name=r1_file.split(os.sep)[-1].split("R1")[0][:-1]
    print("Processing sample: "+sample_name)
    sample_folder=os.path.join(results_path,sample_name)
    if not os.path.exists(sample_folder):
        os.system('mkdir -p '+sample_folder)
    fastq_trimmed=trimming(sample_name,sample_folder,r1_file,r2_file)  
    
    ref_index=os.path.join(sample_folder,reference.split(os.sep)[-1])
    run_cmd([os.path.join(soft_path,'third_party_software','smalt'),'index -k 13 -s 6',ref_index,reference])
    run_cmd([os.path.join(soft_path,'third_party_software','smalt'),'map -d -1 -y 0.7 -x -f samsoft -o',os.path.join(sample_folder,sample_name+".sam"),ref_index, fastq_trimmed])
    run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Sh -F 4',os.path.join(sample_folder,sample_name+".sam"),'>',os.path.join(sample_folder,sample_name+".F4.sam")])
    sam_flag_changed(os.path.join(sample_folder,sample_name+".F4.sam")) #produces .60.sam
    run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Shu',os.path.join(sample_folder,sample_name+'.F4.sam'),'>',os.path.join(sample_folder,sample_name+'.bam')]) 
    run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'sort',os.path.join(sample_folder,sample_name+'.bam'),os.path.join(sample_folder,sample_name+'.sorted')])       
    run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'index',os.path.join(sample_folder,sample_name+'.sorted.bam')])
    run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'faidx -o',os.path.join(sample_folder,reference.split(os.sep)[-1]+'.faidx'),reference])
    run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'mpileup','-uf',reference,os.path.join(sample_folder,sample_name+'.sorted.bam'),'>',os.path.join(sample_folder,sample_name+'.mpileup.bcf')])    
    run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -cg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.mpileup.vcf')])
    run_cmd(['perl',os.path.join(soft_path,'third_party_software','vcfutils.pl'),'vcf2fq',os.path.join(sample_folder,sample_name+'.mpileup.vcf'),'>',os.path.join(sample_folder,sample_name+'.fq')])   
    run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -vcg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.snp')])
    snps_filter(150,1,2,os.path.join(sample_folder,sample_name+'.snp'),os.path.join(sample_folder,sample_name+'_SN.csv'))
    pileup_stats_file=pileupStats(os.path.join(sample_folder,sample_name+".mpileup.vcf"),os.path.join(sample_folder,sample_name+".fq"),reference,3)
    run_cmd(['python',os.path.join(soft_path,'calculatePresentGenes3.py'),os.path.join(sample_folder,sample_name+'.fq'),reference,pileup_stats_file,'150','0.2','2','yes'])

    compare_file=find_file('*_CompareTo_*',sample_folder)[0]
    print(compare_file)
    print(efsa_dict)
    run_cmd(['python',os.path.join(soft_path,'good_snps_12_2019.py'),os.path.join(sample_folder,compare_file),str(percentageID),str(numofsnps),efsa_dict])
    
  
    
###################################
###################################
###################################
###################################

#Command line arguments parshing
#Two arguments expected when calling WGS_mapper.py
soft_path=""
reference=""
data_path=""
results_path=""
percentageID=0
numofsnps=0
efsa_dict=""

args=sys.argv
if len(args)>1:
    arguments_file=args[1]

else:
    #Just to run the script from a GUI
    arguments_file='/home/javi/APHASeqFinder/template_arguments_file.args'

#Loading other arguments from the arguments file
print('Reading arguments from file: '+arguments_file)
with open(arguments_file,'r') as f:
    for line in f:
        if line[0] not in ['\n',' ','#']:
            print('Loading argument: '+ line.strip())
            exec(line.strip())
'''  
soft_path="/home/javi/APHASeqFinder"
reference="/home/javi/APHASeqFinder/references/AMR/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fna"
data_path="/home/javi/WGS_Data/Project_1/fastq"
results_path="/home/javi/WGS_Results/Project_1"
percentageID=70
numofsnps=5
efsa_dict="/home/javi/APHASeqFinder/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"
'''

try:
    value=int(input('Number of cores to use? (There are '+str(mp.cpu_count())+' available:\n'))
    if value<=mp.cpu_count():
        ncores=value
except:
    ncores=1
    

if not os.path.exists(results_path):
    os.system('mkdir -p '+results_path)

fils=find_file("*R1*.fastq.gz", data_path)
fils=[os.path.join(data_path,fil) for fil in fils]

#for fil in fils:
#    one_sample(fil)

pool=mp.Pool(ncores)
result=pool.map_async(one_sample,fils)
result.wait()

print('Now combining all the good_snps.csv files')
list_of_files=[p for p in Path(results_path).rglob('*CompareTo*good_snps.csv')]
df = pd.concat((pd.read_csv(f) for f in list_of_files))
df_final=df.drop(['real-len','other'],1)
df_final.to_csv(os.path.join(results_path,"seqfinder_compilation.csv"),index=False)
print('Done! Output written to seqfinder_compilation.csv')

print('Now combining all the good_snps_chr.csv files')
list_of_files=[p for p in Path(results_path).rglob('*CompareTo*good_snps_only_chromosomal.csv')]
df = pd.concat((pd.read_csv(f) for f in list_of_files))
df_final=df.drop(['real-len','other'],1)
df_final.to_csv(os.path.join(results_path,"seqfinder_chr_compilation.csv"),index=False)
print('Done! Output written to seqfinder_chr_compilation.csv')
    




