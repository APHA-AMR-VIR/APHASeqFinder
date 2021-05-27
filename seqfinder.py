#!/usr/bin/env python

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
from os import listdir
import numpy as np

###################################
# functions used in this script
def find_file(pattern, path):
    result = []
    for f in os.listdir(path):
        if fnmatch.fnmatch(f, pattern):
            result.append(f) #os.path.join(root, name))
    return result

def run_cmd(lis,ver=1):
    if ver==1:
        print("********************************")
        print(" ".join(lis))
        print("********************************")
    os.system(" ".join(lis))
'''    
def trimming(sample_name,sample_folder,r1_file,r2_file):
    r1_file_trimmed_paired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r1_file_trimmed_unpaired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_unpaired.fastq.gz"))
    r2_file_trimmed_paired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r2_file_trimmed_unpaired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_unpaired.fastq.gz"))
    run_cmd(['java -jar',os.path.join(soft_path,'third_party_software','trimmomatic-0.30.jar'),'PE -phred33','-threads 1',r1_file,r2_file,r1_file_trimmed_paired,r1_file_trimmed_unpaired,r2_file_trimmed_paired,r2_file_trimmed_unpaired,'SLIDINGWINDOW:10:20 MINLEN:80'])
    fastq_trimmed=os.path.join(sample_folder,sample_name+'.trimmed.fastq.gz')
    run_cmd(['cat',r1_file_trimmed_paired,r1_file_trimmed_unpaired,r2_file_trimmed_paired,r2_file_trimmed_unpaired,'>',fastq_trimmed])    
    return(fastq_trimmed)    
'''
def trimming(sample_name,sample_folder,r1_file,r2_file):
    r1_file_trimmed_paired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r2_file_trimmed_paired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    #fastq_trimmed=os.path.join(sample_folder,sample_name+'.trimmed.fastq.gz')
    run_cmd(['java -jar',os.path.join(soft_path,'third_party_software','trimmomatic-0.30.jar'),'SE -phred33',r1_file,r1_file_trimmed_paired,'SLIDINGWINDOW:10:20 MINLEN:80'])
    run_cmd(['java -jar',os.path.join(soft_path,'third_party_software','trimmomatic-0.30.jar'),'SE -phred33',r2_file,r2_file_trimmed_paired,'SLIDINGWINDOW:10:20 MINLEN:80'])
    fastq_trimmed=os.path.join(sample_folder,sample_name+'.trimmed.fastq.gz')
    run_cmd(['cat',r1_file_trimmed_paired,r2_file_trimmed_paired,'>',fastq_trimmed])
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
def readfnaFileSeveralSequences(fname):
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
    return(ids,seqs)

def contig_checker(fname):
    ids,seqs=readfnaFileSeveralSequences(fname)
    lens=[len(seqs[i]) for i in range(len(seqs))]
    return([fname.split(os.sep)[-1],len(seqs),round(np.mean(lens)),min(lens)])

def filter_contigs(fin,fout,min_size=300):
    ids,seqs=readfnaFileSeveralSequences(fin)
    lines=[]
    for i in range(len(ids)):
        if len(seqs[i])>=300:
            if i==0:
                lines.append(">"+ids[i])
            else:
                lines.append("\n"+">"+ids[i])
            for j in range(0,len(seqs[i]),60):
                lines.append("\n"+seqs[i][j:j+60])

    fileOut = open(fout, 'w')
    fileOut.writelines(lines)
    fileOut.close()
    return(contig_checker(fout))



def delete_files(results_path,strin):
    print('Deleting files ending in '+strin+' in '+results_path)
    files_to_delete=[str(path) for path in Path(results_path).rglob('*'+strin)]
    run_cmd(['rm']+files_to_delete,ver=0)
    #for fil in files_to_delete:
    #    print('Deleting file: '+fil)
    #    run_cmd(['rm',fil])
        
   
def one_sample(r1_file):
    try:
        r2_file=r1_file.replace("R1","R2")
        sample_name=r1_file.split(os.sep)[-1].split("_")[0]
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
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Shu',os.path.join(sample_folder,sample_name+'.60.sam'),'>',os.path.join(sample_folder,sample_name+'.bam')]) 
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'sort',os.path.join(sample_folder,sample_name+'.bam'),os.path.join(sample_folder,sample_name+'.sorted')])       
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'index',os.path.join(sample_folder,sample_name+'.sorted.bam')])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'faidx',reference])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'mpileup','-q -1 -uf',reference,os.path.join(sample_folder,sample_name+'.sorted.bam'),'>',os.path.join(sample_folder,sample_name+'.mpileup.bcf')])    
        run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -cg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.mpileup.vcf')])
        run_cmd(['perl',os.path.join(soft_path,'third_party_software','vcfutils.pl'),'vcf2fq',os.path.join(sample_folder,sample_name+'.mpileup.vcf'),'>',os.path.join(sample_folder,sample_name+'.fq')])   
        run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -vcg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.snp')])
        snps_filter(150,1,2,os.path.join(sample_folder,sample_name+'.snp'),os.path.join(sample_folder,sample_name+'_SN.csv'))
        pileup_stats_file=pileupStats(os.path.join(sample_folder,sample_name+".mpileup.vcf"),os.path.join(sample_folder,sample_name+".fq"),reference,3)
        run_cmd(['python',os.path.join(soft_path,'calculatePresentGenes3.py'),os.path.join(sample_folder,sample_name+'.fq'),reference,pileup_stats_file,'150','0.2','2','yes'])
        
        compare_file=find_file('*_CompareTo_*',sample_folder)[0]
        run_cmd(['python',os.path.join(soft_path,'good_snps_12_2019.py'),os.path.join(sample_folder,compare_file),str(percentageID),str(numofsnps),efsa_dict])
    except:
        print("Something went wrong with the mapping stage.")
        print("Therefore, "+sample_name+" no processed at all.")
        return([sample_name,"NA","NA","NA"])
        
        
    ##### assembler
    if fastas_folder=="":
        assem_folder=os.path.join(sample_folder,sample_name+"_spades")
        k='-k 25'
        run_cmd(['spades.py','--careful',k,'-1',r1_file_trimmed_paired,'-2',r2_file_trimmed_paired,'-o',assem_folder])
        assem_file=os.path.join(assem_folder,"contigs.fasta")
        fasta_file=os.path.join(sample_folder,sample_name+"_spades.fasta")
        if os.path.isfile(assem_file):
            [fasta_name,number_contigs,mean_size_bp,min_size_bp]=filter_contigs(assem_file,fasta_file)
            spades_stats_table=[["fasta_name","number_contigs","mean_size_bp","min_size_bp"],[fasta_name,number_contigs,mean_size_bp,min_size_bp]]
            writeCSV(fasta_file.replace("_spades.fasta","_spades_stats.csv"),spades_stats_table)
        else:
            print("Something went wrong with the spades assembly.")
            print("Therefore, no abricate run")
            return([sample_name,"NA","NA","NA"])
    else:
        fasta_file=[f for f in listdir(fastas_folder) if f[:len(sample_name)]==sample_name and f.split(".")[-1] in ["fasta","fa"]]
        if len(fasta_file)==1:
            fasta_file=os.path.join(fastas_folder,fasta_file[0])
        else:
            print("No corresponding fasta file detected in folder "+fastas_folder)
            print("Therefore, no abricate run")
            return([sample_name,"NA","NA","NA"])
    
    ###### abricate
    #abricate --setupdb
    reference_name=reference.split(os.sep)[-1].replace(".fna","")
    abricate_file_name=os.path.join(sample_folder,sample_name+".abricate")
    run_cmd(['abricate','--datadir',str(Path.home()),'--db',reference_name,fasta_file,'>',abricate_file_name])
    
    if os.path.isfile(abricate_file_name):
        ####### combination seqfinder abricate
        seqfinder_file_name=os.path.join(sample_folder,compare_file.replace('.csv','_good_snps.csv'))
        run_cmd(['python',os.path.join(soft_path,'abricate_combine_with_seqfinder_v1.py'),abricate_file_name,seqfinder_file_name])
                    
    ########## deleting unwanted files
    delete_files(sample_folder,'.fastq.gz')
    delete_files(sample_folder,'.sam')
    delete_files(sample_folder,'.bam')
    delete_files(sample_folder,'.bai')
    delete_files(sample_folder,'.sma')
    delete_files(sample_folder,'.smi')
    delete_files(sample_folder,'.bcf')
    delete_files(sample_folder,'.vcf')
    delete_files(sample_folder,'.snp')
    delete_files(sample_folder,'_alignment_stats.csv')
    run_cmd(['rm','-r',assem_folder])               

    
###################################
###################################
###################################
###################################

########### Global variables
percentageID=70
numofsnps=100

###########Command line arguments parshing
soft_path=""
reference=""
data_path=""
R1_pattern=""
results_path=""
efsa_dict=""
fastas_folder=""
ncores=0

args=sys.argv
if len(args)>1:
    arguments_file=args[1]

else:
    arguments_file="/home/javi/APHASeqFinder/template_arguments_file.args"
    print("Argument file not given or doesn't exist. Please re run with /full/path/to/arguments/file/arguments_file.args")
    #sys.exit()

########Loading other arguments from the arguments file
print('Reading arguments from file: '+arguments_file)
with open(arguments_file,'r') as f:
    for line in f:
        if line[0] not in ['\n',' ','#']:
            print('Loading argument: '+ line.strip())
            exec(line.strip())
            
            
''' Explanation of the loaded arguments 
soft_path="/home/user/APHASeqFinder"
reference="/home/user/APHASeqFinder/references/AMR/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fna"
data_path="/home/user/WGS_Data/Project_1/fastq"
results_path="/home/user/WGS_Results/Project_1"
efsa_dict="/home/user/APHASeqFinder/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"
fastas_folder=""
R1_pattern
'''

####### ask user how many cores to use
#try:
#    if ncores==0:
#        value=int(input('Number of cores to use? (There are '+str(mp.cpu_count())+' available:\n'))
#        if value<=mp.cpu_count():
#            ncores=value
#except:
#    ncores=1
 
####### results directory
if not os.path.exists(results_path):
    run_cmd(['mkdir','-p',results_path])
run_cmd(['cp ',arguments_file,os.path.join(results_path,arguments_file.split(os.sep)[-1])])

########## abricate database creation
reference_name=reference.split(os.sep)[-1].split(".")[0]
abricate_ref_folder=os.path.join(Path.home(),reference_name)
if not os.path.exists(abricate_ref_folder):
    run_cmd(['mkdir','-p',abricate_ref_folder])
run_cmd(['cp',reference,os.path.join(abricate_ref_folder,"sequences")])
run_cmd(['makeblastdb','-in',os.path.join(abricate_ref_folder,"sequences"),'-dbtype','nucl','-hash_index'])


###################################################         
fils_R1=find_file("*"+R1_pattern+"*.fastq.gz", data_path)
fils_R2=find_file("*"+R1_pattern.replace("R1","R2")+"*.fastq.gz", data_path)
fils=[os.path.join(data_path,fil) for fil in fils_R1 if fil.replace("R1","R2") in fils_R2]
print("Processing "+str(len(fils))+" samples. These are:")
for fil in fils:
    print(fil.split(os.sep)[-1].split("R1")[0][:-1])


###############sequencial
#for fil in fils[:1]:
#    one_sample(fil)

############### parallel processing of the samples
pool=mp.Pool(ncores)
result=pool.map_async(one_sample,fils)
result.wait()

#### combine results
print('Now combining all the good_snps.csv files')
list_of_files=[p for p in Path(results_path).rglob('*CompareTo*good_snps.csv')]
df = pd.concat((pd.read_csv(f) for f in list_of_files))
df_final=df.drop(['real-len','other'],1)
df_final.to_csv(os.path.join(results_path,results_path.split(os.sep)[-1]+"__"+reference_name+"__seqfinder_compilation.csv"),index=False)
print('Done! Output written to seqfinder_compilation.csv')

print('Now combining all the good_snps_chr.csv files')
list_of_files=[p for p in Path(results_path).rglob('*CompareTo*good_snps_only_chromosomal.csv')]
df = pd.concat((pd.read_csv(f) for f in list_of_files))
df_final=df.drop(['real-len','other'],1)
df_final.to_csv(os.path.join(results_path,results_path.split(os.sep)[-1]+"__"+reference_name+"__seqfinder_chr_compilation.csv"),index=False)
print('Done! Output written to seqfinder_chr_compilation.csv')

run_cmd(['python',os.path.join(soft_path,'compilation_for_abricate_seqfinder.py'),results_path,reference_name])

########## deleting abricate database
run_cmd(['rm','-r',abricate_ref_folder])

    
