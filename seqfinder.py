#!/usr/bin/env python

'''
APHASeqfinder
version 4.0.6
submitted to github on 04/10/2022
from pathlib import Path
from os import listdir
from datetime import datetime
Javier Nunez, AMR Team, Bacteriology
Animal and Plant Health Agency
'''

### packages
import sys,os,fnmatch,csv,os.path
from multiprocessing import Pool,current_process
from datetime import datetime
from pathlib import Path
from os import listdir
import pandas as pd
import re
import time

try:
    import numpy as np
    print("Python module 'numpy' is installed, continuing\n")
except ModuleNotFoundError:
    print("Python module 'numpy' is not installed.\nPlease install numpy before running this pipeline")
    print("\nTo install numpy: sudo pip3 install numpy or conda install -c anaconda numpy")
    exit()
try:
    import pandas
    print("Python module 'pandas' is installed, continuing\n")
except ModuleNotFoundError:
    print("Python module 'pandas' is not installed.\nPlease install pandas before running this pipeline")
    print("\nTo install pandas: sudo pip3 install pandas or conda install -c anaconda pandas")
    exit()

###################################
# functions used in this script

def write_csv(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")

def read_csv(fname,ch=','):
    infile=open(fname,"r")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return(dataOut)


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


def trimming(sample_name,sample_folder,r1_file,r2_file):
    r1_file_trimmed_paired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r2_file_trimmed_paired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    unpaired_trimmed=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_unpaired.fastq.gz"))
    run_cmd([os.path.join(soft_path,'third_party_software','fastp'),' -i '+r1_file+' -I '+r2_file+' -o '+r1_file_trimmed_paired+' -O '+r2_file_trimmed_paired+' --unpaired1 '+unpaired_trimmed+' --unpaired2 '+unpaired_trimmed+' -l 80 -r --cut_right_window_size 10 -w 1 -j /dev/null -h /dev/null'])
    fastq_trimmed=os.path.join(sample_folder,sample_name+'.trimmed.fastq.gz')
    run_cmd(['cat',r1_file_trimmed_paired,r2_file_trimmed_paired,unpaired_trimmed,'>',fastq_trimmed])
    run_cmd(['rm',r1_file_trimmed_paired])
    run_cmd(['rm',r2_file_trimmed_paired])
    run_cmd(['rm',unpaired_trimmed])
    return(fastq_trimmed)    


def sam_flag_changed(file_name_in,mapqth=60):
    print("Changing flags to sam file: "+file_name_in)
    file_name_out=file_name_in.replace(".F4.sam",".60.sam")
    file_in = open(file_name_in, 'r')
    file_out = open(file_name_out,'w')
    for line in file_in:
        if not line[0]=="@":
            line = line.strip()
            line = line.split("\t")
            if int(line[4])<mapqth:
                line[4]=str(mapqth) 
            line = "\t".join(line)+"\n"
        file_out.write(line)
    file_in.close()
    file_out.close()
    
    
def snps_filter(th_qual,th_prop,th_min,file_name_in,file_name_out):
    file_out=open(file_name_out,"w")
    writer = csv.writer(file_out, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
    with open(file_name_in) as infile:
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
                if len(alt)==1 and qual>=th_qual and gcov[0]<=th_prop*gcov[2] and gcov[1]<=th_prop*gcov[3] and (gcov[2]>th_min or gcov[3]>th_min):
                    writer.writerow([refGenome,pos,ref,alt,qual,cov]+gcov)
    file_out.close()



def read_fq(file_name_in):
    file_in = open(file_name_in, 'r')
    lines= file_in.readlines()
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


def read_fasta(fname):
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
    ids,seqs=read_fasta(fname)
    lens=[len(seqs[i]) for i in range(len(seqs))]
    return([fname.split(os.sep)[-1],len(seqs),round(np.mean(lens)),min(lens)])

def filter_contigs(fin,fout,min_size=300):
    ids,seqs=read_fasta(fin)
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
'''
def combine_tables_of_results(mother_path,pattern,out_file):
    print('Now combining files in '+mother_path+' that contains '+pattern)
    list_of_files=[str(p) for p in Path(output_path).rglob(pattern)]
    out_table=read_csv(list_of_files[0])
    for fil in list_of_files[1:]:
        out_table=out_table+[x for x in read_csv(fil)[1:]]
    write_csv(out_file,out_table)
'''
def combine_tables_of_results(mother_path, pattern, out_file):
    print('Now combining files in ' + mother_path + ' that contains ' + pattern)
    try:
        list_of_files = [str(p) for p in Path(output_path).rglob(pattern)]
        out_table = read_csv(list_of_files[0])
        for fil in list_of_files[1:]:
            out_table = out_table + [x for x in read_csv(fil)[1:]]
        write_csv(out_file, out_table)
    except IndexError as z:
        print(f"An index error occurred, likely because there are no files to combine: {z}")
        print("\nContinuing with final checks then exiting\n")

def delete_files(output_path,strin):
    print("********************************")
    print('Deleting files ending in '+strin+' in '+output_path)
    files_to_delete=[str(path) for path in Path(output_path).rglob('*'+strin)]
    for fil in files_to_delete:
        print('Deleting file: '+fil)
        run_cmd(['rm',fil])
    print("********************************")
        
   
def one_sample(file_to_process):
    try:
        r1_file=file_to_process[0]
        r2_file=file_to_process[1]
        fasta_file=file_to_process[2]
        sample_name=r1_file.split(os.sep)[-1].split("_")[0]
        current_pid = os.getpid()  # Get the current process's PID
        with open(os.path.join(output_path, "processed.csv"), "a") as f:
            f.write(sample_name + " processed by PID "+str(current_pid)+"\n")
        print("Processing sample: "+sample_name)
        sample_folder=os.path.join(output_path,sample_name)
        if not os.path.exists(sample_folder):
            os.system('mkdir -p '+sample_folder)
        
        fastq_trimmed=trimming(sample_name,sample_folder,r1_file,r2_file)
    
        ref_index=os.path.join(sample_folder,ref_copy.split(os.sep)[-1])
        run_cmd([os.path.join(soft_path,'third_party_software','smalt'),'index -k 13 -s 6',ref_index,ref_copy])
        run_cmd([os.path.join(soft_path,'third_party_software','smalt'),'map -d -1 -y 0.7 -x -f samsoft -o',os.path.join(sample_folder,sample_name+".sam"),ref_index, fastq_trimmed])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Sh -F 4',os.path.join(sample_folder,sample_name+".sam"),'>',os.path.join(sample_folder,sample_name+".F4.sam")])
        sam_flag_changed(os.path.join(sample_folder,sample_name+".F4.sam")) #produces .60.sam
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Shu',os.path.join(sample_folder,sample_name+'.60.sam'),'>',os.path.join(sample_folder,sample_name+'.bam')]) 
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'sort',os.path.join(sample_folder,sample_name+'.bam'),os.path.join(sample_folder,sample_name+'.sorted')])       
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'index',os.path.join(sample_folder,sample_name+'.sorted.bam')])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'faidx',ref_copy])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'mpileup','-q -1 -uf',ref_copy,os.path.join(sample_folder,sample_name+'.sorted.bam'),'>',os.path.join(sample_folder,sample_name+'.mpileup.bcf')])    
        run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -cg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.mpileup.vcf')])
        run_cmd(['perl',os.path.join(soft_path,'third_party_software','vcfutils.pl'),'vcf2fq',os.path.join(sample_folder,sample_name+'.mpileup.vcf'),'>',os.path.join(sample_folder,sample_name+'.fq')])   
        run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -vcg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.snp')])
        snps_filter(th_qual,th_prop,th_min,os.path.join(sample_folder,sample_name+'.snp'),os.path.join(sample_folder,sample_name+'_SN.csv'))
        #pileup_stats_file=pileupStats(os.path.join(sample_folder,sample_name+".mpileup.vcf"),os.path.join(sample_folder,sample_name+".fq"),ref_copy,3)
        run_cmd(['python',os.path.join(soft_path,'calculate_genes_presence.py'),sample_folder,ref_copy,mlst_fasta,str(th_qual),str(th_prop),str(th_min)])
        
        found_compare_file = False

        for filename in os.listdir(sample_folder):
            if "CompareTo" in filename:
                compare_file = find_file('*_CompareTo_*', sample_folder)[0]
                found_compare_file = True
                try:
                    run_cmd(['python', os.path.join(soft_path, 'good_snps_filtering.py'), os.path.join(sample_folder, compare_file), str(percentageID), str(numofsnps), efsa_dict, database_type, vir_dict, reference_name, bio_metal_dict])
                    default_output_file = os.path.join(sample_folder, compare_file)
                except:
                    print("\ngood_snps_filtering.py encountered an error, trying to move onto the next sample\n")
                    raise ValueError("Exiting")
                    
        if not found_compare_file:
            raise ValueError("Exiting")
        
            
    except:
        print("Something went wrong with the mapping stage.")
        print("Therefore, "+sample_name+" wasn't processed at all.")
        os.system("find "+sample_folder+" ! -name mlst_error.txt -type f -exec rm {} +")
        #Making output for files that filed the MLST mapping stage
        mlst_error=os.path.join(sample_folder,"mlst_error.txt")
        if os.path.isfile(mlst_error):
            columns = ['warnings','strain', 'id', 'gene', 'antimicrobial', 'class', 'ref_len', 'mapped_len', 'mean_depth', 'norm_depth', 'non_calls', 'perc_mapped', 'snps', 'other', 'good_snps', 'syn', 'non', 'annotation', 'reference', 'EFSA_dict', 'result_abr', 'contig_abr', 'plasmid_abr', 'coverage_abr', 'identity_abr', 'location_abr']
            failed_isolate_df = pd.DataFrame(columns=dict.fromkeys(columns, []))
            failed_isolate_df.loc[0,'strain']=sample_name
            failed_isolate_df = failed_isolate_df.fillna('No MLST genes map to sample')
            output_filename = sample_name+"_CompareTo_"+reference_name+".csv"
            output_filename_good_snps=output_filename.replace(".csv","_good_snps.csv")
            output_filename_abricate_seqfinder=output_filename.replace('.csv','_good_snps_abricate_seqfinder.csv')
            output_file = os.path.join(sample_folder,output_filename)
            output_file_good_snps = os.path.join(sample_folder,output_filename_good_snps)
            output_file_abricate_seqfinder = os.path.join(sample_folder,output_filename_abricate_seqfinder)
            failed_isolate_df.to_csv(output_file, index=False)
            failed_isolate_df.to_csv(output_file_good_snps, index=False)
            failed_isolate_df.to_csv(output_file_abricate_seqfinder, index=False)
        raise ValueError("Exiting")
    
    
    ###### abricate
    
    abricate_file_name=os.path.join(sample_folder,sample_name+".abricate")
    seqfinder_file_name=find_file('*_good_snps.csv',sample_folder)[0]
    seqfinder_file_name=os.path.join(sample_folder,seqfinder_file_name)
    default_output_file=os.path.join(sample_folder,compare_file)

    with open(seqfinder_file_name,"r") as f:
        contents = f.read()
        if "No genes passed filter" in contents:
            os.system("find {} ! -name '*CompareTo*' -type f -exec rm {} +".format(sample_folder, "{}"))
            failed_seqfinder_file_name=seqfinder_file_name.replace('_good_snps.csv','_good_snps_abricate_seqfinder.csv')
            with open(failed_seqfinder_file_name, "w") as f:
                f.write(contents)
            raise Exception('Exiting function, no genes passed filter in '+sample_name)
            #exit()
        else:
            os.system("abricate --datadir {} --mincov 10 --minid 10 --db {} {} > {}".format(str(Path.home()), reference_name, fasta_file, abricate_file_name))
            abricate_version_command = os.popen('abricate --version').read()
            abricate_version = abricate_version_command.split(" ")[1]
            abricate_version=abricate_version.replace('\n', '')
            
            with open(fasta_file) as f:
                first_line = f.readline().strip()
                first_line = first_line.replace('>','')                            
    #seqfinder_file_name=find_file('*_good_snps.csv',sample_folder)[0]
    #seqfinder_file_name=os.path.join(sample_folder,seqfinder_file_name)
            if os.path.isfile(abricate_file_name) and os.path.isfile(seqfinder_file_name):
                ####### combination seqfinder abricate
                try:
                    run_cmd(['python',os.path.join(soft_path,'abricate_combine_with_seqfinder.py'),abricate_file_name,seqfinder_file_name,default_output_file,efsa_dict,first_line,abricate_version,mlst_fasta])
                except:
                    print('\nabricate_combine_with_seqfinder.py failed. Exiting\n')
                    raise ValueError('exiting')
                    
                    
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

    
###################################
###################################
### Starting
start_time=str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
###################################
###################################

########### Global variables
version='4.0.6'
date='04/08/2023'

## good snps thresholds
th_qual=150
th_prop=0.2
th_min=2
## good_snps filtering stage thresholds 
percentageID=70
numofsnps=100

###########Command line arguments parshing
soft_path=""
reference=""
database_type="AMR"
data_path=""
R1_pattern=""
results_path=""
efsa_dict=""
vir_dict=""
bio_metal_dict=""
fastas_folder=""
ncores=1
sample_list_file=""
mlst_fasta=""
###########Additional variables for folder naming
today= datetime.now()
date_of_run=today.strftime('%Y%m%d')
reference_name=reference.split(os.sep)[-1].split(".")[0]
output_path=os.path.join(results_path,reference_name,date_of_run)


### checking the arguments file
args=sys.argv
if len(args)>1:
    arguments_file=args[1]
    if not os.path.exists(arguments_file):    
        print("Argument file not given or doesn't exist. Please re run with /full/path/to/arguments/file/arguments_file.args")
        sys.exit()
else:
    print("Please provide an arguments file. Check readme.txt for more info")
    sys.exit()

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
database_type="AMR"
data_path="/home/user/WGS_Data/Project_1/fastq"
results_path="/home/user/WGS_Results/Project_1"
efsa_dict="/home/user/APHASeqFinder/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"
fastas_folder=""
R1_pattern
sample_list
'''
 
####### results directory
#######HASHED OUT 07092022 DUE TO FOLDER CREATION ERROR
#if not os.path.exists(output_path):
#    run_cmd(['mkdir','-p',output_path])
#run_cmd(['cp ',arguments_file,os.path.join(output_path,arguments_file.split(os.sep)[-1])])


########## checking that fastas folder exists
if not os.path.exists(fastas_folder):
    print("fastas_folder doesn't exist "+fastas_folder)
    sys.exit()
    
########## checking that paths provided in the template argument exist 
if not os.path.exists(soft_path):
    print('\nIncorrect path to SeqFinder provided in the arguments file, make sure this is correct above.')
    sys.exit()
if not os.path.isfile(reference):
    print('\nReference file cannot be located, make sure you have provided the correct path above.')
    sys.exit()
if not os.path.isfile(mlst_fasta):
    print('\nMLST file cannot be located, make sure you have provided the correct path above.')
    sys.exit()
if not os.path.isfile(efsa_dict):
    print('\nEFSA dictionary cannot be located, make sure you have provided the correct path above.')
    sys.exit()
if not os.path.isfile(vir_dict):
    print('\nVirulence dictionary cannot be located, make sure you have provided the correct path above.')
    sys.exit()
if not os.path.isfile(bio_metal_dict):
    print('\nBiocide and metal resistance dictionary cannot be located, make sure you have provided the correct path above.')
    sys.exit()
 
###### Checking the reference file is formatted correctly
print("\n*****Checking your reference file\n")
error_found = False

with open(reference) as file:
    for line in file:
        if re.search(r"^>[^a-zA-Z0-9]|[^a-zA-Z0-9>\-_.]|>{2,}|^\s*$", line.rstrip('\n')):
            error_found = True
            if re.search(r"^\s*$", line.rstrip('\n')):
                print("There are lines in your reference file that are empty.\n")
            else:
                print(line)
                print("This line in your reference file contains illegal characters or does not have an alphanumeric character directly following '>'.\n")
                print("Only '>', '_', '-', and '.' are permitted in headers, and headers must have an alphanumeric character directly following '>'.\n")
                print("A simple way to fix this for your headers would be: sed -e '/^>/ s/>[^[:alnum:]]/>_/g' -e '/^>/ s/__*/_/g' your_reference_file > your_fixed_reference_file")
            exit()

if not error_found:
    print("*****No errors found in the reference file\n")
    
########## checking that R2 and assembly files exist

fastq_to_process=[]
sample_names = set()  # Store sample names to check for duplicates
R2_pattern=R1_pattern.replace("R1","R2")
today= datetime.now()
date_of_run=today.strftime('%Y%m%d')
reference_name=reference.split(os.sep)[-1].split(".")[0]
output_path=os.path.join(results_path,reference_name,date_of_run)
print("***** Checking samples to be run")
if sample_list_file=="":
    if not os.path.exists(data_path) or not os.path.exists(fastas_folder):
        print("Check your paths arguments data_path and fastas_folder: "+data_path+" and "+fastas_folder)
        sys.exit()
    else:
        summary=[["R1_fastq","R2_status","fasta_status","size"]]
        fastq_R1s=find_file("*"+R1_pattern+"*.fastq.gz", data_path)
        for fastq_R1 in fastq_R1s:
            sample_name = fastq_R1.split(os.sep)[-1].split("_")[0]
            if sample_name in sample_names:
                # Sample with the same name already encountered, ask the user which one to keep
                duplicate_samples = [item for item in fastq_to_process if item[0].split(os.sep)[-1].split("_")[0] == sample_name]
                print(f"A sample with the name '{sample_name}' has been encountered again. "
                      "Do you want to keep the current one (c) or the new one (n)?")
                print("Current sample (c):")
                for duplicate in duplicate_samples:
                    print(f"R1: {duplicate[0]}")
                user_input = input("New sample (n):\nR1: "+data_path+"/"+fastq_R1+"\n")
                
                if user_input.lower() == 'n':
                    # Remove the previous occurrence of the sample
                    fastq_to_process = [item for item in fastq_to_process if item[0].split(os.sep)[-1].split("_")[0] != sample_name]
                    sample_names.remove(sample_name)
                    print(f"Removed the previous occurrence of '{sample_name}'.")
                else:
                    print(f"Keeping the current occurrence of '{sample_name}'.")
            else:
                # New sample name, add it to the set
                sample_names.add(sample_name)
            fastq_R2=fastq_R1.replace(R1_pattern,R2_pattern)
            if os.path.isfile(os.path.join(data_path,fastq_R2)):
                R2_ok="found"
                if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000 or os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:
                    file_size="small"
                    if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000:
                        print(os.path.join(data_path, fastq_R1) + " is less than 10MB. This might not provide an output.")
                    if os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:
                        print(os.path.join(data_path, fastq_R2) + " is less than 10MB. This might not provide an output.")
                else:
                    file_size="good"
                if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000 or os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:   
                    remove_small_files=""
                    while remove_small_files.lower() not in ['y', 'n']:
                        remove_small_files=input("Do you want to remove this isolate from your SeqFinder run (y/n)?\n")
                    if remove_small_files.lower() == "y":
                        R2_ok="too_small_excluded_by_user"
                        print("\nRemoved\n")
                    elif remove_small_files.lower() == "n":
                        print("\nKept\n")   
                        
                
                fasta_file=[f for f in listdir(fastas_folder) if f[:len(sample_name)]==sample_name and f.split(".")[-1] in ["fasta","fa",".fna"]]    
                #fasta_file = [f for f in os.listdir(fastas_folder) if f.startswith(sample_name + "_") and not f[len(sample_name) + 1].isdigit() and f.split(".")[-1] in ["fasta", "fa", ".fna"]]
                if len(fasta_file)==1:
                    fasta_ok="found"
                    fasta_file=fasta_file[0]
                elif len(fasta_file)==0:
                    fasta_ok="not_found"
                    print("Missing assembly file for file: "+fastq_R1)
                else:
                    fasta_ok="several_found"
                    print("More than one assembly for: "+fastq_R1)
            else:
               R2_ok="not_found"
               print("Missing R2 for file: "+fastq_R1)
            
            if R2_ok=="found" and fasta_ok=="found":
                fastq_to_process.append([os.path.join(data_path,fastq_R1),os.path.join(data_path,fastq_R2),os.path.join(fastas_folder,fasta_file),file_size])
                #This sorts the input into the script by the third column so that if small files are chosen to be included by the user they are processed last
                fastq_to_process.sort(key=lambda x: x[3])
                
            summary.append([fastq_R1,R2_ok,fasta_ok,file_size])
else:
    if not os.path.isfile(sample_list_file):
        print("Cannot find sample_list_file at "+sample_list_file)
        sys.exit()
    else:
        summary=[["R1_file","fasta_file","R1_status","R2_status","fasta_status"]]
        samples=read_csv(sample_list_file)
        for sample in samples:
            fastq_R1=sample[0]
            fastq_R2=fastq_R1.replace(R1_pattern,R2_pattern)
            fasta_file=sample[1]
            if os.path.isfile(fastq_R1):
                R1_ok="found"
            else:
                R1_ok="not_found"
                print("Missing file: "+fastq_R1)
            if os.path.isfile(fastq_R2):
                R2_ok="found"
            else:
                R2_ok="not_found"
                print("Missing file: "+fastq_R2)
                
            if os.path.isfile(fasta_file):
                fasta_ok="found"
            else:
                fasta_ok="not_found"
                print("Missing file: "+fasta_file)
            if R1_ok=="found" and R2_ok=="found" and fasta_ok=="found":
                fastq_to_process.append([fastq_R1,fastq_R2,fasta_file])
            
            summary.append(sample+[R1_ok,R2_ok,fasta_ok])

run_cmd(['mkdir -p',output_path])
write_csv(os.path.join(output_path,"initial_file_testing.csv"),summary)
print('********\nOutput will be stored at:'+output_path+'\n********')    
print("Check file "+os.path.join(output_path,"summary.csv")) 
print("***** Processing "+str(len(fastq_to_process))+" samples.")

value=input('happy to go ahead (y/n)?\n')
if value!='y':
    sys.exit()

########## abricate database creation and mlst inclusion
reference_name=reference.split(os.sep)[-1].split(".")[0]
abricate_ref_folder=os.path.join(Path.home(),reference_name)
if not os.path.exists(abricate_ref_folder):
    run_cmd(['mkdir','-p',abricate_ref_folder])

ref_copy=os.path.join(abricate_ref_folder,reference_name+".fasta")   


if os.path.isfile(mlst_fasta):
    print("MLST file found and added")
    run_cmd(['cat',reference,mlst_fasta,'>',ref_copy])

else:
    print("MLST file not found. Please check the MLST file")
    sys.exit()


run_cmd(['cp',ref_copy,os.path.join(abricate_ref_folder,"sequences")])
run_cmd(['makeblastdb','-in',os.path.join(abricate_ref_folder,"sequences"),'-dbtype','nucl','-hash_index'])
'''
############### parallel processing of the samples
pool=Pool(ncores)
try:
    result=pool.map(one_sample,fastq_to_process)
except (ValueError, Exception) as e:
    print(e)
    pool.close()
    pool.join()
'''

# Initialize a dictionary to store the updated information
updated_info_dict = {}

samples_to_reprocess = list(range(len(fastq_to_process)))
pool = Pool(ncores)

try:
    while samples_to_reprocess:
        print("\n***Starting sample processing***\n")
        results = []
        # Create a new list to store samples that should be kept
        remaining_samples = []
        processed_samples = 0  # Define processed_samples here
        
        for i in samples_to_reprocess:
            sample_name = fastq_to_process[i][0].split(
                os.sep)[-1].split("_")[0]
            sample_folder = os.path.join(output_path, sample_name)

            abri_file = os.path.join(sample_folder, sample_name + ".abricate")
            if os.path.isfile(abri_file):
                if os.path.getsize(abri_file):
                    abri = "ok"
                else:
                    abri = "no_contents"
            else:
                abri = "not_present"

            seq_file = os.path.join(
                sample_folder, sample_name + "_CompareTo_" + reference_name + "_good_snps.csv")
            if os.path.isfile(seq_file):
                if os.path.getsize(seq_file):
                    seq = "ok"
                else:
                    seq = "no_contents"
            else:
                seq = "not_present"

            comp_file = os.path.join(sample_folder, sample_name + "_CompareTo_" +
                                     reference_name + "_good_snps_abricate_seqfinder.csv")
            if os.path.isfile(comp_file):
                if os.path.getsize(comp_file):
                    comp = "ok"
                else:
                    comp = "no_contents"
            else:
                comp = "not_present"

            fastq_to_process[i] = fastq_to_process[i] + [seq, abri, comp]

            if seq == "ok":
                processed_samples += 1
            else:
                # If not "ok", process the sample using the one_sample function
                try:
                    result = pool.apply_async(
                        one_sample, args=(fastq_to_process[i],))
                    results.append(result)  # Store the result object
                    remaining_samples.append(i)  # Add the sample to the list of remaining samples

                except (ValueError, Exception) as e:
                    print(e)

        # Wait for the results of all submitted tasks
        for result in results:
            updated_info = result.get()
            # Store the updated info in the dictionary
            sample_name = fastq_to_process[i][0].split(os.sep)[-1].split("_")[0]
            updated_info_dict[sample_name] = updated_info

        # Update the samples_to_reprocess list with the remaining samples
        samples_to_reprocess = remaining_samples

        # If no samples were processed in this iteration, break out of the loop
        if processed_samples == 0:
            break

        time.sleep(10)

    # After processing, update the fastq_to_process list with the updated info
    '''
    for i, sample_info in enumerate(fastq_to_process):
        sample_name = sample_info[0].split(os.sep)[-1].split("_")[0]
        if sample_name in updated_info_dict:
            updated_info = updated_info_dict[sample_name]
            fastq_to_process[i] = sample_info + updated_info
    '''
except Exception as e:
    print(f"An exception occurred: {e}")

finally:
    # Close and join the pool to ensure proper cleanup
    pool.close()
    pool.join()
    
# Re-check the presence of "abri," "seq," and "comp" files for already processed samples
for i in range(len(fastq_to_process)):
    sample_name = fastq_to_process[i][0].split(os.sep)[-1].split("_")[0]
    sample_folder = os.path.join(output_path, sample_name)

    abri_file = os.path.join(sample_folder, sample_name + ".abricate")
    if os.path.isfile(abri_file):
        if os.path.getsize(abri_file):
            abri = "ok"
        else:
            abri = "no_contents"
    else:
        abri = "not_present"

    seq_file = os.path.join(
        sample_folder, sample_name + "_CompareTo_" + reference_name + "_good_snps.csv")
    if os.path.isfile(seq_file):
        if os.path.getsize(seq_file):
            seq = "ok"
        else:
            seq = "no_contents"
    else:
        seq = "not_present"

    comp_file = os.path.join(sample_folder, sample_name + "_CompareTo_" +
                             reference_name + "_good_snps_abricate_seqfinder.csv")
    if os.path.isfile(comp_file):
        if os.path.getsize(comp_file):
            comp = "ok"
        else:
            comp = "no_contents"
    else:
        comp = "not_present"

    fastq_to_process[i][-3:] = [seq, abri, comp]    
    

intro = [
    ["version", version, "", "", ""],
    ["version_date", date, "", "", ""],
    ["Start_time", start_time, "", "", ""],
    ["End_time", str(datetime.now().strftime("%d/%m/%Y %H:%M:%S")), "", "", ""],
    ["", "", "", "", ""],
    ["R1", "R2", "fasta", "file_size_check", "Seqfinder", "Abricate", "Compilation"]
]

write_csv(os.path.join(output_path, "seqfinder_summary.csv"), intro + fastq_to_process)

########## deleting abricate database
run_cmd(['rm','-r',abricate_ref_folder])

#### combine results
out_file_prefix=output_path.split(os.sep)[-2]+"_"+output_path.split(os.sep)[-1]
out_file=os.path.join(output_path,out_file_prefix+"__seqfinder_compilation.csv")
combine_tables_of_results(output_path,'*CompareTo*good_snps.csv',out_file)


out_file=os.path.join(output_path,out_file_prefix+"__abricate_seqfinder_compilation.csv")
combine_tables_of_results(output_path,'*_abricate_seqfinder.csv',out_file)

if database_type=="AMR":
    out_file=os.path.join(output_path,out_file_prefix+"__seqfinder_chr_compilation.csv")
    combine_tables_of_results(output_path,'*CompareTo*good_snps_only_chromosomal.csv',out_file)


summary_file_path= os.path.join(output_path,"seqfinder_summary.csv")
# Check if the file contains "not_present"
with open(summary_file_path) as file:
    if "not_present" in file.read():
        print("\n***Warning some of your files did not work***\n Please check the summary file for details")