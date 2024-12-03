#!/usr/bin/env python

'''
APHASeqfinder
version 4.1.0
submitted to github on 04/10/2022
from pathlib import Path
from os import listdir
from datetime import datetime
Javier Nunez and Nick Duggett, AMR Team, Bacteriology
Animal and Plant Health Agency
'''

### packages
import sys,os,fnmatch,csv,os.path
from multiprocessing import Pool, Manager
from datetime import datetime
from pathlib import Path
from os import listdir,strerror
import re
import time
import subprocess
import shutil
import gzip
from errno import ENOTDIR

# ANSI escape codes for text Colours
class TextColours:
    HEADER = '\033[1;95m'
    OKBLUE = '\033[94m'
    WARNING = '\033[93m'
    QUESTION = '\033[36;1m'
    FAIL = '\033[1;91m'
    ENDC = '\033[0m'  # Reset Colour


print(TextColours.HEADER + "Checking python modules installed\n" + TextColours.ENDC)

class NotADirectoryError(OSError):
    def __init__(self, filename):
        self.filename = filename
        self.errno = ENOTDIR
        self.strerror = strerror(ENOTDIR)

def check_and_install(package_name, not_installed):
    try:
        __import__(package_name)
        print(TextColours.OKBLUE +f"Python module '{package_name}' is installed."+TextColours.ENDC)
        return True
    except ModuleNotFoundError:
        print(TextColours.FAIL + f"Required Python module '{package_name}' is not installed." + TextColours.ENDC)
        not_installed.append(package_name)
        return False

package_mapping = {'biopython': 'Bio'}

# List of required packages
required_packages = ['numpy', 'pandas', 'tqdm', 'Bio']

# Check each package
not_installed = []
for package in required_packages:
    mapped_package = package_mapping.get(package, package)
    if not check_and_install(mapped_package, not_installed):
        # Append the original package name to not_installed list for consistent display
        not_installed.append(package)

# If any required package is not installed, exit the script
if not_installed:
    print(TextColours.WARNING + "\nSome required packages are not installed. Exiting." + TextColours.ENDC)
    print(TextColours.WARNING + "To install the required packages: conda install -c anaconda", "numpy pandas tqdm biopython"+ TextColours.ENDC)
    exit()

import pandas as pd
from Bio import SeqIO
import numpy as np 
from tqdm import tqdm

def check_pandas_version():
    try:
        # Get the pandas version
        pandas_version = pd.__version__
        # Split the version string into major, minor, and micro versions
        major, minor, micro = map(int, pandas_version.split('.')[:3])  # Take only the first three elements
        return major, minor, micro
    except Exception as e:
        print("An error occurred while checking the pandas version:", e)
        return None, None, None

major, _, _ = check_pandas_version()
if major is not None and major < 2:
    print(TextColours.WARNING+f"Error: Your pandas version ({pd.__version__}) is not supported by this script. Please use pandas version 2.x."+TextColours.ENDC)
    print(TextColours.WARNING+"You may have to create a new python environment >=3.7 e.g.: conda create -n seqfinder python=3.10"+TextColours.ENDC)
    sys.exit(1)
else: 
    print(TextColours.OKBLUE +"\nPandas version check passed. Proceeding with the script.\n"+TextColours.ENDC)

def is_tool_installed(name):
    """Check whether 'name' is in PATH and marked as executable."""
    try:
        subprocess.run([name, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, NotADirectoryError):
        return False

def get_tool_version(name):
    """Get the version of the tool."""
    try:
        result = subprocess.run([name, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        version_output = result.stdout.decode('utf-8').strip()  # Decode bytes to string
        # Parsing the version from the output
        version = version_output.split()[1]
        return version
    except subprocess.CalledProcessError:
        return None


print(TextColours.HEADER + "Checking accessory programs are installed\n" + TextColours.ENDC)
tools = ['mlst', 'abricate']
required_versions = {'mlst': '2.23.0', 'abricate': '1.0.1'}

not_installed_tools = []
incorrect_version_tools = []

for tool in tools:
    if not is_tool_installed(tool): 
        not_installed_tools.append(tool)
    else:
        # Get the version of the tool
        tool_version = get_tool_version(tool)
        if tool_version:
            if tool_version != required_versions[tool]:
                incorrect_version_tools.append(tool)
        else:
            print(TextColours.FAIL+f"Unable to determine the version of {tool}."+TextColours.ENDC)

if not_installed_tools or incorrect_version_tools:
    print(TextColours.FAIL+f"The following tools are not installed or the incorrect version is installed: {', '.join(tools)}"+TextColours.ENDC)
    print(TextColours.WARNING+"Please install them with:\nSlower -->  conda install -c bioconda " + ' '.join([f"{tool}={required_versions[tool]}" for tool in tools])+TextColours.ENDC)
    print(TextColours.WARNING+"***OR***\nQuicker --> mamba install -c bioconda " + ' '.join([f"{tool}={required_versions[tool]}" for tool in tools])+TextColours.ENDC)
    exit()

# If all required packages are installed, proceed with the rest of your script
print(TextColours.OKBLUE +"All pre-requirements have been fulfilled. Continuing with the pipeline.\n"+TextColours.ENDC)

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

def run_cmd_silent(lis):
    os.system(" ".join(lis))

def run_subprocess_silent(cmd):
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")

### This section is all for checking that the average read length does not exceed the maximum permitted read length for the shortest reference sequence
### If the average read length exceeds the maximum, the reads are split in half and recombined. 
def find_shortest_reference_sequence(reference):
    # Initialise variables to store the shortest sequence information
    shortest_length = float('inf')
    shortest_name = None

    # Iterate over the sequences in the reference file
    for record in SeqIO.parse(reference, "fasta"):
        sequence_length = len(record.seq)

        # Update shortest sequence info ation if the current sequence is shorter
        if sequence_length < shortest_length:
            shortest_length = sequence_length
            shortest_name = record.name

    return shortest_name, shortest_length


def calculate_average_read_length(input_file, num_reads=1000):
    # Check if the file exists before opening
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        return 0

    total_length = 0
    read_count = 0

    with gzip.open(input_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Lines containing the sequence
                total_length += len(line.strip())
                read_count += 1

            if read_count == num_reads:
                break

    if read_count == 0:
        return 0  # Avoid division by zero

    average_length = total_length / read_count
    return average_length

def process_file(index, file):
    input_file = file[0]
    average_length = calculate_average_read_length(input_file)
    return {'Index': index, 'Filename': input_file, 'AverageLength': average_length}

def split_fastq(fastx_location, input_file, output_file_1, output_file_2):
    # Correctly generate the decompressed file path based on input_file
    decompressed_file = os.path.join(fastx_location, os.path.basename(input_file).replace('.gz', ''))
    print(f"Splitting {input_file} in half.\n")
    isolate_name = os.path.basename(input_file).split("_")[0]
    isolate_name_split = isolate_name + "_split"
    concatenated_output_file_name = os.path.basename(decompressed_file).replace(isolate_name,isolate_name_split)
    concatenated_output_file = os.path.join(fastx_location, concatenated_output_file_name)

    # Decompress the input gzip file to the specified output location
    with gzip.open(input_file, 'rt') as f_in, open(decompressed_file, 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)

    # Your existing implementation of calculate_average_read_length
    average_length = calculate_average_read_length(input_file)
    trim_length = average_length // 2

    # Use fastx_trimmer on the decompressed file
    run_cmd([os.path.join(soft_path,'third_party_software','fastx_trimmer'), '-i', decompressed_file, '-o', output_file_1, '-f', str(trim_length + 1)])
    run_cmd([os.path.join(soft_path,'third_party_software','fastx_trimmer'), '-i', decompressed_file, '-o', output_file_2, '-l', str(trim_length)])
    
    os.remove(decompressed_file)
    
    run_cmd(['cat',output_file_1, output_file_2,'>',concatenated_output_file])
    os.remove(output_file_1)
    os.remove(output_file_2)
    
    run_cmd(['gzip', concatenated_output_file])

    return trim_length, concatenated_output_file_name

def process_positive_index(index, shared_fastq_to_process, concatenated_output_file):
    positive_indexes_data = []
    sublist = shared_fastq_to_process[index]

    if sublist and len(sublist) >= 2:
        first_element = sublist[0]
        second_element = sublist[1]

        # Process the positive index within this loop
        input_file_forward = first_element
        input_file_reverse = second_element
        forward_filename = os.path.basename(first_element)
        reverse_filename = os.path.basename(second_element)
        output_file_forward_1 = os.path.join(fastx_location, f"{forward_filename}_output_file_1.fastq")
        output_file_forward_2 = os.path.join(fastx_location, f"{forward_filename}_output_file_2.fastq")
        output_file_reverse_1 = os.path.join(fastx_location, f"{reverse_filename}_output_file_1.fastq")
        output_file_reverse_2 = os.path.join(fastx_location, f"{reverse_filename}_output_file_2.fastq")
        
        trim_length_forward, concatenated_output_file_forward_name = split_fastq(fastx_location, input_file_forward, output_file_forward_1, output_file_forward_2)
        trim_length_reverse, concatenated_output_file_reverse_name = split_fastq(fastx_location, input_file_reverse, output_file_reverse_1, output_file_reverse_2)
        
        # Update the fastq_to_process list with the concatenated_output_file names
        isolate_name = os.path.basename(forward_filename).split("_")[0]
        isolate_name_split = isolate_name + "_split"
        concatenated_output_file_forward = input_file_forward.replace(isolate_name, isolate_name_split)
        concatenated_output_file_reverse = input_file_reverse.replace(isolate_name, isolate_name_split).replace('_R1', '_R2')
        # Concatenate fastx_location with the generated file names
        concatenated_output_file_forward = os.path.join(fastx_location, os.path.basename(concatenated_output_file_forward))
        concatenated_output_file_reverse = os.path.join(fastx_location, os.path.basename(concatenated_output_file_reverse))
    
        # Update the sublist
        sublist[0] = concatenated_output_file_forward
        sublist[1] = concatenated_output_file_reverse
        
        # Re-assign the list after updating the sublist
        shared_fastq_to_process[index]=sublist

        positive_indexes_data.append(sublist)

    return positive_indexes_data

def process_files_read_length(files_to_process, reference, fastx_location):
    manager = Manager()
    shared_fastq_to_process = manager.dict(enumerate(files_to_process))
    
    results = []

    with Pool() as pool:
        results = pool.starmap(process_file, [(index, file) for index, file in shared_fastq_to_process.items()])

    df = pd.DataFrame(results)

    # Find the shortest reference sequence
    shortest_reference_name, shortest_reference_length = find_shortest_reference_sequence(reference)
    max_length= int(shortest_reference_length / 0.75)
    print(TextColours.HEADER + f"\nChecking if average read-length for any isolate is 75% larger than the shortest reference sequence: {shortest_reference_name} -> {shortest_reference_length}bp\n" +TextColours.ENDC)
    # Check if each average length is 75% larger than the shortest reference sequence
    larger_than_75_percent = df['AverageLength'] > shortest_reference_length / 0.75
    
    # Print the number of positive indexes before processing
    num_positive_indexes = larger_than_75_percent.sum()
    
    if num_positive_indexes > 0:
        print(TextColours.WARNING+ f"{num_positive_indexes} isolates have been found with average read-lengths over the maximum of {max_length} and will be split in half.\n" +TextColours.ENDC)
        print(TextColours.WARNING+ "The files for these isolates will go through decompressing --> splitting in half with fastx --> concatenating the output --> compressing\n" +TextColours.ENDC)
        os.system('mkdir -p '+fastx_location)
    else:
        print(TextColours.OKBLUE+ "All isolates passed, continuing to next step.\n" +TextColours.ENDC)
    # Use Pool.starmap to parallelise the processing of positive indexes
    with Pool() as pool:
        concatenated_output_files = []
        for index in df[larger_than_75_percent].index:
            sublist = shared_fastq_to_process[index]
            first_element = sublist[0]
            input_file_forward = first_element
            concatenated_output_file = os.path.join(fastx_location, os.path.basename(input_file_forward).replace('.gz', ''))
            concatenated_output_files.append(concatenated_output_file)

        positive_indexes_data_list = pool.starmap(process_positive_index, [(index, shared_fastq_to_process, concatenated_output_file) for index, concatenated_output_file in zip(df[larger_than_75_percent].index, concatenated_output_files)])   # Combine results from all processes
    positive_indexes_data = [item for sublist in positive_indexes_data_list for item in sublist]

    # Add a short sleep to allow processes to complete
    time.sleep(1)

    # Convert shared_fastq_to_process dictionary back to a list of lists
    updated_fastq_to_process = [value for _, value in sorted(shared_fastq_to_process.items())]

    return df, updated_fastq_to_process



def trimming(sample_name,sample_folder,r1_file,r2_file):
    r1_file_trimmed_paired=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    r2_file_trimmed_paired=os.path.join(sample_folder,r2_file.split(os.sep)[-1].replace(".fastq.gz","_trimmed_paired.fastq.gz"))
    unpaired_trimmed=os.path.join(sample_folder,r1_file.split(os.sep)[-1].replace(".fastq.gz","_unpaired.fastq.gz"))
    run_cmd([os.path.join(soft_path,'third_party_software','fastp'),' -i '+r1_file+' -I '+r2_file+' -o '+r1_file_trimmed_paired+' -O '+r2_file_trimmed_paired+' --unpaired1 '+unpaired_trimmed+' --unpaired2 '+unpaired_trimmed+' -l 80 -r --cut_right_window_size 10 -w '+str(ncores_per_sample)+' -j /dev/null -h /dev/null'])
    fastq_trimmed=os.path.join(sample_folder,sample_name+'.trimmed.fastq.gz')
    run_cmd(['cat',r1_file_trimmed_paired,r2_file_trimmed_paired,unpaired_trimmed,'>',fastq_trimmed])
    run_cmd(['rm',r1_file_trimmed_paired])
    run_cmd(['rm',r2_file_trimmed_paired])
    run_cmd(['rm',unpaired_trimmed])
    if os.path.exists(fastx_location):
        run_cmd(['rm -r', fastx_location])
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
    print('\nNow combining files in ' + mother_path + ' that contains ' + pattern)
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
    

def get_mlst(fastq_to_process):
    try:
        fasta_file = fastq_to_process[2]
        mlst_output = os.path.join(output_path, "mlst.csv")
        run_cmd_silent(["mlst","--quiet", "--csv", fasta_file, ">>", mlst_output])
        return True
    except Exception as e:
        print("Something went wrong with the mlst stage:", str(e))
        return False
        
def failures(removed_isolates):
    try:
        r1_file=removed_isolates[0]
        sample_name=r1_file.split(os.sep)[-1].split("_")[0]
        sample_folder=os.path.join(output_path,sample_name)
        if not os.path.exists(sample_folder):
            os.system('mkdir -p '+sample_folder)
        columns = ['warnings','strain', 'id', 'gene', 'antimicrobial', 'class', 'ref_len', 'mapped_len', 'mean_depth', 'norm_depth', 'non_calls', 'perc_mapped', 'snps', 'other', 'good_snps', 'syn', 'non', 'annotation', 'reference', 'EFSA_dict', 'result_abr', 'contig_abr', 'plasmid_abr', 'coverage_abr', 'identity_abr', 'location_abr']
        failed_isolate_df = pd.DataFrame(columns=dict.fromkeys(columns, []))
        failed_isolate_df.loc[0,'strain']=sample_name
        failed_isolate_df = failed_isolate_df.fillna('No MLST scheme could be applied')
        output_filename = sample_name+"_CompareTo_"+reference_name+".csv"
        output_filename_good_snps=output_filename.replace(".csv","_good_snps.csv")
        output_filename_abricate_seqfinder=output_filename.replace('.csv','_good_snps_abricate_seqfinder.csv')
        output_file = os.path.join(sample_folder,output_filename)
        output_file_good_snps = os.path.join(sample_folder,output_filename_good_snps)
        output_file_abricate_seqfinder = os.path.join(sample_folder,output_filename_abricate_seqfinder)
        failed_isolate_df.to_csv(output_file, index=False)
        failed_isolate_df.to_csv(output_file_good_snps, index=False)
        failed_isolate_df.to_csv(output_file_abricate_seqfinder, index=False)
    except Exception as e:
        print("Something went wrong when dealing with failures:", str(e))

           
def build_reference(fastq_to_process):
    try:
        r1_file=fastq_to_process[0]
        sample_name=r1_file.split(os.sep)[-1].split("_")[0]
        sample_folder=os.path.join(output_path,sample_name)
        reference_folder=os.path.join(sample_folder,"reference_folder")
        abricate_folder=os.path.join(reference_folder,reference_name)
        if not os.path.exists(sample_folder):
            os.system('mkdir -p '+sample_folder)
        if not os.path.exists(reference_folder):
            os.system('mkdir -p '+reference_folder)
        if not os.path.exists(abricate_folder):
            os.system('mkdir -p '+abricate_folder)
        #ref_copy=os.path.join(abricate_ref_folder,reference_name+".fasta")
        run_cmd_silent(['cp',reference,os.path.join(reference_folder)])
        run_cmd_silent(['cp',reference,os.path.join(abricate_folder,"sequences")])
        run_subprocess_silent(['makeblastdb','-in',os.path.join(abricate_folder,"sequences"),'-dbtype','nucl','-hash_index'])
        return True
    except Exception as e:
        print("Something went wrong with the mlst stage:", str(e))
        return False

def fetch_error_info(comp_file):
    # Assuming the information to be added to "Errors" column is in the first row of the first column
    if os.path.exists(comp_file):
        with open(comp_file, 'r') as file:
            reader = csv.reader(file)
            next(reader, None)
            for row in reader:
                if row:
                    return row[0]
        
def one_sample(file_to_process):
    try:
        r1_file=file_to_process[0]
        r2_file=file_to_process[1]
        fasta_file=file_to_process[2]
        mlst_fasta=file_to_process[4]
        sample_name=r1_file.split(os.sep)[-1].split("_")[0]
        sample_folder=os.path.join(output_path,sample_name)
        reference_folder=os.path.join(sample_folder,"reference_folder")
        ref_index=reference.split(os.sep)[-1].split("/")[-1]
        ref_copy=os.path.join(reference_folder,ref_index)
        mlst_file_path = file_to_process[4]
        ref_copy_file_path = ref_copy
        # Read contents of mlst_fasta
        with open(mlst_file_path, 'r') as mlst_file:
            mlst_contents = mlst_file.read()
        # Read contents of ref_copy
        with open(ref_copy_file_path, 'r') as ref_copy_file:
            ref_copy_contents = ref_copy_file.read()       
        # Concatenate contents
        concatenated_contents = mlst_contents +'\n'+ ref_copy_contents     
        # Output path for the concatenated file
        concatenated_file_path = os.path.join(ref_copy)     
        # Write the concatenated contents to a new file
        with open(concatenated_file_path, 'w') as concatenated_file:
            concatenated_file.write(concatenated_contents)        
        current_pid = os.getpid()  # Get the current process's PID
        with open(os.path.join(output_path, "processed.csv"), "a") as f:
            f.write(sample_name + " processed by PID "+str(current_pid)+"\n")
        with open(os.path.join(output_path, "processed.csv"), "r") as f:
           lines=len(f.readlines())
        print(TextColours.HEADER + f"Processing sample: {sample_name} ({lines}/{len(fastq_to_process)})" + TextColours.ENDC)
        sample_folder=os.path.join(output_path,sample_name)
        if not os.path.exists(sample_folder):
            os.system('mkdir -p '+sample_folder)
        fastq_trimmed=trimming(sample_name,sample_folder,r1_file,r2_file)
        #ref_index=os.path.join(sample_folder,ref_copy.split(os.sep)[-1])
        run_cmd([os.path.join(soft_path,'third_party_software','smalt'),'index -k 13 -s 6',ref_copy,concatenated_file_path])
        run_cmd([os.path.join(soft_path,'third_party_software','smalt'),f'map -n {ncores_per_sample} -d -1 -y 0.7 -x -f samsoft -o',os.path.join(sample_folder,sample_name+".sam"),ref_copy, fastq_trimmed])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Sh -F 4',os.path.join(sample_folder,sample_name+".sam"),'>',os.path.join(sample_folder,sample_name+".F4.sam")])
        sam_flag_changed(os.path.join(sample_folder,sample_name+".F4.sam")) #produces .60.sam
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'view -Shu',os.path.join(sample_folder,sample_name+'.60.sam'),'>',os.path.join(sample_folder,sample_name+'.bam')]) 
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'sort',os.path.join(sample_folder,sample_name+'.bam'),os.path.join(sample_folder,sample_name+'.sorted')])       
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'index',os.path.join(sample_folder,sample_name+'.sorted.bam')])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'faidx',concatenated_file_path])
        run_cmd([os.path.join(soft_path,'third_party_software','samtools_0.1.19'),'mpileup','-q -1 -uf',concatenated_file_path,os.path.join(sample_folder,sample_name+'.sorted.bam'),'>',os.path.join(sample_folder,sample_name+'.mpileup.bcf')])    
        run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -cg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.mpileup.vcf')])
        run_cmd(['perl',os.path.join(soft_path,'third_party_software','vcfutils.pl'),'vcf2fq',os.path.join(sample_folder,sample_name+'.mpileup.vcf'),'>',os.path.join(sample_folder,sample_name+'.fq')])   
        run_cmd([os.path.join(soft_path,'third_party_software','bcftools_0.1.19'),'view -vcg',os.path.join(sample_folder,sample_name+'.mpileup.bcf'),'>',os.path.join(sample_folder,sample_name+'.snp')])
        snps_filter(th_qual,th_prop,th_min,os.path.join(sample_folder,sample_name+'.snp'),os.path.join(sample_folder,sample_name+'_SN.csv'))
        #pileup_stats_file=pileupStats(os.path.join(sample_folder,sample_name+".mpileup.vcf"),os.path.join(sample_folder,sample_name+".fq"),ref_copy,3)
        run_cmd(['python',os.path.join(soft_path,'calculate_genes_presence.py'),sample_folder,concatenated_file_path,mlst_fasta,str(th_qual),str(th_prop),str(th_min)])
        
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
        print(TextColours.FAIL +"Something went wrong with the mapping stage."+ TextColours.ENDC)
        print(TextColours.FAIL + "Therefore, "+sample_name +
              " wasn't processed." + TextColours.ENDC)
        #os.system("find "+sample_folder +
        #          " ! -name mlst_error.txt -type f -exec rm {} +")
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
            run_cmd(["rm -r "+reference_folder])
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
            run_cmd(["rm -r "+reference_folder])
            with open(failed_seqfinder_file_name, "w") as f:
                f.write(contents)
            raise Exception(TextColours.FAIL +'Exiting function, no genes passed filter in '+sample_name+ TextColours.ENDC)
            #exit()
        else:
            os.system("abricate --datadir {} --mincov 10 --minid 10 --db {} {} > {}".format(reference_folder, reference_name, fasta_file, abricate_file_name))
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
    run_cmd(["rm -r "+reference_folder])
    
###################################
###################################
### Starting
start_time=str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
###################################
###################################

########### Global variables
version='4.1.0'
date='23/05/2024'

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
mlst_scheme_dictionary = ""
mlst_picker=""
results_path=""
efsa_dict=""
vir_dict=""
bio_metal_dict=""
fastas_folder=""
ncores_per_sample=1
ncores=1
sample_list_file=""
#mlst_fasta=""
###########Additional variables for folder naming
today= datetime.now()
date_of_run=today.strftime('%Y%m%d')
reference_name=reference.split(os.sep)[-1].split(".")[0]
output_path=os.path.join(results_path,reference_name,date_of_run)
output_path_logging = os.path.join(output_path, "logging.txt")

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
#if not os.path.isfile(mlst_fasta):
#    print('\nMLST file cannot be located, make sure you have provided the correct path above.')
#    sys.exit()
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
print(TextColours.HEADER +"\nChecking your reference file\n"+ TextColours.ENDC)
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
    print(TextColours.OKBLUE +"No errors found in the reference file\n"+ TextColours.ENDC)
'''
# Call the function to find the shortest sequence
shortest_reference_name, shortest_reference_length = find_shortest_reference_sequence(reference)

# Print the Name and length of the shortest reference sequence
if shortest_reference_name:
    print(f"Shortest Reference Sequence: {shortest_reference_name}")
    print(f"Length: {int(shortest_reference_length)}")
else:
    print("No sequences found in the reference file.")
'''



########## checking that R2 and assembly files exist

fastq_to_process=[]
sample_names = set()  # Store sample names to check for duplicates
R2_pattern=R1_pattern.replace("R1","R2")
today= datetime.now()
date_of_run=today.strftime('%Y%m%d')
reference_name=reference.split(os.sep)[-1].split(".")[0]
output_path=os.path.join(results_path,reference_name,date_of_run)
print(TextColours.HEADER +"\nChecking samples to be run\n"+ TextColours.ENDC)
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
                        print(TextColours.WARNING+os.path.join(data_path, fastq_R1) + " is less than 10MB. This might not provide an output."+TextColours.ENDC)
                    if os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:
                        print(TextColours.WARNING+os.path.join(data_path, fastq_R2) + " is less than 10MB. This might not provide an output."+TextColours.ENDC)
                else:
                    file_size="good"
                if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000 or os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:   
                    remove_small_files=""
                    while remove_small_files.lower() not in ['y', 'n']:
                        remove_small_files=input(TextColours.QUESTION+"Do you want to remove this isolate from your SeqFinder run (y/n)?\n"+TextColours.ENDC)
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
                    print(TextColours.FAIL+"Missing assembly file for file: "+sample_name+TextColours.ENDC)
                else:
                    fasta_ok="several_found"
                    print(TextColours.FAIL+"More than one assembly for: "+sample_name+TextColours.ENDC)
            else:
               R2_ok="not_found"
               print(TextColours.FAIL+"Missing R2 for file: "+fastq_R1+TextColours.ENDC)
            
            if R2_ok=="found" and fasta_ok=="found":
                fastq_to_process.append([os.path.join(data_path,fastq_R1),os.path.join(data_path,fastq_R2),os.path.join(fastas_folder,fasta_file),file_size])
                #This sorts the input into the script by the third column so that if small files are chosen to be included by the user they are processed last
                fastq_to_process.sort(key=lambda x: x[3])
                
            summary.append([fastq_R1,R2_ok,fasta_ok,file_size])
else:
    if not os.path.isfile(sample_list_file):
        print(TextColours.FAIL+"Cannot find sample_list_file at "+sample_list_file+TextColours.ENDC)
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
                print(TextColours.FAIL+"Missing file: "+fastq_R1+TextColours.ENDC)
            if os.path.isfile(fastq_R2):
                R2_ok="found"
                if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000 or os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:
                    file_size="small"
                    if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000:
                        print(TextColours.WARNING+os.path.join(data_path, fastq_R1) + " is less than 10MB. This might not provide an output."+TextColours.ENDC)
                    if os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:
                        print(TextColours.WARNING+os.path.join(data_path, fastq_R2) + " is less than 10MB. This might not provide an output."+TextColours.ENDC)
                else:
                    file_size="good"
                if os.path.getsize(os.path.join(data_path, fastq_R1)) < 10000000 or os.path.getsize(os.path.join(data_path, fastq_R2)) < 10000000:   
                    remove_small_files=""
                    while remove_small_files.lower() not in ['y', 'n']:
                        remove_small_files=input(TextColours.QUESTION+"Do you want to remove this isolate from your SeqFinder run (y/n)?\n"+TextColours.ENDC)
                    if remove_small_files.lower() == "y":
                        R2_ok="too_small_excluded_by_user"
                        print("\nRemoved\n")
                    elif remove_small_files.lower() == "n":
                        print("\nKept\n")   
            else:
                R2_ok="not_found"
                print(TextColours.FAIL+"Missing file: "+fastq_R2+TextColours.ENDC)
                
            if os.path.isfile(fasta_file):
                fasta_ok="found"
            else:
                fasta_ok="not_found"
                print(TextColours.FAIL+"Missing file: "+fasta_file+TextColours.ENDC)
            if R1_ok=="found" and R2_ok=="found" and fasta_ok=="found":
                fastq_to_process.append([fastq_R1,fastq_R2,fasta_file,file_size])
            
            summary.append([fastq_R1,R2_ok,fasta_ok,file_size])
print(TextColours.OKBLUE+"Check is complete\n"+ TextColours.ENDC)


run_cmd(['mkdir -p',output_path])
write_csv(os.path.join(output_path,"initial_file_testing.csv"),summary)
print('********\nOutput will be stored at:'+output_path+'\n********')    
print("Check file "+os.path.join(output_path,"summary.csv")) 
print("\n***** There are "+str(len(fastq_to_process))+" isolates that APHASeqFinder will process.")
if len(fastq_to_process) == 0:
   print(TextColours.FAIL+"\nNo fastq files were detected. Ensure file path is correct and fastq files end in .fastq.gz"+TextColours.ENDC)
   sys.exit()   

fastx_location=os.path.join(output_path,"fastx")

value=input(TextColours.QUESTION+'\nWould you like to proceed (y/n)?\n'+TextColours.ENDC)
if value!='y':
    sys.exit()

result_dataframe, updated_fastq_to_process = process_files_read_length(fastq_to_process, reference, fastx_location)

#Updating fastq_to_process if required
fastq_to_process = updated_fastq_to_process

print(TextColours.HEADER + "\nRunning mlst for each isolate to determine best scheme for normalisation" + TextColours.ENDC)

pool = Pool(ncores)

try:
    # Use tqdm to create a progress bar
    results = list(tqdm(pool.imap(get_mlst, fastq_to_process), total=len(fastq_to_process)))
except (ValueError, Exception) as e:
    print(e)
finally:
    pool.close()
    pool.join()

print(TextColours.OKBLUE + "Finished mlst" + TextColours.ENDC)

    

############### parallel processing of mlst

mlst_path=os.path.join(output_path,"mlst.csv")


# Function to replace commas within parentheses with semicolons when there are multiple calls for an allele
def replace_commas_in_parentheses(line):
    return re.sub(r'\(([^)]+)\)', lambda m: '(' + m.group(1).replace(',', ';') + ')', line)

# Update the mlst file with the commas replaced by semicolons where required
corrected_lines = []
max_columns = 0

with open(mlst_path, 'r') as file:
    for line in file:
        corrected_line = replace_commas_in_parentheses(line.strip())
        corrected_lines.append(corrected_line)
        #Determine number of columns in file after replacing commas with semi-colons
        num_columns = len(corrected_line.split(','))
        if num_columns > max_columns:
            max_columns = num_columns

# Write the corrected lines back to the same file
with open(mlst_path, 'w') as file:
    for line in corrected_lines:
        file.write(line + '\n')

# Generate the column names with the nuber of allele columns elastic based on highest number of columns in any line
mlst_results_column_names = ['Fasta_file', 'scheme', 'ST'] + [f'allele_{i+1}' for i in range(max_columns - 3)]

# Read the corrected data into a pandas DataFrame
mlst_results = pd.read_csv(mlst_path, header=None, names=mlst_results_column_names)
mlst_results = mlst_results.apply(lambda x: x.replace('-', 'unknown'))
mlst_results = mlst_results.apply(lambda x: x.fillna('unknown'), axis=1)

#Making MLST scheme dictionary work for all users
mlst_scheme_dictionary_df = pd.read_csv(mlst_scheme_dictionary)
#Take everything before last backslash in arguments file input
# Extract the user-specific part from the mlst_scheme_dictionary variable
user_specific_part_mlst_scheme_dictionary = mlst_scheme_dictionary.rsplit('/', 1)[0]
# Replace everything before the last backslash in the filepaths
mlst_scheme_dictionary_df['filepath'] = mlst_scheme_dictionary_df['filepath'].apply(lambda x: user_specific_part_mlst_scheme_dictionary + '/' + x.rsplit('/', 1)[-1])

# Now your modified scheme_dict is ready for use
scheme_dict = mlst_scheme_dictionary_df.set_index('scheme')['filepath'].to_dict()

#scheme_dict = pd.read_csv(mlst_scheme_dictionary, header=None, index_col=0).to_dict()[1]
mlst_results['path'] = mlst_results.iloc[:, 1].map(scheme_dict)
mlst_results = mlst_results.iloc[:, [0,1, -1]]

to_remove = []  # List to store indices of sublists to be removed
removed_isolates = []  # List to store removed isolates

for sublist in fastq_to_process:
    match_found = False
    
    if mlst_picker!="":
        print("reading mlst picker")
        mlst_picker=pd.read_csv(mlst_picker)
        for _, row in mlst_picker.iterrows():
            if sublist[2] == row['Fasta_full_path']:
                sublist.append(row["Scheme"])
                match_found = True
                break

    if not match_found:
        matching_row = mlst_results[mlst_results['Fasta_file'] == sublist[2]]
        if not matching_row.empty:
            corresponding_path = matching_row['path'].iloc[0]
            if pd.notna(corresponding_path):  # Check for NaN values
                sublist.append(corresponding_path)
                match_found = True

    if not match_found:
        # If still no match is found, mark for removal
        to_remove.append(sublist)
        removed_isolates.append(sublist)
        print(TextColours.FAIL + f"Isolate {sublist[0]} failed to be designated an MLST scheme.\n" + TextColours.ENDC)

# Remove sublists marked for removal
fastq_to_process = [sublist for sublist in fastq_to_process if sublist not in to_remove]

#If there are any isolates that failed to be designated an MLST scheme then process them here so they are included in downstream functions
if removed_isolates:
    print("Dealing with failed isolates.\n")
    pool = Pool(ncores)
    try:
        pool.map(failures,removed_isolates)
    except(Exception) as e:
        print(e)

print(TextColours.HEADER+"\nBuilding reference database for each isolate"+ TextColours.ENDC)

pool = Pool(ncores)
try:
    results = list(tqdm(pool.imap(build_reference, fastq_to_process), total=len(fastq_to_process)))
except (ValueError, Exception) as e:
     print(e)
finally:
    pool.close()
    pool.join()
print(TextColours.OKBLUE+"Finished building reference databases"+ TextColours.ENDC)


#Remove processing file if it already exists   
processed_csv_path = os.path.join(output_path, "processed.csv")
if os.path.isfile(processed_csv_path):
    run_cmd_silent(["rm "+processed_csv_path])


# Initialise a dictionary to store the updated information
updated_info_dict = {}

samples_to_reprocess = list(range(len(fastq_to_process)))
pool = Pool(ncores)

try:
    while samples_to_reprocess:
        print(TextColours.HEADER+"\n***Starting sample processing***\n"+ TextColours.ENDC)
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
    if "not_present" in [seq, abri, comp]:
        error_info = fetch_error_info(comp_file)
        fastq_to_process[i].append(error_info)

intro = [
    ["version", version, "", "", ""],
    ["version_date", date, "", "", ""],
    ["Start_time", start_time, "", "", ""],
    ["End_time", str(datetime.now().strftime("%d/%m/%Y %H:%M:%S")), "", "", ""],
    ["", "", "", "", ""],
    ["R1", "R2", "fasta", "file_size_check", "mlst scheme","Seqfinder", "Abricate", "Compilation","Issues"]
]

write_csv(os.path.join(output_path, "seqfinder_summary.csv"), intro + fastq_to_process)


#### combine results
out_file_prefix=output_path.split(os.sep)[-2]+"_"+output_path.split(os.sep)[-1]
out_file_seqfinder_compilation=os.path.join(output_path,out_file_prefix+"__seqfinder_compilation.csv")
combine_tables_of_results(output_path,'*CompareTo*good_snps.csv',out_file_seqfinder_compilation)


out_file_abricate_seqfinder_compilation=os.path.join(output_path,out_file_prefix+"__abricate_seqfinder_compilation.csv")
combine_tables_of_results(output_path,'*_abricate_seqfinder.csv',out_file_abricate_seqfinder_compilation)

out_file=os.path.join(output_path,out_file_prefix+"__abricate_compilation.csv")
combine_tables_of_results(output_path,'*.abricate',out_file)

if database_type=="AMR":
    out_file=os.path.join(output_path,out_file_prefix+"__seqfinder_chr_compilation.csv")
    combine_tables_of_results(output_path,'*CompareTo*good_snps_only_chromosomal.csv',out_file)



summary_file_path= os.path.join(output_path,"seqfinder_summary.csv")
isolates_without_genes = "No genes passed filter"
isolate_names_without_genes = []

# Check if the file contains "not_present" or "No genes passed filter"
with open(summary_file_path) as file:
    file_content = file.read()

    if "not_present" in file_content:
        print(TextColours.FAIL + "\n***Warning - some of your files did not work***\nPlease check the summary file for details: "+summary_file_path+ TextColours.ENDC)

    # Find isolates without genes
    for line in file_content.split('\n'):
        if isolates_without_genes in line:
            isolate_name = line.split(",")[0].split("/")[-1].split("_")[0]
            isolate_names_without_genes.append(isolate_name)

try:
    abri_seq=pd.read_csv(out_file_abricate_seqfinder_compilation)
    genes_to_check = [
        "colis-g2197_mcr1", "colis-g2231_mcr2", "colis-g2238_mcr1-2", "colis-g2239_mcr1-3",
        "colis-g2240_mcr1-4", "colis-g2241_mcr1-5", "colis-g2242_mcr1-6", "colis-g2243_mcr1-7",
        "colis-g2244_mcr1-8", "colis-g2245_mcr1-9", "colis-g2246_mcr1-10", "colis-g2248_mcr2-2",
        "colis-g2249_mcr6-1", "colis-g2250_mcr3", "colis-g2254_mcr4", "colis-g2256_mcr4-3",
        "colis-g2257_mcr5-1", "colis-g2258_mcr7-1", "colis-g2317_mcr8-1", "colis-g2333_mcr9-1",
        "colis-g2334_mcr10-1", "colis-g2831_mcr1-11", "colis-g2832_mcr1-12", "colis-g2833_mcr1-13",
        "colis-g2834_mcr1-14", "colis-g2835_mcr1-26", "colis-g2836_mcr1-27", "colis-g2837_mcr1-9",
        "colis-g2838_mcr3-10", "colis-g2839_mcr3-11", "colis-g2840_mcr3-12", "colis-g2841_mcr3-13",
        "colis-g2842_mcr3-14", "colis-g2843_mcr3-15", "colis-g2844_mcr3-16", "colis-g2845_mcr3-17",
        "colis-g2846_mcr3-18", "colis-g2847_mcr3-19", "colis-g2848_mcr3-2", "colis-g2849_mcr3-20",
        "colis-g2850_mcr3-21", "colis-g2851_mcr3-22", "colis-g2852_mcr3-23", "colis-g2853_mcr3-24",
        "colis-g2854_mcr3-25", "colis-g2855_mcr3-3", "colis-g2856_mcr3-4", "colis-g2857_mcr3-5",
        "colis-g2858_mcr3-6", "colis-g2859_mcr3-7", "colis-g2860_mcr3-8", "colis-g2861_mcr3-9",
        "colis-g2863_mcr4-2", "colis-g2864_mcr4-4", "colis-g2865_mcr4-5", "colis-g2866_mcr4-6",
        "colis-g2867_mcr5-2","betaL-g1033_OXA-48",
        "Vir0018_stx2B_1", "Vir0093_stx2A_1", "Vir0378_stx2A_2", "Vir0410_stx2A_3", "Vir1017_stx1B_1",
        "Vir1258_stx1A_1", "Vir1319_stx2A_4", "Vir1512_stx2B_2", "Vir2206_stx2e-A_2", "Vir2376_stx2A_5",
        "Vir2449_stx1A_2", "Vir2775_stx1A_3", "Vir3432_stx1A_4", "Vir3691_stx1a-A_1", "Vir3693_stx1c-A_1",
        "Vir3694_stx1c-B_1", "Vir3695_stx1d-A_1", "Vir3696_stx1d-B_1", "Vir3697_stx1a-A_2", "Vir3699_stx2b-A_1",
        "Vir3700_stx2b-B_1", "Vir3701_stx2c-A_1", "Vir3702_stx2c-B_1", "Vir3703_stx2d-A_1", "Vir3704_stx2d-B_1",
        "Vir3705_stx2e-A_1", "Vir3706_stx2e-B_1", "Vir3707_stx2f-A_1", "Vir3708_stx2f-B_1", "Vir3709_stx2g-A_1",
        "Vir3710_stx2g-B_1", "Vir3726_stx2A_6", "Vir3740_stx2A_7", "Vir3751_stx2A_8", "Vir3754_stx2A_9",
        "Vir3758_stx2A_10", "Vir3759_stx1B_2", "Vir3762_stx2A_11", "Vir3766_stx2A_12", "Vir3767_stx2A_13",
        "Vir3779_stx2A_14", "Vir3783_stx2A_15", "Vir3785_stx1B_3", "Vir3799_stx1B_4", "Vir3801_stx2A_16",
        "Vir3813_stx2A_17", "Vir3815_stx2A_18", "Vir3821_stx1A_5", "Vir3822_stx1B_5", "Vir3825_stx2A_19",
        "Vir3833_stx2A_20", "Vir3834_stx2A_21", "Vir3844_stx2A_22", "Vir3845_stx1B_6", "Vir3860_stx1B_7",
        "Vir3863_stx2A_23", "Vir3866_stx2A_24", "Vir3870_stx1A_6", "Vir3872_stx2A_25", "Vir3883_stx2A_26",
        "Vir3889_stx2A_27", "Vir3891_stx2A_28", "Vir3893_stx1A_7", "Vir3894_stx1A_8", "Vir3896_stx1A_9",
        "Vir3900_stx2A_29", "Vir3903_stx2A_30", "Vir3904_stx1A_10", "Vir3908_stx2A_31", "Vir3910_stx2A_32",
        "Vir3914_stx2A_33", "Vir3918_stx2A_34", "Vir3920_stx2A_35", "Vir3923_stx1A_11", "Vir3924_stx2A_36",
        "Vir3925_stx2A_37", "Vir3926_stx2A_38", "Vir3930_stx2A_39", "Vir3932_stx2A_40", "Vir3934_stx2A_41",
        "Vir3936_stx2A_42", "Vir3942_stx2A_43", "Vir3944_stx1B_8", "Vir3948_stx1A_12", "Vir3953_stx2A_44",
        "Vir3956_stx2A_45", "Vir3957_stx2A_46", "Vir3966_stx2A_47", "Vir3967_stx2A_48", "Vir3968_stx2A_49",
        "Vir3971_stx2A_50", "Vir3972_stx2A_51", "Vir3974_stx2A_52", "Vir3983_stx2A_53", "Vir3984_stx2A_54",
        "Vir3986_stx2A_55", "Vir3987_stx2A_56", "Vir4002_stx2A_57", "Vir4003_stx2A_58", "Vir4011_stx1A_13",
        "Vir4015_stx2A_59", "Vir4027_stx2A_60", "Vir4030_stx2A_61", "Vir4035_stx1B_9", "Vir4041_stx1B_10",
        "Vir4042_stx1A_14", "Vir4044_stx2A_62", "Vir4045_stx2A_63", "Vir4047_stx2A_64", "Vir4049_stx1A_15",
        "Vir4050_stx2A_65", "Vir4051_stx2A_66", "Vir4052_stx2A_67", "Vir4057_stx2A_68", "Vir4058_stx2A_69",
        "Vir4060_stx1B_11", "Vir4061_stx2A_70", "Vir4062_stx2A_71", "Vir4065_stx2A_72", "Vir4077_stx1A_16",
        "Vir4082_stx1B_12", "Vir4114_stx2A_73", "Vir4116_stx1A_17", "Vir4119_stx2A_74", "Vir4121_stx2A_75",
        "Vir4126_stx2A_76", "Vir4129_stx2A_77", "Vir4134_stx2A_78", "Vir4138_stx2A_79", "Vir4145_stx1A_18"
    ]    
    for idx, row in abri_seq.iterrows():
        found_genes = [gene for gene in genes_to_check if gene in row['id']]
        if found_genes:
            print(TextColours.FAIL+"\n Warning "+row['strain'] + " contains " + ', '.join(found_genes)+TextColours.ENDC)
except Exception as e:
    print("Error:", e)
    print("An error occurred whilst looking for specific genes")
        
    


# Report isolates without genes
if isolate_names_without_genes:
    print(TextColours.FAIL + "\nIsolates with no genes passing filter:" + TextColours.ENDC)
    for isolate_name in isolate_names_without_genes:
        print(TextColours.WARNING + isolate_name + "\n" + TextColours.ENDC)

if removed_isolates:
    print(TextColours.FAIL+"\nIsolates not processed because appropriate MLST scheme could not be established:"+TextColours.ENDC)
    for isolate in removed_isolates:
        modified_info = isolate[0].split("/")[-1].split("_")[0]
        print(TextColours.WARNING+modified_info+"\n"+TextColours.ENDC)
        

print(TextColours.HEADER+"\nFinished"+TextColours.ENDC)
