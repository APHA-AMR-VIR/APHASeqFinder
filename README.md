# APHA SeqFinder
This is a tool that determines the presence of genes in bacteria from illumina sequencing data.

Requirements:
Forward and reverse reads that follow the illumina format of *R1_001.fastq.gz | *R2_001.fastq.gz

Edit file "template_arguments_file.args" and modify the global variables as appropiate. Following is an example:

soft_path="/home/javi/APHASeqFinder"
reference="/home/javi/APHASeqFinder/references/AMR/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fna"
data_path="/home/javi/WGS_Data/Project_1/fastq"
results_path="/home/javi/WGS_Results/Project_1"
percentageID=70
numofsnps=5
efsa_dict="/home/javi/APHASeqFinder/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"

run: "python seqfinder.py" if you saved your arguments file as  "template_arguments_file.args"
or run: "python seqfinder.py my_arguments_file.args" if you changed the arguments file's name.




