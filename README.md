# APHA SeqFinder
This is a tool that determines the presence of genes in bacteria from illumina sequencing data.

Requirements: Forward and reverse reads that follow the illumina format of *R1_001.fastq.gz | *R2_001.fastq.gz

Edit file "template_arguments_file.args" and modify the global variables as appropiate. Following is an example:

#path to the seqFinder
soft_path="/home/javi/APHASeqFinder"
 
#database path and name. You can use any fasta file you wish. We this pipeline we provide several databases which you can find at /soft_path/references
reference="/home/javi/APHASeqFinder/references/AMR/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fna"

#path to the fastq files. Notice the name's format that the pipeline would be expecting (see above). 
data_path="/home/javi/WGS_Data/Project_1/fastq"

#path to the directory where the results will be placed
results_path="/home/javi/WGS_Results/Project_1"

#percentage of similarity between a sequence in the database and the sample being analysed to mark the gene as present.  
percentageID=70

#maximum number of snps in a sequence w.r.t. to the reference allowed to be called present.
numofsnps=5

#EFSA antimicrobial panel dictionary. Thi spipeline includes several versions at /soft_path/EFSA_panel
efsa_dict="/home/javi/APHASeqFinder/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"

You can save the new arguments with a different file name. For example: "/my/path/to/my_arguments_file_for_my_new_run.args"

To run the pipeline: python seqfinder.py /my/path/to/my_arguments_file_for_my_new_run.args

You will be asked how many cores you want to use. If you enter a wrong option, only 1 core will be used. For each sample, a folder with the same name as the sample will be created. The main results are summarised in a CSV file called "sample_name_CompareTo_database_name.csv".
Summary tables (CSV) for all the samples analysed can be found in folder results_path.
