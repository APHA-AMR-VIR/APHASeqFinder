# APHA SeqFinder
This is a tool that determines the presence of genes in bacteria from illumina sequencing data.
It runs first SeqFinder, then spades to create the assemblies (this step is optional, see below argument fastas_folder), then abricate and finally it combines both results.

The output files for each samples are in results_path/sample_name (see below how to set results_path)
Also in results_path, all the output files combining all the sample will be dropped.

Requirements: Forward and reverse reads. For example sample_R1_001.fastq.gz and sample_R2_001.fastq.gz
Dependencies: 

	- abricate (https://github.com/tseemann/abricate)
	- SPAdes (https://github.com/ablab/spades) 

How to run SeqFinder:

Edit file "template_arguments_file.args" and modify the global variables as appropiate. Following is an example, please read:

#path to the seqFinder
soft_path="/home/javi/APHASeqFinder"
 
#database path and name. You can use any fasta file you wish. We this pipeline we provide several databases which you can find at /soft_path/references
reference="/home/javi/APHASeqFinder/references/AMR/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fna"

#path to the fastq files. 
data_path="/home/javi/WGS_Data/Project_1/fastq"

#pattern followed by the R1 fastq file so seqFinder can find all the sample to run, e.i. some pattern unique in the R1 fastq files names that differentiate form other files. Examples:     
R1_pattern="_R1_"   if the R1 files is named sample_R1_001.fastq.gz

or

R1_pattern="_R1."  if the R1 files is named sample_R1.fastq.gz

#path to the directory where the results will be placed
results_path="/home/javi/WGS_Results/Project_1"

#EFSA antimicrobial panel dictionary. Thi spipeline includes several versions at /soft_path/EFSA_panel
efsa_dict="/home/javi/APHASeqFinder/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"

#fastas_folder. If the assemblies for the samples have already been done, specify the path where they are. If fasta_folder is left equal to "" seqfinder will run spades

fastas_folder=""           (if you want to run spades to create the assemblies)

or

fastas_folder="/home/javi/WGS_Results/Project_1/fastas"  (if the assemblies are already available in /home/javi/WGS_Results/Project_1/fastas)

#Number of cores to be used to run the samples in parallel. If you want to know how many cores your instance has got, type the command: htop 
ncores=2

#
#
#

You must save the new arguments with a different file name in any folder you want. For example: "/my/path/to/my_arguments_file_for_my_new_run.args"

To run the pipeline: ./seqfinder.py /my/path/to/my_arguments_file_for_my_new_run.args

For each sample, a folder with the same name as the sample will be created. The main results are summarised in a CSV file called "sample_name_CompareTo_database_name.csv".
Summary tables (CSV) for all the samples analysed can be found in folder results_path. To keep track of how you run seqfinder, a copy of the arguments file will be place in the results path.
