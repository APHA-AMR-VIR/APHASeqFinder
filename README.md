# APHA SeqFinder
This is a tool that determines the presence of AMR genes in bacteria from illumina sequencing data.
It runs first SeqFinder and then abricate and finally it combines both results.

The output files for each samples are in results_path/sample_name (see below how to set results_path)
Also in results_path, all the output files combining all the sample will be dropped.

Requirements: Forward and reverse reads. For example sample_R1_001.fastq.gz and sample_R2_001.fastq.gz and assembliy fasta files. The fasta files' names have to match with the fastq R1 file up to the first "_" which is taken as the sample_name. Also the fasta files extension must be ".fasta" or ".fa"
Dependencies: 

	- abricate (https://github.com/tseemann/abricate)
	- SPAdes (https://github.com/ablab/spades) 

How to run SeqFinder:

Edit file "template_arguments_file.args" and modify the global variables as appropiate. Following is an example, please read:

```

### This is the template file that you need to alter. You must provide the direct paths you want to use. A quick way to do this would be to do a find a replace of "nickduggett" with your own name. This can be done in the command-line with sed (sed 's/nickduggett/yournamegoeshere/g' template_arguments_file.args > template_arguments_file_yourname.args)

#path to the seqFinder
soft_path="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2"
 
#database path and name. You can use any fasta file you wish. We this pipeline we provide several databases which you can find at /soft_path/references. Delete "#" to select that the reference database you want

reference="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/references/AMR/AMRDatabase_20200729_and_EnteroPLasmids_20190514_short-tetA6.fna"
#reference="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/references/virulence/VirulenceFactors_2022_05_26.fna"
#reference="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/references/disinfectant/disinfectant_2022_06_23.fna"

#database type; only four possible values "AMR", "VIR", "DIS" or "other". Delete "#" to select that database type
database_type="AMR"
#database_type="VIR"
#database_type="DIS"
#database_type="other"

## MLST fasta for depth normalization
mlst_fasta="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/references/mlst/ECO-MLST-MG1655-alleles.fna"

#EFSA antimicrobial panel dictionary. This pipeline includes several versions at /soft_path/EFSA_panel
efsa_dict="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/EFSA_panel/EFSA_antimcriobial_panel_dictionary_191219.csv"

#Virulence dictionary containing information about gene function and associated pathotypes where available
vir_dict="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/references/virulence/vir_dict_2022_06_17.csv"

#Disinfectant dictionary containing information about gene class, phenotype and function where available
dis_dict="/home/nickduggett/mnt/fsx-044/Bactipipes/BactiPipes_SCE3/AMR_github/APHASeqFinder_4.0.2/references/disinfectant/disinfectant_dictionary_2022_06_23.csv"

#path to the fastq files. Notice the name's format that the pipeline would be expecting (see above). 
data_path='/home/nickduggett/mnt/fsx-044/RDVM0529/VM0529/2016_2018_2020_paper/poppunk/2020_turkey/enterobase_download/matches'

#fastq R1 pattern
R1_pattern='_R1_'

#fastas_folder
fastas_folder='/home/nickduggett/mnt/fsx-044/RDVM0529/VM0529/2016_2018_2020_paper/poppunk/2020_turkey/enterobase_download/matches'

## comma separated table containing a list of R1_fastq and corresponding fasta (same line, one sasmple each line)
## This argument overwrite data_path and fastas_folder.
#sample_list_file='/home/nickduggett/wgs_data/test/samples_list.csv'
sample_list_file=''

#path to the directory where the results will be placed
results_path='/home/nickduggett/mnt/fsx-044/RDVM0529/VM0529/2016_2018_2020_paper/poppunk/2020_turkey/enterobase_download/matches'

#ncores
ncores=4

```

#
#
#

You must save the new arguments with a different file name in any folder you want. For example: "/my/path/to/my_arguments_file_for_my_new_run.args"

To run the pipeline: ./seqfinder.py /my/path/to/my_arguments_file_for_my_new_run.args

For each sample, a folder with the same name as the sample will be created. The main results are summarised in a CSV file called "sample_name_CompareTo_database_name.csv".
Summary tables (CSV) for all the samples analysed can be found in folder results_path. To keep track of how you run seqfinder, a copy of the arguments file will be place in the results path. Aso a seqfinder_summary.csv table is generated to let you know if the samples (and what samples) have been succesfully run.
