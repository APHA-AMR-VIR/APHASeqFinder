# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
#

# This bash script will run the AMR gene analysis for any number of squencing runs.
#
# go to:
# 	cd ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1/software/WGS_GeneMapper
#Running from kraken machine to ~ included in path. if running from own machine remove the ~ in the path.

# and run:
#	sh WGS_GeneMapper_AMRdatabase_20170815.sh
# 	sh /mnt/Molsig_VM0518/Miranda_Ecoli_Aug17/WGS_GeneMapper_AMRdatabase_20170815.sh

# Every line calls WGS_GeneMapper_Run providing the following arguments:
# 	arg1 the reference fasta database containing the AMR genes sequences in fasta format
# 	arg2 the name of the run (folder name where the fastq files are stored).
# 	arg3 number of cpus to be used
#
#sleep 5h
#  /mnt/Molsig_VM0518/Programs/SeqFinder/software/WGS_GeneMapper
# the path/link will depend on you machine!  eg softpath="/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder/software/WGS_GeneMapper"
softpath="/home/p170122/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1/software/WGS_GeneMapper"
refpath="/home/p170122/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1/software/references"
#refpath="/mnt/Molsig_VM0518/Programs/Seqfinder_v1/software/references" 


#	the line below must contain a folder labelled fastqs within folder called fastqs
#datapath="/mnt/Molsig_VM0518/Raw_Seq_Reads_DO.NOT.EDIT/Seqfinder_example_folder" 
#datapath="/mnt/Molsig_VM0518/Raw_Seq_Reads_DO.NOT.EDIT/Ecoli5" 
datapath="/mnt/Molsig_VM0518/Miranda_Ecoli_Aug17" 


#resultspath="/mnt/Molsig_VM0518/Raw_Seq_Reads_DO.NOT.EDIT/Seqfinder_example_folder/Seqfinder_results"
#resultspath="/mnt/Molsig_VM0518/Raw_Seq_Reads_DO.NOT.EDIT/Ecoli5/fastqs/Seqfinder_results"
resultspath="/mnt/Molsig_VM0518/Miranda_Ecoli_Aug17/Seqfinder_output"
cores=1

#sh WGS_GeneMapper_Run.sh database_of_genes.fasta name_of_your_project $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh AMRDatabase_Colistin_mcr_variants_20170815.fna fastqs1 $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh AMRDatabase_20170815.fna fastqs $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh AMR-20180308_and_EnterobPlasmidsReps_20180309.fas fastqs $cores $softpath $refpath $datapath $resultspath
sh WGS_GeneMapper_Run.sh AMRDatabase_20180308.fna fastqs $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh Enterobacteriaceae_Plasmids_RepliconType_20180309.fas fastqs $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh Heavy_metals_180604.fas fastqs $cores $softpath $refpath $datapath $resultspath
