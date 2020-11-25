# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
#
#
# This bash script will run the AMR gene analysis for any number of squencing runs.
#
# go to:
# 	cd /home/nick/software/WGS_GeneMapper
# and run:
#	sh WGS_GeneMapper_mini_campy.sh

# Every line calls WGS_GeneMapper_Run providing the following arguments:
# 	arg1 the reference fasta database containing the AMR genes sequences in fasta format
# 	arg2 the name of the run (folder name where the fastq files are stored).
# 	arg3 number of cpus to be used
#
#sleep 2h

softpath="/home/nick/software/WGS_GeneMapper"
refpath="/home/nick/software/references"
datapath="/home/nick/WGS_Data"
resultspath="/home/nick/WGS_results"
cores=2

sh WGS_GeneMapper_Run.sh AMR_genes_20150612.fna broiler_salmonella $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh BacMet_metals_nucleotide.fna ESBL $cores $softpath $refpath $datapath $resultspath
