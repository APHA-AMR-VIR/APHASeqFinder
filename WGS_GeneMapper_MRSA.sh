# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
#
#
# This bash script will run the AMR gene analysis for any number of squencing runs.
#
# go to:
# 	cd /media/Second3TB/Work/WorkVLA/Software/WGS_GeneMapper
# and run:
#	sh WGS_GeneMapper_MRSA.sh

# Every line calls WGS_GeneMapper_Run providing the following arguments:
# 	arg1 the reference fasta database containing the AMR genes sequences in fasta format
# 	arg2 the name of the run (folder name where the fastq files are stored).
# 	arg3 number of cpus to be used

#sleep 30h

softpath="/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1/software/WGS_GeneMapper"
refpath="/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1/software/references"
datapath="/home/p175776/zinc_test"
resultspath="/home/p175776/WGS_results"

cores=7

#sh WGS_GeneMapper_Run.sh GeneList_20151014.fas APHA $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh GeneList_20151014.fas APHA_2 $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh GeneList_20151014.fas Belgium $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh GeneList_20151014.fas Bruno $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh GeneList_20151014.fas NC_017333 $cores $softpath $refpath $datapath $resultspath
#sh WGS_GeneMapper_Run.sh GeneList_20151014.fas PHE $cores $softpath $refpath $datapath $resultspath
sh WGS_GeneMapper_Run.sh staph_czr.fna MRSA $cores $softpath $refpath $datapath $resultspath 



