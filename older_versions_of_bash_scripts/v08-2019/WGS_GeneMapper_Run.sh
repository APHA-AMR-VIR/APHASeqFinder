# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
#
#
# This bash script will run the AMR gene analysis for a squencing run.
#
# go to:
# 	cd /media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMR_pipe
# and run:
#	sh AMR_run.sh "databasename.fna" "datafolder" "# cpu"

#	The arguments needed are:
# 	$1 the reference fasta database containing the AMR genes sequences in fasta format
# 	$2 the name of the run (folder name where the fastq files are stored).
# 	$3 number of cpus to be used


# The following are variables of the location of the files involved in the analysis.
# This only need to be changed if the pipeline is installed on another location/computer. 

echo "Running sequencing run "$2" for database "$1

refname=$(echo $1 | cut -d'.' -f1)

softpath=$4
refpath=$5
datapath=$6 
resultspath=$7"/"$refname

# The results folder is created
mkdir -p $resultspath
mkdir $resultspath/filtered_results
#echo "Results stored at "$resultspath
#echo "Data path is "$datapath
#echo "Softpath is "$softpath
#echo "refpath is "$refpath
# The statistics document is initiated with the name of the columns
statsFile=$2"_Mapped_"$refname"_SamplesInfo.txt"
echo -n "fileName	runName	fastq1_reads	med	fastq2_reads	med	total	f_fastq1_reads	med	f_fastq2_reads	med	total	t_fastq1_reads	med	t_fastq2_reads	med	total	mapped_reads	percentage	med\n" > $resultspath/$statsFile 


cd $datapath
#echo "Data path is "$datapath

# For each sample in the data folder bash script AMR_sample.sh is called with the following arguments:
# arg1 software path
# arg2 reference database
# arg3 data path
# arg4 stats file name
# arg5 results path
# arg6 reference path
# arg7 run name

# if $3 (number of cpus) is an empty, the execution is sequential otherwise is in parallel 
echo "Using $3 cores"
#echo $3
if [ "$3" = "1" ]; then
	for file in `ls *1_001.fastq.gz`
	do
	echo $file
	sh $softpath/WGS_GeneMapper_Sample.sh $softpath $1 $datapath $statsFile $resultspath $refpath $2 $file  
	done
else
	ls *1_001.fastq.gz  | parallel -j $3 sh $softpath/WGS_GeneMapper_Sample.sh $softpath $1 $datapath $statsFile $resultspath $refpath $2 {}
fi

# The python script "CompileResultsFRomGeneLists.py" is called to compile the result tables for all the samples in the run
# into a single table for a given column.  

#python $softpath/CompileResultsFRomGeneLists.py $resultspath '_CompareTo_' '% mapped-gaps/real' $7
#python $softpath/CompileResultsFRomGeneLists.py $resultspath '_CompareTo_' 'goodSnps' $7
#python $softpath/CompileResultsFRomGeneLists.py $resultspath '_CompareTo_' '#syn' $7
#python $softpath/CompileResultsFRomGeneLists.py $resultspath '_CompareTo_' '#non' $7
cp $resultspath/*/*Compare* $resultspath/filtered_results
cd $resultspath/filtered_results && sh $softpath/APHA_condense.sh
mkdir $resultspath/filtered_results/original_output
mv $resultspath/filtered_results/*Compare* $resultspath/filtered_results/original_output
