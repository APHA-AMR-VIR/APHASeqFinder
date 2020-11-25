# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
#
# need to have pandas    'conda install pandas==0.23.4'
# 'conda install -c conda-forge biopython'
#
# This bash script will run the AMR gene analysis for any number of squencing runs.
#
# go to:
# 	cd ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.2/seqfinder/
# and run:
#	sh seqfinder_MA1.2.sh path_to_fastq_files
#	sh seqfinder_MA1.2.sh ~/mnt/Non-Defra_Projects/EU_Projects/IMPART_2018/RAW_WGS/Staphylococcus_hyicus/ 

# Every line calls WGS_GeneMapper_Run providing the following arguments:
# 	arg1 the reference fasta database containing the AMR genes sequences in fasta format
# 	arg2 the name of the run (folder name where the fastq files are stored).
# 	arg3 number of cpus to be used
#
#sleep 24h

softpath=$(pwd)
refpath=$(pwd)/references
systemcores=$(grep -c ^processor /proc/cpuinfo)
reference_key=$(cat reference_key.txt)
amr=$(echo $(ls -al references/AMR | grep '^-' | awk '{print $9}'))
amrpath=$(pwd)/references/AMR/
vir=$(echo $(ls -al references/virulence | grep '^-' | awk '{print $9}'))
virpath=$(pwd)/references/virulence/
metal=$(echo $(ls -al references/metal | grep '^-' | awk '{print $9}'))
metalpath=$(pwd)/references/metal/
#echo "Do you want to download the lastest databases from github y/n?"
#read -t 30 latest_refs
#if [[ $latest_refs == "y" ]]; then
#    wget https://github.com/PyroEndy/APHASeqFinder/raw/master/references.tar.gz references.tar.gz && tar -xf references.tar.gz && rm #references.tar.gz*;elif [[ $latest_refs == "n" ]]; then echo "Using pre-downloaded databases..."
#fi
#if [[ ! $latest_refs =~ ^(y|n)$ ]]; then echo "'$latest_refs' is not a valid response. Please choose y or n."; exit 1
#fi
echo "You have $systemcores cores, how many do you want to use?"
read -t 30 varcore
if [[ $varcore == "" ]]; then
    echo "Exiting because number of cores not chosen!"; exit 1
fi
if [[ $varcore == 0 ]]; then
    echo "Exiting because number of cores is zero!"; exit 1
fi
if [[ -n ${varcore//[0-9]/} ]]; then
    echo "Exiting because '$varcore' contains non-numerical characters!"; exit 1
fi
#if [[ $varcore > $systemcores ]]; then
#    echo "Exiting because you've chosen more cores than you have!"; exit 1
#fi
echo "You have chosen $varcore cores, so I will try and run up to $varcore samples in parallel"
echo "Pick the reference you want to run from this list:
$reference_key"
read -t 30 reference_pick
if [[ $reference_pick == "1" ]]; then
reference=$amr 
refpath=$amrpath
elif [[ $reference_pick == "2" ]]; then
reference=$vir
refpath=$virpath
elif [[ $reference_pick == "3" ]]; then
reference=$metal
refpath=$metalpath
elif [[ $reference_pick == "4" ]]; then
echo "Pick from this reference list"
ls -p $refpath | grep -v /
read -t 30 other
reference=$other
elif [[ $reference_pick > 4 ]]; then echo "Exiting because reference number not in the list!"; exit 1; elif [[ -n ${reference_pick//[0-9]/} ]]; then echo "Exiting because input contained non-numerical characters!"; exit 1; elif [[ $reference_pick == "" ]]; then echo "Exiting because database not chosen"; exit 1
fi
echo "What is the minimum percent mapped to a gene you want to use for filtering? Recommended >95%"
read -t 30 PerID
if [[ $PerID == "" ]]; then
    echo "Exiting because percent ID not chosen!"; exit 1
fi
if [[ $PerID == 0 ]]; then
    echo "Exiting because percent ID is zero!"; exit 1
fi
if [[ -n ${PerID//[0-9]/} ]]; then
    echo "Exiting because '$PerID' contains non-numerical characters!"; exit 1
fi
echo "What is the maximum number of SNPs you want in a gene to still be called present? Recommended <100"
read -t 30 NumSNPs
if [[ $NumSNPs == "" ]]; then
    echo "Exiting because number of number of SNPs not chosen!"; exit 1
fi
if [[ -n ${NumSNPs//[0-9]/} ]]; then
    echo "Exiting because '$NumSNPs' contains non-numerical characters!"; exit 1
fi
echo "Running your samples with $reference as your reference"
#datapath= $1
#resultspath="/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1/nick/WGS_results"
cores=$varcore
sh WGS_GeneMapper_Run.sh $reference seqfinder $cores $softpath $refpath $1 $1 $PerID $NumSNPs
