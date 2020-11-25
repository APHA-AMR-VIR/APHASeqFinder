# June 2015, Central Sequencing Unit, APHA
# Author: Javier Nunez
#
#
# This bash script will run the AMR gene analysis for a squencing run.
#
# go to:
# 	cd /media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMR_pipe
# and run:
# 	sh AMR_sample.sh $1 $2 $3 $4 $5 $6 $7 $8 
# The arguments are as follow:
#
# 	$1 software path
# 	$2 reference database
# 	$3 data path
# 	$4 stats file name
# 	$5 results path
# 	$6 reference path
# 	$7 run name
# 	$8 name of the R1 fastq file

echo "Starting sample processing"
softpath=$1

of=`basename $8 "1_001.fastq.gz"` # of is the common name for both pair fastq files


nf=$(echo $8 | cut -d '_' -f1) # nf is the name of the sample

isoldir=$5/$nf # isolate results folder
mkdir $isoldir


cp $6/$2 $isoldir/$2 # copy the database into the isolate folder for future usage

cd $isoldir

# unzip fastq.gz
gunzip -c $3/$of"1_001.fastq.gz" > $nf"_R1.fastq"
gunzip -c $3/$of"2_001.fastq.gz" > $nf"_R2.fastq"

mv $nf"_R1.fastq" $nf".nodup_R1.fastq"
mv $nf"_R2.fastq" $nf".nodup_R2.fastq"
#rm $nf"_R1.fastq"
#rm $nf"_R2.fastq"
# The python script "duprm10.py" removes repeated reads (both pairs at the same time)
#python $1/duprm10.py $nf"_R1.fastq" $nf"_R2.fastq" $nf".nodup_R1.fastq" $nf".nodup_R2.fastq" $nf".dup_R1.fastq" $nf".dup_R2.fastq"

# Trimming the reads 
java -jar $1/trimmomatic-0.30.jar SE -phred33 $nf".nodup_R1.fastq" $nf".TrimIn_R1.fastq" SLIDINGWINDOW:10:20 MINLEN:80

java -jar $1/trimmomatic-0.30.jar SE -phred33 $nf".nodup_R2.fastq" $nf".TrimIn_R2.fastq" SLIDINGWINDOW:10:20 MINLEN:80

rm $nf.nodup_R*

cat $nf.TrimIn_R1.fastq $nf.TrimIn_R2.fastq > $nf.TrimIn.fastq
rm $nf.TrimIn_R*.fastq
# mapping step

echo "&&&&&&&&&&&&&&&& Mapping sample "$nf" to reference "$2

#echo $softpath/smalt

$softpath/smalt index -k 13 -s 6 $2 $2
$softpath/smalt map -d -1 -y 0.7 -x -f samsoft -o $nf".sam" $2  $nf".TrimIn.fastq" 
rm $nf.TrimIn.fastq
$softpath/samtools view -Sh -F 4 $nf".sam" > $nf".F4.sam"
rm $nf.sam
#perl $1/xa2multi.pl $nf".F4.sam" > $nf".xa.sam" 

python2 $1/SamFlagChanged.py $nf".F4.sam" $nf".60.sam" 60
rm $nf".F4.sam"
echo "&&&&&&&&&&&&&&&&&&&&&& samtools sam to bam"

$softpath/samtools view -Shu $nf".60.sam" > $nf".bam"
rm "$nf.60.sam"
echo "&&&&&&&&&&&&&&&& samtools sorting bam"

$softpath/samtools sort $nf".bam" $nf".sorted"
rm $nf.bam
#java -Xmx8g -jar $1/SamToFastq.jar INPUT=$nf".sorted.bam" F=$nf".sorted.fastq"
bedtools bamtofastq -i $nf".sorted.bam" -fq $nf".sorted.fastq"
#rm $isoldir/$nf".sorted.fastq"
echo "&&&&&&&&&&&&&&&& Pileup bcf"

$softpath/samtools index $nf".sorted.bam"
$softpath/samtools faidx $2
$softpath/samtools mpileup -q -1 -uf $2 $nf".sorted.bam" > $nf".pileup.bcf" 

########## consensus sequence
echo "&&&&&&&&&&&&&&&& Calculating VCF and fq consensus"

$softpath/bcftools view -cg $nf".pileup.bcf" > $nf".pileup.vcf"
perl $softpath/vcfutils.pl vcf2fq $nf".pileup.vcf" > $nf".fq"

########## Alignment stats
# 	arg1 $nf or name of the strain
# 	arg2 or $7 is the run name
# 	arg3 or $2 is the reference database
# 	arg4 or $4 is the stats file name
# 	arg5 or $of is the base name of the mates fastq files of the sample
# 	arg6 or $isoldir is the directory of the sample
# 	arg7 or $5 is the results path

echo "&&&&&&&&&&&&&&&& Collecting information related to the mapping process from the SAM and VCF files"
echo "run name is $7"
python2 $1/SMA-VCF_Stats.py $nf $7 $2 $4 $of $isoldir $5

echo "&&&&&&&&&&&&&&&& Applying snp filter min of 180 quality, max of 0.05 ratio wrong good coverage"

########## variant calling
$softpath/bcftools view -vcg $nf".pileup.bcf" > $nf".snp"
perl $softpath/vcfutils.pl varFilter $nf".snp" > $nf".flt.snp"
python2 $1/snpsFilter.py 150 1 2 $nf".snp" $nf"_SN.csv" #180 0.5
rm $nf".pileup.bcf"
echo "&&&&&&&&&&&&&&&& AMR Genes"

python2 $1/calculatePresentGenes3.py $nf".fq" $2 $nf"_alignment_stats.csv" 150 0.2 2 yes

# FASTQC for quality control
#for file in `ls *.fastq`
#do
#echo "Doing FASTQC for "$file 
#java -Xmx250m -classpath $1/FastQC:$1/FastQC/sam-1.32.jar:$1/FastQC/jbzip2-0.9.jar -Dfastqc.output_dir=$isoldir uk.ac.babraham.FastQC.FastQCApplication $file
#done

gzip $nf".pileup.vcf"
gzip $nf"_alignment_stats.csv"
#rm $nf".pileup.vcf"
#rm $nf"_alignment_stats.csv"
rm *fastq
rm *snp
rm *smi
rm *sma
tar -zcvf accessory_files.tgz ./ --exclude=./*Compare*
rm *.bam
rm *.gz
rm *goodsnps.csv
rm *SN.csv
rm *Deleted*
rm *.fai
rm *.bai
rm *.f*


