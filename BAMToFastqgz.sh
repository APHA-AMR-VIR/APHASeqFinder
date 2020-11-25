
echo $1

cd $1

for filename in `ls *".bam"`
do
echo "file="$filename

of=`basename $filename ".bam"`

#java -Xmx6g -jar /media/Second3TB/Work/WorkVLA/Projects/NSOR1041/AMR_pipe/SamToFastq.jar INCLUDE_NON_PF_READS=true INPUT=$filename F=$of"_R1_001.fastq" F2=$of"_R2_001.fastq"

picard-tools SamToFastq I=$filename F=$of"_R1_001.fastq" F2=$of"_R2_001.fastq"

#samtools view -h -o $of".sam" $filename
#cat $of".sam" | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > $of"_R1_001.fastq"
#cat $of"_R1_001.fastq" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $of"_sorted_R1_001.fastq"

#cat $of".sam" | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > $of"_R2_001.fastq"
#cat $of"_R2_001.fastq" | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $of"_sorted_R2_001.fastq"

gzip $of"_R1_001.fastq"
gzip $of"_R2_001.fastq"

#rm *.sam

done


