#!/usr/bin/bash
echo 'Started running script at:'; date
echo 'Number of samples to run:'; ls *.csv | wc -l
for file in *Compare*.csv
do python2 ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.2/seqfinder/good_snps_12_2019.py $file $1 $2
done
python2 ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.2/seqfinder/making_a_compilation_for_seqfinder_output.py
python2 ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.2/seqfinder/making_a_Chr_compilation_for_seqfinder_output.py 

