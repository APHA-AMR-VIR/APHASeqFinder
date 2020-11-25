#!/usr/bin/bash
echo 'Started running script at:'; date
echo 'Number of samples to run:'; ls *.csv | wc -l
for file in *Compare*.csv
do python ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.1/seqfinder/good_snps_08_2019.py $file
done
python ~/mnt/BactiPipes_MA2017/APHA_pipelines/SeqFinder_v1.1/seqfinder/making_a_compilation_for_seqfinder_output.py
