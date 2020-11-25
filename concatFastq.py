#python $HOME/Desktop/Software/joinFastq.py $HOME/Desktop/data/20131023_HiSeq

from os import listdir
import os, sys, os.path, csv, operator
from operator import itemgetter
import subprocess, os
from subprocess import check_call

args=sys.argv
path = sys.argv[1]

samples =[]
fnames = listdir(path)
for fname in fnames:
  if fname.split("_")[0] not in samples:
  	samples=[fname.split("_")[0]]+samples
path=path+'/'
for sample in samples:
	print 'Concatenating files for sample ' + sample
	for Rn in ['R1','R2']:
   		fends = [x for x in fnames if sample==x.split("_")[0] and Rn+'_001.fastq.gz' in x]
                fends.sort()
   		if len(fends) > 1 :
			print fends
   			for fend in fends:
   				cmd = 'gunzip -c '+ path + fend + ' > ' + path + fend.split(".")[0] +'.fastq'
   				check_call(cmd,shell=True)
   			cats = ' '.join([path + x[:-3] for x in fends])
   			cmd = 'cat ' + cats + ' > ' + path + sample + '_' + Rn + '_001.fastq'
   			check_call(cmd,shell=True)
   			cmd = 'gzip ' + path + sample + '_' + Rn + '_001.fastq'
   			check_call(cmd,shell=True)
			check_call('rm '+ path + '*.fastq',shell=True)
			for fname in fends:
 				os.remove(path + fname)

