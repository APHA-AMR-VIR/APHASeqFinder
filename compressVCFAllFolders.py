import os,sys,re,csv,itertools
import os, fnmatch

patho="/media/Second3TB/Work/WorkVLA/Data/MRSA/Results/WGS_Results"
words=[".pileup.vcf","_alignment_stats.csv"]

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

for word in words:
    print "Doing word: "+word
    for fil in find('*'+word, patho):
        print "Compressing file: "+fil
        os.system("gzip "+fil)