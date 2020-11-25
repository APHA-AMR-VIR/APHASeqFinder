#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 11:43:19 2019

@author: p312523
"""

import pandas as pd
import glob
print('Now combining all the good_snps_chr.csv files')
path = r"./"
list_of_files= glob.glob(path + "*good_snps_only_chromosomal.csv")
df = pd.concat((pd.read_csv(f) for f in list_of_files))
df_final=df.drop(['real-len','other'],1)
outputname="seqfinder_chr_compilation.csv"
df_final.to_csv(outputname,index=False)
print('Done! Output written to seqfinder_chr_compilation.csv')
