#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 11:43:19 2019

@author: p312523
"""

import pandas as pd
import os
import sys
from pathlib import Path

print('Now combining all the abricate_seqfinder.csv files')
results_path=sys.argv[1]
reference_name=sys.argv[2]
list_of_files=[p for p in Path(results_path).rglob('*_abricate_seqfinder.csv')]
#list_of_files= glob.glob(results_path + "*_abricate_seqfinder.csv")
print(list_of_files)
df = pd.concat((pd.read_csv(f) for f in list_of_files))
df_final=df.drop(['real-len','mapped-len','other'],1)
outputname=os.path.join(results_path,results_path.split(os.sep)[-1]+"__"+reference_name+"__abricate_seqfinder_compilation.csv")
df_final.to_csv(outputname,index=False)
print('Done! Output written to abricate_seqfinder_compilation.csv')

