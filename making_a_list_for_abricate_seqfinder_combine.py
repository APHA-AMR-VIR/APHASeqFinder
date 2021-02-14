#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:49:15 2019

@author: p312523
"""
import pandas as pd
import fnmatch
import os

path ="./"
list_of_seqfinder=[x for x in os.listdir(path) if fnmatch.fnmatch(x, '*good_snps.csv')]

list_of_strains=[]
list_of_seqfinder_df=pd.DataFrame(list_of_seqfinder)
for i in list_of_seqfinder:
    test=i.split("_CompareTo_")
    list_of_strains.append(test[0])
list_of_seqfinder_df["strain"]=list_of_strains
list_of_abricate_files=[]
for i in list_of_seqfinder_df["strain"]:
    list_of_abricate=[x for x in os.listdir(path) if fnmatch.fnmatch(x, i+'.abricate')]
    list_of_abricate_files.append(list_of_abricate)
list_of_seqfinder_df['abricate']=list_of_abricate_files
list_of_seqfinder_df['abricate']=list_of_seqfinder_df['abricate'].map(lambda x: str(x)[2:-2])
list_of_seqfinder_df=list_of_seqfinder_df.drop(['strain'], axis=1)
list_of_seqfinder_df.to_csv('./samples.csv', header=False, index=False)

