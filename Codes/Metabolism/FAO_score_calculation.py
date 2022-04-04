# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 17:34:26 2022
@author: Srinath
"""

#Script for FAO score calculation as described in Jia D, Paudel B, et.al. Front. Oncol. 2020
"""
input
    File containing gene expression matrix (with genes as rows and samples as columns)
    Gene signature file for FAO (Fatty acid oxidation) pathway (contains one column) 
return
    FAO scores for each sample written onto a text file 
"""
import os
import pandas as pd 

os.chdir('Datasets')
filelist = os.listdir() 
GSEID = [element.strip('.csv') for element in filelist] 
num = len(GSEID)
sigs_list = list(pd.read_csv('Gene signatures/FAO_geneset.txt', header='infer')['Gene']) 

for i in range(num):  
    data = pd.read_csv(GSEID[i]+'.csv', header='infer') 
    data = data.drop(data.columns[[1]], axis = 1) 
    temp = list(data.columns); temp[0] = 'Gene'; data.columns = temp  #First column labelled as 'Gene'
    data = data.loc[data['Gene'].isin(sigs_list)]
    data = data.set_index(keys='Gene')   #Column name in signature file must match that of first column in gene exp file, i.e, 'Gene'
    GSE = GSEID[i]  
    s = data.subtract(data.mean(axis=1), axis=0).truediv(data.std(axis=1, ddof=0), axis=0).sum()/14
    s.name = 'FAO Scores'
    faodf = pd.DataFrame(s).reset_index()
    faodf = faodf.rename(columns={'index': 'Sample'})
    faodf.to_csv(f'/Compiled Scores/Metabolism/FAO Output/{GSE}_faoscores.txt', index=False, sep='\t') 