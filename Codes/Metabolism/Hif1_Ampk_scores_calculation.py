# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 21:32:32 2022
@author: Srinath
"""

#Script for HIF-1 and AMPK scores calculation as described in Yu L, Lu M, et.al. Cancer Res. 2017
"""
input
    File containing gene expression matrix (with genes as rows and samples as columns)
    Gene signature files for HIF-1 and AMPK enzyme activities (each file contains one column containing downstream targets of respective enzyme)   
return
    HIF-1 scores for each sample written onto a text file 
    AMPK scores for each sample written onto another text file
"""

import os 
import pandas
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from singscore.singscore import score

os.chdir('Datasets')
FileList = os.listdir() 
GSEID = [element.strip('.csv') for element in FileList] 
num = len(GSEID) 
sigs_list = ['AMPK_geneset.txt', 'HIF1_geneset.txt']

for i in range(num): 
    data = pandas.read_csv(GSEID[i]+'.csv', header='infer')
    data = data.drop(data.columns[[1]], axis = 1) 
    temp = list(data.columns); temp[0] = 'Gene'; data.columns = temp  #First column labelled as 'Gene'
    data = data.set_index(keys='Gene')   
    GSE = GSEID[i]   
    for j in range(2): 
        sigs = pandas.read_csv(open('../Gene signatures/' + sigs_list[j], 'r'), header = 'infer') 
        sigs = list(sigs['Gene'])   #Column name in signature file must match that of first column in gene exp file, i.e, 'Gene'
        scored_data_single = score(up_gene = sigs, sample=data, norm_method='theoretical', full_data=True)
        fh = open('/Compiled Scores/MEtabolism/Sing scores output/' + GSE + '_' + sigs_list[j], 'w')  
        fh.write(scored_data_single.to_csv(sep = '\t', header = True))
        fh.close()