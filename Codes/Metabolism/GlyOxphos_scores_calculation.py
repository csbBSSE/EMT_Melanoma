# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""

#Script for OXPHOS and glycolysis score calculation using Single Sample Gene Set Enrichment Analysis (SSGSEA) 
"""
input
    File containing gene expression matrix (with genes as rows and samples as columns)
    Gene signature file for OXPHOS and Glycolysis pathways (hallmark gene sets obtained from MSigDB) 
return
    For both pathways, normalised enrichment scores obtained from SSGSEA written onto a text file for each sample 
"""

import os
import pandas
import gseapy as gp
import matplotlib.pyplot as plt 

os.chdir('Datasets')
FileList = os.listdir() 
GSEID = [element.strip('.csv') for element in FileList] 
num = len(GSEID) 

for i in range(num):
    df = pandas.read_csv(GSEID[i] + '.csv', header = 'infer')  
    df = df.drop(df.columns[[1]], axis = 1)           
    GMT_file = 'Gene signatures/GSEA_geneset.gmt'  
    ss = gp.ssgsea(data= df,
         gene_sets= GMT_file,
         outdir = '/Compiled Scores/Metabolism/GSEA_Output/' + GSEID[i],  
         sample_norm_method='rank',
         no_plot = True, processes=4, format='png')
    
    