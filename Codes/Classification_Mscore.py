#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""
Created on Mon Nov 1 2021
Last Edited on Fri Apr 1 2022

@author: Nilay

Description: Counting Cell Lines per Melonoma Phenotype based on chosen Range and M/E-Score. Written as a part of the
project 'Mapping phenotypic heterogeneity in melanoma onto the epithelial-hybrid-mesenchymal axis'

"""


# In[2]:


#Importing Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# In[3]:


'''
Setting some arrays for ease in sorting and presenting data later. These are tailored based on the data provided
since the data format provided is different for different datasets
'''
Files = ['GSE4843', 'GSE7127', 'GSE81383', 'GSE80829', 'GSE10916', 'GSE115978', 'GSE134432']
Phenotypes = ['Proliferative', 'Transitory', 'NCSC', 'Invasive'] #In Cluster Order
Score = ['KS_M', 'KS_E']
Range = ['Bottom', 'Intermediate', 'Top']


# In[4]:


'''
NOTE:
The following code is tailored around the provided 8 datasets, and the format in which they were provided.
It is possible to increase efficiency using os and glob libraries, for much larger number of datasets in a 
suitable format.
'''


# In[5]:


#User Input
ScoreType = Score[0] #Choose 0: M-Score 1: E-Score
RangeType = Range[0] #Choose 0: Bottom Third 1: Intermediate 2: Top Third


# In[6]:


#Sorting the scores as per requirement
for FileNum in range(len(Files)):
    
    #Choosing the 8 datasets from a large number. Change the path based on where the files are being stored.
    if FileNum == 0:
        Data = pd.read_csv(r'Compiled Scores\GSE4843.csv')
        Clust = pd.read_csv(r'Clusters\4843_clusters.txt', sep='\t')
    
    if FileNum == 1:
        Data = pd.read_csv(r'Compiled Scores\GSE7127.csv')
        Clust = pd.read_csv(r'Clusters\7127_clusters.txt', sep='\t')
    
    if FileNum == 2:
        Data = pd.read_csv(r'Compiled Scores\GSE81383.csv')
        Clust = pd.read_csv(r'Clusters\81383_clusters.txt', sep='\t')
    
    if FileNum == 3:
        Data = pd.read_csv(r'Compiled Scores\GSE80829.csv')
        Clust = pd.read_csv(r'Clusters\80829_clusters.txt', sep='\t')
    
    if FileNum == 4:
        Data = pd.read_csv(r'Compiled Scores\GSE10916.csv')
        Clust = pd.read_csv(r'Clusters\10916_clusters.txt', sep='\t')
    
    if FileNum == 5:
        Data = pd.read_csv(r'Compiled Scores\GSE115978.csv')
        Clust = pd.read_csv(r'Clusters\115978_clusters.txt', sep='\t')
        
    if FileNum == 6:
        Data = pd.read_csv(r'Compiled Scores\GSE134432.csv')
        Clust = pd.read_csv(r'Clusters\134432_clusters.txt', sep='\t')
    
    #Sorting
    Data[ScoreType] = (Data[ScoreType] - np.mean(Data[ScoreType]))/np.std(Data[ScoreType])
    minimum = min(Data[ScoreType])
    maximum = max(Data[ScoreType])
    lim = (maximum-minimum)/3
    
    #Some processing to account for cluster file format
    if FileNum < 6:
        Clust.columns = ['Clusters']
    else:
        Clust.columns = ['Val', 'Clusters'] #Because of an extra column
        
    Clust['Index'] = range(len(Data))
    Clust.set_index('Index', inplace=True, drop = True)
    Data['Cluster'] = Clust['Clusters']

    InCounts = np.empty(len(Phenotypes))
    OutCounts = np.empty(len(Phenotypes))
    FracIn = np.empty(len(Phenotypes))
    for i in range(len(Phenotypes)):
        
        if RangeType == 'Bottom':
            TempCount = Data[(Data[ScoreType] < (minimum + lim)) & (Data.Cluster== i + 1)]
        elif RangeType == 'Intermediate':
            TempCount = Data[((Data[ScoreType] <= (maximum - lim)) & (Data[ScoreType] >= (minimum + lim))) & (Data.Cluster== i + 1)]
        elif RangeType == 'Top':
            TempCount = Data[(Data[ScoreType] > (maximum - lim)) & (Data.Cluster== i + 1)]
        else:
            print('Error: Please check if you have selected a Range Type')
            
        InCounts[i] = len(TempCount)
        OutCounts[i] = len(Data[Data.Cluster == i + 1]) - InCounts[i]
        FracIn[i] = InCounts[i]/(OutCounts[i]+InCounts[i])
        
    print(Files[FileNum])
    print(ScoreType)
    print(RangeType)
    print(Phenotypes)
    print('Total cell lines inside the chosen range')
    print(InCounts)
    print('Total cell lines outside the chosen range')
    print(OutCounts)
    print('Fraction of cell lines inside')
    print(FracIn)
    print('\n')


# In[7]:


#GSE116237 needs to be done seperately because its data format is not consistent with the rest

Data = pd.read_csv(r'Compiled Scores\GSE116237.csv')
Clust = pd.read_csv(r'Clusters\116237_clusters.txt', sep='\t')

Data[ScoreType] = (Data[ScoreType] - np.mean(Data[ScoreType]))/np.std(Data[ScoreType])
minimum = min(Data[ScoreType])
maximum = max(Data[ScoreType])
lim = (maximum-minimum)/3

Clust['Index'] = range(len(Data))
Clust['Index2'] = range(len(Data))
Clust.set_index('Index', inplace=True, drop = True)
Clust = Clust[(Clust.x == '1') | (Clust.x == '2') | (Clust.x == '3') | (Clust.x == '4')]
Data = Data[Data.index.isin(Clust['Index2'])]
Data['Cluster'] = Clust['x'].astype(int)

InCounts = np.empty(len(Phenotypes))
OutCounts = np.empty(len(Phenotypes))
FracIn = np.empty(len(Phenotypes))
for i in range(len(Phenotypes)):
        
    if RangeType == 'Bottom':
        TempCount = Data[(Data[ScoreType] < (minimum + lim)) & (Data.Cluster== i + 1)]
    elif RangeType == 'Intermediate':
        TempCount = Data[((Data[ScoreType] <= (maximum - lim)) & (Data[ScoreType] >= (minimum + lim))) & (Data.Cluster== i + 1)]
    elif RangeType == 'Top':
        TempCount = Data[(Data[ScoreType] > (maximum - lim)) & (Data.Cluster== i + 1)]
    else:
        print('Error: Please check if you have selected a Range Type')
            
    InCounts[i] = len(TempCount)
    OutCounts[i] = len(Data[Data.Cluster == i + 1]) - InCounts[i]
    FracIn[i] = InCounts[i]/(OutCounts[i]+InCounts[i])
        
print('GSE116237')
print(ScoreType)
print(RangeType)
print(Phenotypes)
print('Total cell lines inside the chosen range')
print(InCounts)
print('Total cell lines outside the chosen range')
print(OutCounts)
print('Fraction of cell lines inside')
print(FracIn)


# In[ ]:




