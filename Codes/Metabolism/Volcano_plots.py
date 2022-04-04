# -*- coding: utf-8 -*-
"""
Created on Tue May 18 16:58:38 2021
@author: 1
"""

#Script for generating a volcano plot between negative log to the base 10 of p-value and correlation coefficient (R-value) 
#Each dot in the plot represents a dataset that has a p-value and r-value associated with it
"""
input
    Text file containing r-value in the first column and p-value in the second column for all the datasets analysed
return
    Volcano plot for -log10(p-value) vs r-value. 
    #Datasets with positive correlation depicted in red color and those with negative correlation in blue color. 
    #Threshold p-values and r-values depicted using horizontal and vertical lines respectively
"""

import os
import matplotlib.pyplot as plt    
import math

os.chdir('Correlations/Metabolism')  
FileList = os.listdir()    
num = len(FileList) 
var_name = [element.split('_') for element in FileList] 
title_list = [(element2[1] + ' vs ' + element2[2]) for element2 in var_name]   
master_corr = []; master_pval = []    
for i in range(num):   
   file = FileList[i]   
   fh= open(file,'r')   
   lines = fh.readlines()   
   Corr =[((element.strip('\n')).split('\t'))[1] for element in lines]; Corr_list = [float(element1) for element1 in Corr]
   P_val =[((element.strip('\n')).split('\t'))[2] for element in lines]; Pval_list = [(math.log10(float(element1)) * -1) for element1 in P_val]
   master_corr.append(Corr_list); master_pval.append(Pval_list)         
   fh.close() 

#%%
fig = plt.figure(figsize=(12,15))    
rect = fig.patch   
ax1 = fig.add_subplot(4,3,1); ax2 = fig.add_subplot(4,3,2); ax3 = fig.add_subplot(4,3,3)
ax4 = fig.add_subplot(4,3,4); ax5 = fig.add_subplot(4,3,5); ax6 = fig.add_subplot(4,3,6)
ax7 = fig.add_subplot(4,3,7); ax8 = fig.add_subplot(4,3,8); ax9 = fig.add_subplot(4,3,9)
ax10 = fig.add_subplot(4,3,10); ax11 = fig.add_subplot(4,3,11); ax12 = fig.add_subplot(4,3,12) 
axes_list = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12]   
size = 12 

def volcano_plot(ax, x_axis1, y_axis1, x_axis2, y_axis2, x_axis3, y_axis3, title):
    ax.clear()
    ax.scatter(x_axis1, y_axis1, size, color='red', label = 'Positive')
    ax.scatter(x_axis2, y_axis2, size, color='blue', label = 'Negative') 
    ax.scatter(x_axis3, y_axis3, size, color='grey', label = 'No correlation')  
    ax.yaxis.label.set_color('k') 
    ax.xaxis.label.set_color('k')
    ax.axvline(x=-0.3, color = 'k')  
    ax.axvline(x=0.3, color = 'k')  
    yval = -1 * math.log10(0.05) 
    ax.axhline(y = yval, color = 'k')  
    ax.set_xlim([-1,1]) 
    ymax = max( [max(y_axis1), max(y_axis2), max(y_axis3)] )    
    ax.text(0.7, ymax/1.5, str(len(x_axis1)), color = 'k', fontstyle = 'italic', fontsize = 14, ha = 'right', va = 'top') 
    ax.text(-0.7, ymax/1.5, str(len(x_axis2)), color = 'k', fontstyle = 'italic', fontsize = 14, ha = 'left', va = 'top')  
    ax.set_xlabel('R value', fontsize = 13, fontstyle = 'italic') 
    ax.set_ylabel('-log10 (P value)', fontsize = 13, fontstyle = 'italic')   
    ax.set_title(title, fontsize = 15, fontstyle = 'italic')   
    #ax.legend(loc = 'best')     
     
for i in range(num): 
    ax = axes_list[i] 
    x = master_corr[i]; y = master_pval[i] 
    title = title_list[i]
    x_axis1 = []; x_axis2 = []; x_axis3 = [] 
    y_axis1 = []; y_axis2 = []; y_axis3 = [] 
    for j in range(len(x)):
        if y[j]> (-1 * math.log10(0.05)) and x[j]> 0.3:
            y_axis1.append(y[j]); x_axis1.append(x[j])
        elif y[j]> (-1 * math.log10(0.05)) and x[j]< -0.3: 
            y_axis2.append(y[j]); x_axis2.append(x[j])   
        else:
            y_axis3.append(y[j]); x_axis3.append(x[j])
            
    volcano_plot(ax, x_axis1, y_axis1, x_axis2, y_axis2, x_axis3, y_axis3, title)  
fig.subplots_adjust(hspace = 0.65, wspace = 0.4) 
fig.savefig('Figures/Fig.5/PI_metabolism_volcano_plot_new.png', bbox_inches='tight')   
fig.clf()         