#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 08:16:23 2019

Investigate ancestry in SPARK dataset from PCA results.

@author: nbaya
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/'

## View results from Ricopili PCA
df = pd.read_csv(wd+'pca_preimp3.menv.mds',delim_whitespace=True) #NOTE: number of samples is 7173, 37 fewer than ADMIXTURE results

df['anc_true'] = df['FID'].str.split('_',expand=True)[3] #decided to call it "anc" not pop, because pop is a command in pandas
df.loc[df['anc_true'] == 'mix', 'anc_true'] = 'unknown'

ancestry_ls = ['unknown','eur','afr','amr','asn']

fig,ax=plt.subplots(figsize=(6*1.2,4*1.2))
for anc in ancestry_ls:
    plt.scatter(df[df.anc_true==anc].C1,df[df.anc_true==anc].C2,alpha=0.5,s=5,marker='.')
plt.legend(ancestry_ls)


## View results from ADMIXTURE
admix_Q = pd.read_csv(wd+'merged_v2.4.Q',delim_whitespace=True,header=None,names=[f'Q{i}' for i in range(4)]) # ancestry fractions
admix_P = pd.read_csv(wd+'merged_v2.4.P',delim_whitespace=True,header=None,names=[f'Q{i}' for i in range(4)]) # allele frequencies
'''
Pop0 = African (afr)
Pop1 = Native American (amr)
Pop2 = Asian (asn)
Pop3 = European (eur)
'''
fam = pd.read_csv(wd+'merged_v2.fam',delim_whitespace=True,header=None,names=['FID','IID','PAT','MAT','SEX','PHEN'])
fam = fam.join(admix_Q)


#fam[fam['anc_true']=='amr'][[f'Q{i}' for i in range(4)]]



##Combine PCA results with ADMIXTURE results
merged = df.merge(fam,on=['FID','IID']) #inner join, therefore limited to the sample size of df, 7173 samples remaining, 37 removed from fam

merged['anc_est'] = 'unknown'
min_Q = 0.5
merged.loc[merged.Q0 > min_Q, 'anc_est'] = 'afr' #ancestry estimate
merged.loc[merged.Q1 > min_Q, 'anc_est'] = 'amr'
merged.loc[merged.Q2 > min_Q, 'anc_est'] = 'asn'
merged.loc[merged.Q3 > min_Q, 'anc_est'] = 'eur'

print('\n'.join([f'Individuals of {x} ancestry: {len(merged[merged.anc_est==x])} ({round(len(merged[merged.anc_est==x])/len(merged)*100,2)}%)' for x in ancestry_ls]))
print(f'Total: {len(merged)}')

spark = merged[merged.anc_true=='unknown']
print('\n'.join([f'Individuals of {x} ancestry: {len(spark[spark.anc_est==x])} ({round(len(spark[spark.anc_est==x])/len(spark)*100,2)}%)' for x in ancestry_ls]))
print(f'Total: {len(spark)}')


ancestry_ls = ['unknown','eur','afr','amr','asn']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

pcs = [1,2]
fig,ax=plt.subplots(figsize=(6*1.5,4*1.5))
for anc_idx, anc in enumerate(ancestry_ls):
    df_tmp1 = merged[(merged.anc_est==anc)&(merged.anc_true=='unknown')] #spark individuals
    plt.scatter(df_tmp1[f'C{pcs[0]}'],df_tmp1[f'C{pcs[1]}'],alpha=0.5,s=10,marker='o', c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
for anc_idx, anc in enumerate(ancestry_ls): 
    df_tmp2 = merged[(merged.anc_est==anc)&(merged.anc_true!='unknown')] #reference panel individuals
    plt.scatter(df_tmp2[f'C{pcs[0]}'],df_tmp2[f'C{pcs[1]}'],alpha=0.5,s=20,marker='x', c = colors[anc_idx])
legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]+
                    [Line2D([0],[0],lw=0,markerfacecolor='k',marker='o',label='SPARK',markeredgecolor='w'),
                     Line2D([0],[0],lw=0,markeredgecolor='k',marker='x',label='ref')])
#legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
ax.legend(handles =legend_elements)
plt.xlabel(f'PC{pcs[0]}')
plt.ylabel(f'PC{pcs[1]}')



#single out one ancestry
ancestry = 'unknown' #options: eur, afr, amr, asn, unknown
pcs = [1,2]
fig,ax=plt.subplots(figsize=(6*1.5,4*1.5))
sorting_vals = ([0]*len(ancestry_ls))
sorting_vals[ancestry_ls.index(ancestry)] = 1
ancestry_sort = list(zip(ancestry_ls,sorting_vals))
ancestry_sort = sorted(ancestry_sort, key=lambda x: x[1])
for anc in [x[0] for x in ancestry_sort]:
    df_tmp1 = merged[(merged.anc_est==anc)&(merged.anc_true=='unknown')] #spark individuals
    plt.scatter(df_tmp1[f'C{pcs[0]}'],df_tmp1[f'C{pcs[1]}'],alpha=(0.5 if anc==ancestry else 1),s=10,marker='o', c = (colors[ancestry_ls.index(anc)] if anc==ancestry else '#d8dcd6'))#,edgecolors='k',linewidths=0.1)
    df_tmp2 = merged[(merged.anc_est==anc)&(merged.anc_true!='unknown')] #reference panel individuals
    plt.scatter(df_tmp2[f'C{pcs[0]}'],df_tmp2[f'C{pcs[1]}'],alpha=(0.5 if anc==ancestry else 1),s=20,marker='x', c = (colors[ancestry_ls.index(anc)] if anc==ancestry else '#d8dcd6'))
    
legend_elements = [Line2D([0],[0],lw=0,markerfacecolor=colors[ancestry_ls.index(ancestry)],marker='o',label=ancestry+' (SPARK)',markeredgecolor='w'),
                     Line2D([0],[0],lw=0,markeredgecolor=colors[ancestry_ls.index(ancestry)],marker='x',label=ancestry+' (ref)'),
                     Patch(facecolor='#d8dcd6',label='other')]
ax.legend(handles =legend_elements,title='Ancestry')
plt.xlabel(f'PC{pcs[0]}')
plt.ylabel(f'PC{pcs[1]}')
