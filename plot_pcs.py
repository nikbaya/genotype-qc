#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 08:16:23 2019

Investigate ancestry in SPARK dataset from PCA results.

@author: nbaya
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from datetime import datetime as dt

spark_wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/'
preimp7_wd = '/Users/nbaya//Documents/lab/genotype-qc/spark/preimp7_imus/'


## View results from Ricopili PCA
#df = pd.read_csv(spark_wd+'preimp3/SPARK.parents.IMUS.menv.mds',delim_whitespace=True) #NOTE: number of samples is 6082, 1128 fewer than ADMIXTURE results (compared to 7173 from the previous non-liftover, autosomes+non-autosomes run)
df = pd.read_csv(spark_wd+'preimp3/SPARK.parents.IMUS_v2.menv.mds',delim_whitespace=True) #NOTE: Includes MAF filter >0.5%, number of samples is 6081, 1129 fewer than ADMIXTURE results (compared to 7173 from the previous non-liftover, autosomes+non-autosomes run)

len(df)

df['anc_true'] = df['FID'].str.split('_',expand=True)[3] #decided to call it "anc" not pop, because pop is a command in pandas
df.loc[df.FID.str.contains('mix'), 'anc_true'] = 'unknown'

len(df[df.anc_true!='unknown'])
len(df[df.anc_true=='unknown'])

ancestry_ls = ['unknown','eur','afr','amr','asn']

fig,ax=plt.subplots(figsize=(6*1.5,4*1.5))
for anc in ancestry_ls:
    plt.scatter(df[df.anc_true==anc].C1,df[df.anc_true==anc].C2,alpha=0.5,s=5,marker='.')
plt.legend(ancestry_ls)


fam = pd.read_csv(spark_wd+'preimp3/SPARK.parents.IMUS_v2.menv.fam',delim_whitespace=True,header=None,names=['FID','IID','PAT','MAT','SEX','PHEN'])

df = df.merge(fam, on=['FID','IID'])


## View results from ADMIXTURE
admix_Q = pd.read_csv(spark_wd+'merged_v2.4.Q',delim_whitespace=True,header=None,names=[f'Q{i}' for i in range(4)]) # ancestry fractions
admix_P = pd.read_csv(spark_wd+'merged_v2.4.P',delim_whitespace=True,header=None,names=[f'Q{i}' for i in range(4)]) # allele frequencies
'''
Pop0 = African (afr)
Pop1 = Native American (amr)
Pop2 = Asian (asn)
Pop3 = European (eur)
'''
#fam = pd.read_csv(spark_wd+'merged_v2.fam',delim_whitespace=True,header=None,names=['FID','IID','PAT','MAT','SEX','PHEN'])
#fam = fam.join(admix_Q)
#print(len(fam))

#fam[fam['anc_true']=='amr'][[f'Q{i}' for i in range(4)]]



##Combine PCA results with ADMIXTURE results
merged = df.merge(fam,on=['IID']) #inner join, therefore limited to the sample size of df, 7173 samples remaining, 37 removed from fam

merged['anc_est'] = 'unknown'
min_Q = 0.9
merged.loc[((merged.Q0 > min_Q)&
            (merged.Q1 < min_Q)&
            (merged.Q2 < min_Q)&
            (merged.Q3 < min_Q)), 'anc_est'] = 'afr' #ancestry estimate
merged.loc[((merged.Q1 > min_Q)&
            (merged.Q0 < min_Q)&
            (merged.Q2 < min_Q)&
            (merged.Q3 < min_Q)), 'anc_est'] = 'amr'
merged.loc[((merged.Q2 > min_Q)&
            (merged.Q0 < min_Q)&
            (merged.Q1 < min_Q)&
            (merged.Q3 < min_Q)), 'anc_est'] = 'asn'
merged.loc[((merged.Q3 > min_Q)&
            (merged.Q0 < min_Q)&
            (merged.Q1 < min_Q)&
            (merged.Q2 < min_Q)), 'anc_est'] = 'eur'

print(f'For min Q of {min_Q}')
print('\n'.join([f'Individuals of {x} ancestry: {len(merged[merged.anc_est==x])} ({round(len(merged[merged.anc_est==x])/len(merged)*100,2)}%)' for x in ancestry_ls]))
print(f'Total: {len(merged)}')

ref = merged[merged.anc_true!='unknown']
print('\n'.join([f'1KG ref individuals of {x} ancestry: {len(ref[ref.anc_est==x])} ({round(len(ref[ref.anc_est==x])/len(ref)*100,2)}%)' for x in ancestry_ls]))
print(f'Total: {len(ref)}')


spark = merged[merged.anc_true=='unknown']
print('\n'.join([f'SPARK individuals of {x} ancestry: {len(spark[spark.anc_est==x])} ({round(len(spark[spark.anc_est==x])/len(spark)*100,2)}%)' for x in ancestry_ls]))
print(f'Total SPARK: {len(spark)}')


ancestry_ls = ['unknown','eur','amr','afr','asn']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


pcs = [1,2]
for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
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
    minPCx = merged[f'C{pcs[0]}'].min()
    maxPCx = merged[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = merged[f'C{pcs[1]}'].min()
    maxPCy = merged[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark)} SPARK + {len(merged[merged.anc_true!="unknown"])} 1KG ref = {len(merged)} total, min Q = {min_Q})')
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
#    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.minQ_{min_Q}.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)
#    plt.close()


# plot males and females
for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(6*1.5,4*1.5))
    spark_males = df[(df.anc_true=='unknown')&(df.SEX==1)] #SPARK males
    spark_females = df[(df.anc_true=='unknown')&(df.SEX==2)] #SPARK females
    plt.scatter(spark_females[f'C{pcs[0]}'],spark_females[f'C{pcs[1]}'],alpha=0.5,s=20,marker='o',c=colors[1])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
    plt.scatter(spark_males[f'C{pcs[0]}'],spark_males[f'C{pcs[1]}'],alpha=0.5,s=20,marker='o',c=colors[0])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
#    ref_males = df[(df.anc_true!='unknown')&(df.SEX==1)] #SPARK males
#    plt.scatter(ref_males[f'C{pcs[0]}'],ref_males[f'C{pcs[1]}'],alpha=0.5,s=20,marker='x',c=colors[0])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
#    ref_females = df[(df.anc_true!='unknown')&(df.SEX==2)] #SPARK females
#    plt.scatter(ref_females[f'C{pcs[0]}'],ref_females[f'C{pcs[1]}'],alpha=0.5,s=20,marker='x',c=colors[1])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
    plt.legend(['spark fmales','spark males'])
#    for anc_idx, anc in enumerate(ancestry_ls):
#        df_tmp1 = df[(df.anc_est==anc)] #spark individuals
#        plt.scatter(df_tmp1[f'C{pcs[0]}'],df_tmp1[f'C{pcs[1]}'],alpha=0.5,s=10,marker='o', c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
    minPCx = df[f'C{pcs[0]}'].min()
    maxPCx = df[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = df[f'C{pcs[1]}'].min()
    maxPCy = df[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark)} SPARK + {len(df[df.anc_true!="unknown"])} 1KG ref = {len(df)} total)')
    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.male_female.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)


# plot SPARK and 1kg
for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(6*1.5,4*1.5))
    spark = df[(df.anc_true=='unknown')] #SPARK 
    ref = df[(df.anc_true!='unknown')] #1kg ref
    plt.scatter(spark[f'C{pcs[0]}'],spark[f'C{pcs[1]}'],alpha=0.5,s=20,marker='o',c=colors[0])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
    plt.scatter(ref[f'C{pcs[0]}'],ref[f'C{pcs[1]}'],alpha=0.5,s=20,marker='o',c=colors[1])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
#    ref_males = df[(df.anc_true!='unknown')&(df.SEX==1)] #SPARK males
#    plt.scatter(ref_males[f'C{pcs[0]}'],ref_males[f'C{pcs[1]}'],alpha=0.5,s=20,marker='x',c=colors[0])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
#    ref_females = df[(df.anc_true!='unknown')&(df.SEX==2)] #SPARK females
#    plt.scatter(ref_females[f'C{pcs[0]}'],ref_females[f'C{pcs[1]}'],alpha=0.5,s=20,marker='x',c=colors[1])#, c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
    plt.legend(['spark','1kg'])
#    for anc_idx, anc in enumerate(ancestry_ls):
#        df_tmp1 = df[(df.anc_est==anc)] #spark individuals
#        plt.scatter(df_tmp1[f'C{pcs[0]}'],df_tmp1[f'C{pcs[1]}'],alpha=0.5,s=10,marker='o', c = colors[anc_idx])#,edgecolors='k',linewidths=0.1)
    minPCx = df[f'C{pcs[0]}'].min()
    maxPCx = df[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = df[f'C{pcs[1]}'].min()
    maxPCy = df[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark)} SPARK + {len(df[df.anc_true!="unknown"])} 1KG ref = {len(df)} total)')
    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.spark_vs_1kg.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)


#single out one ancestry
ancestry = 'unknown' #options: eur, afr, amr, asn, unknown
for ancestry in ancestry_ls:
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
    minPCx = merged[f'C{pcs[0]}'].min()
    maxPCx = merged[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = merged[f'C{pcs[1]}'].min()
    maxPCy = merged[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark)} SPARK + {len(merged[merged.anc_true!="unknown"])} 1KG ref = {len(merged)} total)')
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.minQ_{min_Q}.PC{pcs[0]}PC{pcs[1]}.ancestry_{ancestry}.png',dpi=600)
    plt.close()
    



#plot PC loadings
PC=4
for PC in range(1,11):
    fig,ax = plt.subplots(figsize=(6*1.5,4*1.5))
    plt.plot(df[f'C{PC}']**2,'.',ms=2)
    plt.xlabel('variant index')
    plt.ylabel(f'PC{PC}^2')
    plt.title(f'PC loadings for PC{PC}')
    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.pc_loadings.PC{PC}.png',dpi=600)
    
    
    
#plot mendelian errors across variants
me_l = pd.read_csv(spark_wd+'preimp3/SPARK.27K.genotype.20190501.hg19_preimp3.lmendel',delim_whitespace=True)
len(me_l)
plt.plot(me_l.N, '.',alpha=1)
plt.xlabel('variant index')
plt.ylabel('# of mendelian errors')

me_l.N.mean()
plt.hist(me_l.N,100)

#plot mendelian errors across FIDs
me_i = pd.read_csv(spark_wd+'preimp3/SPARK.27K.genotype.20190501.hg19_preimp3.imendel',delim_whitespace=True)
me_i = me_i.sort_values(by='FID')

plt.plot(me_i.sort_values(by='N').N.values,'.')
plt.xlabel('individual index')
plt.ylabel('# of mendelian errors')

mean_N = me_i.groupby('FID')['N'].mean()
count_FID = me_i.groupby('FID')['N'].count()

# get FIDS of families with mean mendelian errors > 300
df = pd.DataFrame(mean_N[mean_N>300].index)
df.to_csv(spark_wd+'preimp3/FID_avg_mendelian_gt_300.tsv',sep='\t',index=False,header=False)


plt.plot(np.sort(mean_N.values),'.')
plt.plot([0,len(mean_N)],[300, 300],'k--')
plt.xlabel('Family index')
plt.ylabel('Mean number of mendelian errors')


family_th = 4 #threshold for family size

plt.scatter(x=range(len(mean_N)), 
            y=me_i.groupby('FID')['N'].mean(),
            s=((count_FID<family_th)&(count_FID>0))*10,
            alpha=0.5)
plt.scatter(x=range(len(mean_N)), 
            y=me_i.groupby('FID')['N'].mean(),
            s=(count_FID>=family_th)*10,
            alpha=0.5)


plt.scatter(x=me_i.index/len(me_i.index)*len(mean_N), 
            y=me_i.N,
            s=5,
            alpha=0.5)

plt.scatter(x=count_FID, y=mean_N,s=5)
plt.xlabel('Number of individuals in same FID')
plt.ylabel('Mean # of mendelian errors in FID')

count_FID[count_FID.index==me_i.groupby('FID')['N'].mean().idxmax()]
mean_N[mean_N.index==me_i.groupby('FID')['N'].mean().idxmax()]



# get tables of reported ancestry (only really need adults for PCA plot)

df_child = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_child.csv',sep=',')
df_child = df_child.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
df_adult = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_adult.csv',sep=',')
df_sibling = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_sibling.csv',sep=',')

#df_adult = df_adult.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
df_all = df_child.append(df_adult)
df_all = df_all.append(df_sibling)
df_all = df_all.rename(columns={'subject_sp_id':'IID','family_id':'FID'})

eur = df_all[(df_all.race_white==1)&(df_all.race_more_than_one_calc!=1)][['IID','FID']]

preimp3_fam = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])

## get white parents ##
#children_w_white_mother = df_child[(df_child.mother_race_white==1)&(df_child.mother_race_more_than_one_calc!=1)][['IID','FID']]
#children_w_white_father = df_child[(df_child.father_race_white==1)&(df_child.father_race_more_than_one_calc!=1)][['IID','FID']]
#white_mother_IIDs = set(preimp3_fam[preimp3_fam.IID.isin(children_w_white_mother.IID)].MAT.values)
#white_father_IIDs = set(preimp3_fam[preimp3_fam.IID.isin(children_w_white_father.IID)].PAT.values)
#white_mother_IIDs.remove('0')
#white_father_IIDs.remove('0')
#white_parent_IIDs = white_mother_IIDs.union(white_father_IIDs)
#white_parents = preimp3_fam[preimp3_fam.IID.isin(white_parent_IIDs)]
#
## get all individuals who identify as only being white
#white_individuals = preimp3_fam[preimp3_fam.IID.isin(white_parent_IIDs)|(preimp3_fam.MAT.isin(white_parent_IIDs)&preimp3_fam.PAT.isin(white_parent_IIDs))]

# get parents labeled by ancestry
anc_ls = ['asian','african_amer','native_amer','native_hawaiian','white','hispanic']
anc_dict = dict(zip(anc_ls,[None]*len(anc_ls)))
assert len(preimp3_fam.IID)==len(set(preimp3_fam.IID)), 'there are individuals with duplicate IIDs, be careful!'


#label individuals by ancestry
for ancestry in anc_ls:
    if ancestry != 'hispanic':
        children_w_ancestry_mother = df_child[(df_child[f'mother_race_{ancestry}']==1)&(df_child.mother_race_more_than_one_calc!=1)&
                                              (df_child.mother_hispanic!=1)][['IID','FID']]
        children_w_ancestry_father = df_child[(df_child[f'father_race_{ancestry}']==1)&(df_child.father_race_more_than_one_calc!=1)&
                                              (df_child.father_hispanic!=1)][['IID','FID']]
    elif ancestry=='hispanic':
        children_w_ancestry_mother = df_child[(df_child[f'mother_{ancestry}']==1)][['IID','FID']]
        children_w_ancestry_father = df_child[(df_child[f'father_{ancestry}']==1)][['IID','FID']]
    ancestry_mother_IIDs = set(preimp3_fam[preimp3_fam.IID.isin(children_w_ancestry_mother.IID)].MAT.values)
    ancestry_father_IIDs = set(preimp3_fam[preimp3_fam.IID.isin(children_w_ancestry_father.IID)].PAT.values)
    if '0' in ancestry_mother_IIDs:
        ancestry_mother_IIDs.remove('0')
    if '0' in ancestry_father_IIDs:
        ancestry_father_IIDs.remove('0')
    ancestry_parent_IIDs = ancestry_mother_IIDs.union(ancestry_father_IIDs)
    anc_dict[ancestry] = ancestry_parent_IIDs
    df.loc[df.IID.isin(ancestry_parent_IIDs),'anc_reported']  = ancestry #annotate df of SPARK parents + 1kg ref (see top section of code)
df.loc[(df.anc_true=='unknown')&(df.anc_reported.isna()),'anc_reported'] = 'unknown'

# plot PCs with reported ancestry
single_ancestry = (False, 'asian')
spark_reported = df[(df.anc_true=='unknown')&(df.anc_reported!='unknown')] #SPARK with reported ancestry
ref = df[(df.anc_true!='unknown')] #1kg ref
for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    for anc_idx, anc in enumerate(['eur','afr','amr','asn']):
        ref_anc = ref[(ref.anc_true==anc)] 
        plt.scatter(ref_anc[f'C{pcs[0]}'],ref_anc[f'C{pcs[1]}'],alpha=0.5,s=50,marker='x', c = colors[anc_idx])#,edgecolors='k',linewidths=0.1) 
    for anc_idx, anc in enumerate(['white','african_amer','native_amer','asian','native_hawaiian','hispanic']):
        spark_anc = spark_reported[(spark_reported.anc_reported==anc)]
        if single_ancestry[0]==True:
            if anc==single_ancestry[1]:
                plt.scatter(spark_anc[f'C{pcs[0]}'],spark_anc[f'C{pcs[1]}'],alpha=0.8,s=20,marker='o', c = colors[anc_idx],edgecolors='k',linewidths=0.5) 
        else:
            plt.scatter(spark_anc[f'C{pcs[0]}'],spark_anc[f'C{pcs[1]}'],alpha=0.8,s=20,marker='o', c = colors[anc_idx],edgecolors='k',linewidths=0.5) 
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(['white','african_amer','native_amer','asian','native_hawaiin','hispanic'])]+
                            [Line2D([0],[0],lw=0,markerfacecolor='grey',marker='o',label='SPARK',markeredgecolor='k',alpha=0.5),
                             Line2D([0],[0],lw=0,markeredgecolor='k',marker='x',label='ref')])
    ax.legend(handles =legend_elements)
    minPCx = df[f'C{pcs[0]}'].min()
    maxPCx = df[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = df[f'C{pcs[1]}'].min()
    maxPCy = df[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark_reported)} SPARK w/ reported ancestry + {len(ref)} 1KG ref = {len(spark_reported)+len(ref)} total)')
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.reported_ancestry.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)

for anc, iids in anc_dict.items():
#    print(f'number of {anc} parents in fam file: {len(iids)}')
    print(f'number of {anc} parents in PC results: {spark_reported[spark_reported.anc_reported==anc].shape[0]}')

    


preimp3_fam[preimp3_fam.IID.isin(anc_dict['white'])|((preimp3_fam.PAT.isin(anc_dict['white']))&(preimp3_fam.MAT.isin(anc_dict['white'])))].to_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp3.eur.fam', sep='\t',index=False)

#merged1 = mastertable_20190501.merge(eur, on=['IID'])

preimp3_eur = preimp3_fam[preimp3_fam.FID.isin(eur.FID)]

#df_all.subject_sp_id.shape

all_fam = pd.read_csv(spark_wd+'preimp3/SPARK.27K.genotype.20190501.liftover.v2.autosome.fam',
                      names=['FID','IID','PAT','MAT','SEX','PHEN'],
                      delim_whitespace=True)


len(set(df_all.IID).intersection(all_fam.IID))

parents_fam = pd.read_csv(spark_wd+'preimp3/SPARK.27K.genotype.20190501.hg19_preimp3.parents.fam',
                          names=['FID','IID','PAT','MAT','SEX','PHEN'],
                          delim_whitespace=True)


parents_merge_adults = parents_fam.merge(df_adult, on='IID')

df_adult = df_adult[['IID']+[x for x in df_adult.columns.values if 'race_' in x]]
df_adult = df_adult[df_adult.race_more_than_one_calc!=1]
reported = df.merge(df_adult, on ='IID') # intersection of individuals with reported ancestry in parents IMUS with PCs

len(set(df[df.anc_true=='unknown'].IID).intersection(df_adult.IID))

reported.loc[reported.race_asian==1,'anc_reported'] = 'asn'
reported.loc[reported.race_african_amer==1,'anc_reported'] = 'afr'
reported.loc[reported.race_native_amer==1,'anc_reported'] = 'amr'
reported.loc[reported.race_white==1,'anc_reported'] = 'eur'
reported = reported[~reported.anc_reported.isna()]


reported = reported.append(df[df.anc_true!='unknown'])

for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(6*1.5,4*1.5))
    for anc_idx, anc in enumerate(ancestry_ls): 
        df_tmp2 = reported[(reported.anc_true==anc)&(reported.anc_true!='unknown')] #reference panel individuals
        plt.scatter(df_tmp2[f'C{pcs[0]}'],df_tmp2[f'C{pcs[1]}'],alpha=0.5,s=20,marker='x', c = colors[anc_idx])
    for anc_idx, anc in enumerate(ancestry_ls):
        df_tmp1 = reported[(reported.anc_reported==anc)&(reported.anc_true=='unknown')] #spark individuals
        plt.scatter(df_tmp1[f'C{pcs[0]}'],df_tmp1[f'C{pcs[1]}'],alpha=1,marker='o', c = colors[anc_idx],s=100,edgecolors='k',linewidths=2)
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]+
                        [Line2D([0],[0],lw=0,markerfacecolor='k',marker='o',label='SPARK',markeredgecolor='w'),
                         Line2D([0],[0],lw=0,markeredgecolor='k',marker='x',label='ref')])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)
    minPCx = reported[f'C{pcs[0]}'].min()
    maxPCx = reported[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = reported[f'C{pcs[1]}'].min()
    maxPCy = reported[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark)} SPARK + {len(reported[reported.anc_true!="unknown"])} 1KG ref = {len(reported)} total, min Q = {min_Q})')
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(spark_wd+f'plots/spark_preimp3parentsIMUS.spark_vs_1kg.reported_ancestry.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)

    

# get set of 119 individuals included in 27k mastertable from 05/01 but not in SPARK.30K.array_genotype.20190423.fam
# to check if there is overlap
mastertable_20190501 = pd.read_csv(spark_wd+'SPARK.27K.mastertable.20190501.fam',
                                   delim_whitespace=True,
                                   header=None,
                                   names=['FID','IID','PAT','MAT','SEX','PHEN'])
fam_20190423 = pd.read_csv(spark_wd+'SPARK.30K.array_genotype.20190423.fam',
                           delim_whitespace=True,
                           header=None,
                           names=['FID','IID','PAT','MAT','SEX','PHEN'])
len(mastertable_20190501 )
len(fam_20190423)

removed = mastertable_20190501[~mastertable_20190501.IID.isin(fam_20190423.IID)]

bad_array_markers = pd.read_csv(spark_wd+'SPARK.27K.mastertable.20190501.with-bad-array-marker.tsv',
                                delim_whitespace=True)

removed_merged_bad_array = removed.merge(bad_array_markers, on=['FID','IID'])

master_merged_bad_array = mastertable_20190501.merge(bad_array_markers, on=['FID','IID','PAT','MAT','SEX','PHEN'])

master_merged_bad_array = master_merged_bad_array.drop('call_rate_below_0.9',axis=1)

master_merged_bad_array.to_csv(spark_wd+'SPARK.27K.mastertable.20190501.bad_array.fam',
                               sep='\t',
                               index=False)


preimp3_fam = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])

preimp3_fam[preimp3_fam.FID.isin(df.FID)]





# change PAT/MAT for individuals identified as having wrong parents

preimp3_fam = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])

wrong_parents = pd.read_csv(spark_wd+'me_gt_300.wrongparents.txt',
                            delim_whitespace=True)

# remove PAT/MAT for selected individuals
preimp3_fam.loc[preimp3_fam.IID.isin(wrong_parents.IID1)|preimp3_fam.IID.isin(wrong_parents.IID2),['PAT','MAT']] = 0

preimp3_fam[['FID','IID','PAT','MAT']].to_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp3.parents_fixed.fam',
                                              sep='\t',
                                              index=None)


# plot mendelian errors across individuals, AFTER fixing fam file for mendelian errors > 300, and filtering for mendelian errors

me_f = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp6.fmendel',delim_whitespace=True)
me_i = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp6.imendel',delim_whitespace=True) #mendelian errors per individual across variants
me_l = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp6.lmendel',delim_whitespace=True) #mendelian errors per variant across individuals

me_l.N.mean()

me_i = me_i.sort_values(by='FID')

plt.plot(me_i.N.values,'.')
plt.xlabel('individual index')
plt.ylabel('# of mendelian errors')

mean_N = me_i.groupby('FID')['N'].mean()
count_FID = me_i.groupby('FID')['N'].count()

# get FIDS of families with mean mendelian errors > 300
df = pd.DataFrame(mean_N[mean_N>300].index)

me_i[me_i.FID.isin(df.FID)]
#df.to_csv(spark_wd+'preimp3/FID_avg_mendelian_gt_300.tsv',sep='\t',index=False,header=False)
mean_N.mean()

plt.plot([0,len(mean_N)],[300, 300],'k--',alpha=0.5)
plt.plot(mean_N.values,'.')
plt.xlabel('Family index')
plt.ylabel('Mean number of mendelian errors')



#compare EUR maf between SPARK and 1kg

maf = pd.read_csv(spark_wd+'SPARK_vs_1kg.eur_maf.frq',delim_whitespace=True) #first set is from SPARK, second is 1kg eur
bim = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.cluster_hg19_preimp3.bim',delim_whitespace=True,header=None)
bim = bim.rename(columns={0:'CHR',1:'SNP',4:'A1_bim',5:'A2_bim'}) 
bim['make_A1_effect_allele'] = np.random.randint(low=0,high=2,size=bim.shape[0])==1 #randomly determines whether effect allele is A1 or A2
bim.loc[bim.make_A1_effect_allele,'effectallele'] = bim['A1_bim']
bim.loc[~bim.make_A1_effect_allele,'effectallele'] = bim['A2_bim']
merge_maf_bim = maf.merge(bim,on='SNP')

plt.plot(maf['MAF']-maf['MAF.1'],'.',ms=2,alpha=0.5)
plt.xlabel('variant index')
plt.ylabel('SPARK EUR maf - 1kg EUR maf')

plt.plot(np.sort(maf['MAF'].values-maf['MAF.1'].values),'.',ms=2,alpha=1)
plt.xlabel('variant rank (sorted)')
plt.ylabel('SPARK EUR maf - 1kg EUR maf')

print(f"mean: {(maf['MAF']-maf['MAF.1']).mean()}")
print(f"std: {(maf['MAF']-maf['MAF.1']).std()}")
plt.hist(maf['MAF']-maf['MAF.1'],100)
plt.xlabel('SPARK EUR maf - 1kg EUR maf')
plt.ylabel('count')
plt.savefig(spark_wd+'maf_diff_hist.png',dpi=600)

plt.plot(maf['MAF.1'],maf['MAF']-maf['MAF.1'],'.',ms=1,alpha=0.5)
plt.xlabel('1kg EUR maf')
plt.ylabel('SPARK EUR maf - 1kg EUR maf')
plt.savefig(spark_wd+'1kg_eur_maf_vs_maf_diff.png',dpi=600)

merge_maf_bim.loc[merge_maf_bim['A1']==merge_maf_bim['effectallele'],'spark_EAF'] = merge_maf_bim['MAF']
merge_maf_bim.loc[merge_maf_bim['A1']!=merge_maf_bim['effectallele'],'spark_EAF'] = 1-merge_maf_bim['MAF']
merge_maf_bim.loc[merge_maf_bim['A1.1']==merge_maf_bim['effectallele'],'1kg_EAF'] = merge_maf_bim['MAF.1']
merge_maf_bim.loc[merge_maf_bim['A1.1']!=merge_maf_bim['effectallele'],'1kg_EAF'] = 1-merge_maf_bim['MAF.1']

plt.plot(merge_maf_bim['1kg_EAF'],merge_maf_bim['spark_EAF'],'.',ms=5,alpha=1)
plt.plot([0,1],[0,1],'k--',alpha=0.5)
plt.plot([0.2,1],[0,0.8],'k:',alpha=0.2)
plt.plot([0,0.8],[0.2,1],'k:',alpha=0.2)
plt.xlabel('1kg EUR EAF')
plt.ylabel('SPARK EUR EAF')
plt.savefig(spark_wd+'1kg_eur_vs_spark_EAF.png',dpi=600)




maf.loc[maf['A1']==maf['A1.1'],'1kgalignedAF'] = maf['MAF.1']
maf.loc[maf['A1']!=maf['A1.1'],'1kgalignedAF'] = 1-maf['MAF.1']
plt.plot(maf['1kgalignedAF'],maf['MAF'],'.',ms=5,alpha=1)
plt.plot([0,0.5],[0,0.5],'k--',alpha=0.5)
plt.plot([0.2,0.7],[0,0.5],'k:',alpha=0.2)
plt.plot([0,0.3],[0.2,0.5],'k:',alpha=0.2)
plt.xlabel(f'1kg EUR A1 freq')
plt.ylabel(f'SPARK EUR A1 freq')
plt.savefig(spark_wd+'1kg_eur_vs_spark_A1F.png',dpi=600)






#compare EUR AF between SPARK and 1kg

af = pd.read_csv(spark_wd+'SPARK_vs_1kg.eur_maf.frqx',sep='\t') #first set is from SPARK, second is 1kg eur
af.columns.values
af['C(non-missing,spark_eur)'] = 2*(af['C(HOM A1)']+af['C(HOM A2)']+af['C(HET)'])
af['C(non-missing,1kg_eur)'] = 2*(af['C(HOM A1).1']+af['C(HOM A2).1']+af['C(HET).1'])
af['C(A1freq,spark_eur)'] = (2*af['C(HOM A1)']+af['C(HET)'])/af['C(non-missing,spark_eur)']
af['C(A1.1freq,1kg_eur)'] = (2*af['C(HOM A1).1']+af['C(HET).1'])/af['C(non-missing,1kg_eur)']
af[(af['A1.1']!=af['A1'])].shape[0]
af.shape[0]
af.loc[(af['A1.1']==af['A1']),'C(A1freq,1kg_eur)'] = af['C(A1.1freq,1kg_eur)']
af.loc[(af['A1.1']!=af['A1']),'C(A1freq,1kg_eur)'] = 1- af['C(A1.1freq,1kg_eur)'] #flip A1 allele frequency for rows where A1 in SPARK eur != A1 in 1kg




plt.plot(af['C(A1freq,1kg_eur)'], af['C(A1freq,spark_eur)'],'.',ms=5,alpha=0.5)
plt.plot([0,0.5],[0,0.5],'k--',alpha=0.5)
#plt.xlim([0,0.5])
#plt.ylim([0,0.5])
plt.xlabel('A1 frequency (1kg eur ref)')
plt.ylabel('A1 frequency (SPARK eur)')




# get population labels for HGDP


preimp7_wd = '/Users/nbaya//Documents/lab/genotype-qc/spark/preimp7_imus/'

version = 'v3' # options: v1, v2, v3. v1 used GRCh38 HGDP w/ GRCh37 SPARK data. v2 used GRCh37 versions of both. v3 is the same as v2 except with 70 SNPs removed, which had European AF differences between HGDP and SPARK greater than 0.2 

if version =='v1':
    pca_hgdp = pd.read_csv(preimp7_wd+'preimp7.founders.imus.hgdp.menv.mds', delim_whitespace=True)
elif version == 'v2' or version=='v3':
    pca_hgdp = pd.read_csv(preimp7_wd+f'preimp7.founders.imus.hgdp_{version}.menv.mds', delim_whitespace=True)
    
# PC flips to make it easier to compare:
if version == 'v2':
    pcs_to_flip = [3]
    for pc in pcs_to_flip:
        pca_hgdp[f'C{pc}'] = -pca_hgdp[f'C{pc}']
if version== 'v3':
    pass
    
hgdp_labels = pd.read_csv(preimp7_wd+'hgdp_wgs.20190516.metadata.txt', delim_whitespace=True)

hgdp_labels

hgdp_labels[['sample','population','region']]

set(hgdp_labels.region) # {'AFRICA','AMERICA', 'CENTRAL_SOUTH_ASIA', 'EAST_ASIA', 'EUROPE', 'MIDDLE_EAST', 'OCEANIA'}


hgdp_merge = pca_hgdp.merge(hgdp_labels,left_on='IID',right_on='sample')

hgdp_regions = sorted(list(set(hgdp_labels.region)))


preimp3_fam = pd.read_csv(spark_wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])


# get reported ancestry labels for spark
def add_spark_reported_ancestry(df, fam):
    '''
    df  : The dataframe to which you wish to append the reported ancestry label
    fam : The fam file containing all individuals 
    '''
    wd = '/Users/nbaya//Documents/lab/genotype-qc/spark/'
    df_child = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_child.csv',sep=',')
    df_child = df_child.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
    df_sibling = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_sibling.csv',sep=',')
    df_sibling = df_sibling.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
    df_adult = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_adult.csv',sep=',')
    df_adult = df_adult.rename(columns={'subject_sp_id':'IID','family_id':'FID'})

    anc_ls = ['asian','african_amer','native_amer','native_hawaiian','white','hispanic','other']
    anc_dict = dict(zip(anc_ls,[None]*len(anc_ls)))
    
    assert len(df.IID)==len(set(df.IID)), 'there are individuals with duplicate IIDs, be careful!' # needs to test if we can use IIDs as unique identifiers of individuals
    
    #label individuals by ancestry
    df.loc[df.IID.isin(fam.IID.values),'spark_anc'] = 'unreported'
    for ancestry in anc_ls:
        if ancestry != 'hispanic':
            children_w_ancestry_mother = df_child[(df_child[f'mother_race_{ancestry}']==1)&(df_child.mother_race_more_than_one_calc!=1)&
                                                  (df_child.mother_hispanic!=1)].IID.values
            children_w_ancestry_father = df_child[(df_child[f'father_race_{ancestry}']==1)&(df_child.father_race_more_than_one_calc!=1)&
                                                  (df_child.father_hispanic!=1)].IID.values        
        elif ancestry=='hispanic':
            children_w_ancestry_mother = df_child[(df_child[f'mother_{ancestry}']==1)].IID.values
            children_w_ancestry_father = df_child[(df_child[f'father_{ancestry}']==1)].IID.values
        ancestry_mother_IIDs = set(fam[fam.IID.isin(children_w_ancestry_mother)].MAT.values)
        ancestry_father_IIDs = set(fam[fam.IID.isin(children_w_ancestry_father)].PAT.values)
        ancestry_child_IIDs = set(children_w_ancestry_mother).intersection(children_w_ancestry_father)
        for IIDs in [ancestry_mother_IIDs, ancestry_father_IIDs, ancestry_child_IIDs]:
             if '0' in IIDs:
                 IIDs.remove('0')
        ancestry_sibling_IIDs = fam[(~fam.IID.isin(ancestry_child_IIDs))&(fam.MAT.isin(ancestry_mother_IIDs))&(fam.PAT.isin(ancestry_father_IIDs))].IID.values        
        ancestry_parent_IIDs = ancestry_mother_IIDs.union(ancestry_father_IIDs)
        ancestry_child_reported = ancestry_parent_IIDs.union(ancestry_child_IIDs).union(ancestry_sibling_IIDs)
#        if len(df.loc[(df.IID.isin(ancestry_child_reported))&(df.spark_anc!='unreported')])>0:
#            print(f'WARNING: individuals also with {ancestry} ancestry:\n{df.loc[(df.IID.isin(ancestry_child_reported ))&(df.spark_anc!="unreported")][["IID","spark_anc"]]}')
        df.loc[(df.IID.isin(ancestry_child_reported))&(df.spark_anc=='unreported'),'spark_anc']  = ancestry #annotate df of SPARK parents + 1kg ref (see top section of code)
                                                  
    for ancestry in anc_ls: #must be after the loop for ancestry from the child reported data
        if ancestry != 'hispanic':
            adults_w_ancestry = df_adult[(df_adult[f'race_{ancestry}']==1)&(df_adult.race_more_than_one_calc!=1)&
                                                  (df_adult.hispanic!=1)]['IID'].values
        elif ancestry=='hispanic':
            adults_w_ancestry = df_adult[(df_adult.hispanic==1)].IID.values
#        if len(df.loc[(df.IID.isin(adults_w_ancestry))&(df.spark_anc!='unreported')])>0:
#            print(f'WARNING: individuals also with {ancestry} ancestry:\n{df.loc[(df.IID.isin(adults_w_ancestry))&(df.spark_anc!="unreported")][["IID","spark_anc"]]}')
        df.loc[(df.IID.isin(adults_w_ancestry))&(df.spark_anc=='unreported'),'spark_anc']  = ancestry #annotate df of SPARK parents + 1kg ref (see top section of code)
        anc_dict[ancestry] = (adults_w_ancestry, sum(df.spark_anc==ancestry))
    return df, anc_dict

def add_hgdp_ancestry(df):
    wd = '/Users/nbaya//Documents/lab/genotype-qc/spark/preimp7_imus/'
    hgdp_labels = pd.read_csv(wd+'hgdp_wgs.20190516.metadata.txt', delim_whitespace=True)
#    hgdp_anc = set(hgdp_labels.region.str.lower()+'.'+hgdp_labels.population.str.lower())
    hgdp_regions = sorted(list(set(hgdp_labels.region)))
    df.loc[df.IID.isin(hgdp_labels['sample']),'hgdp_anc'] = 'unreported'
    for region in hgdp_regions:
        ids_from_region = hgdp_labels[hgdp_labels.region==region]['sample'].values    
        df.loc[df.IID.isin(ids_from_region),'hgdp_anc'] = region.lower()
    return df

pca_hgdp = add_spark_reported_ancestry(pca_hgdp,preimp3_fam)
pca_hgdp = add_hgdp_ancestry(pca_hgdp)

#spark_anc_ls = sorted(set(pca[~pca.spark_anc.isna()].spark_anc.values))
spark_anc_ls = ['unreported','white','hispanic','african_amer','asian','other','native_amer','native_hawaiian']

for anc in spark_anc_ls :
    print(f'   spark {anc}: {len(pca_hgdp[pca_hgdp.spark_anc==anc])}')
print(f'** spark total: {len(pca_hgdp[~pca_hgdp.spark_anc.isna()])}')
spark_n_reported = len(pca_hgdp[(~pca_hgdp.spark_anc.isna())&(pca_hgdp.spark_anc!="unreported")])
print(f'** spark total reported: {spark_n_reported}\n')

hgdp_anc_ls = sorted(set(pca_hgdp[~pca_hgdp.hgdp_anc.isna()].hgdp_anc.values))

for region in hgdp_anc_ls:
    print(f'   hgdp {region}: {len(pca_hgdp[pca_hgdp.hgdp_anc==region])}')
print(f'** hgdp_{version} total: {len(pca_hgdp[~pca_hgdp.hgdp_anc.isna()])}')
hgdp_n_reported = len(pca_hgdp[(~pca_hgdp.hgdp_anc.isna())&(pca_hgdp.hgdp_anc!="unreported")])
print(f'** hgdp_{version} total reported: {hgdp_n_reported}')

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# color by SPARK self-reported ancestry
reported_spark_anc_ls = [x for x in spark_anc_ls if x is not 'unreported']
for pcs in [[x,y] for x in range(1,9) for y in range(1,9) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    hgdp_tmp = pca_hgdp[~pca_hgdp.hgdp_anc.isna()]
    plt.plot(hgdp_tmp[f'C{pcs[0]}'],hgdp_tmp[f'C{pcs[1]}'],'x',c='k',ms=5,alpha=0.2)
    for anc_idx, anc in enumerate(reported_spark_anc_ls):
        spark_tmp = pca_hgdp[pca_hgdp.spark_anc==anc]
        plt.plot(spark_tmp[f'C{pcs[0]}'],spark_tmp[f'C{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(reported_spark_anc_ls)]+
                        [Line2D([0],[0],lw=0,markerfacecolor='k',marker='o',label='SPARK',markeredgecolor='w'),
                         Line2D([0],[0],lw=0,markeredgecolor='k',marker='x',label='ref')])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)    
    minPCx = pca_hgdp[f'C{pcs[0]}'].min()
    maxPCx = pca_hgdp[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_hgdp[f'C{pcs[1]}'].min()
    maxPCy = pca_hgdp[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({spark_n_reported} SPARK + '
    title_str += f'{hgdp_n_reported} HGDP, {spark_n_reported+hgdp_n_reported} total)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark_preimp7.founders.imus.hgdp{version}.pc{pcs[0]}_pc{pcs[1]}.png',dpi=600)

# Color by SPARK vs. HGDP
for pcs in [[x,y] for x in range(1,9) for y in range(1,9) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    spark_tmp = pca_hgdp[~pca_hgdp.spark_anc.isna()]
    plt.plot(spark_tmp[f'C{pcs[0]}'],spark_tmp[f'C{pcs[1]}'],'o',alpha=0.5,markeredgecolor='None')
    hgdp_tmp = pca_hgdp[~pca_hgdp.hgdp_anc.isna()]
    plt.plot(hgdp_tmp[f'C{pcs[0]}'],hgdp_tmp[f'C{pcs[1]}'],'o',alpha=0.5,markeredgecolor='None')
    ax.legend(['SPARK','HGDP'])    
    minPCx = pca_hgdp[f'C{pcs[0]}'].min()
    maxPCx = pca_hgdp[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_hgdp[f'C{pcs[1]}'].min()
    maxPCy = pca_hgdp[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark_tmp)} SPARK + '
    title_str += f'{len(hgdp_tmp)} HGDP, {len(spark_tmp)+len(hgdp_tmp)} total)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark_preimp7.founders.imus.spark_vs_hgdp{version}.pc{pcs[0]}_pc{pcs[1]}.png',dpi=600)
    

# Color by HGDP ancestry
for pcs in [[x,y] for x in range(1,9) for y in range(1,9) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    spark_tmp = pca_hgdp[~pca_hgdp.spark_anc.isna()]
    plt.plot(spark_tmp[f'C{pcs[0]}'],spark_tmp[f'C{pcs[1]}'],'x',alpha=0.1,markeredgecolor='k')
    for anc_idx, anc in enumerate(hgdp_anc_ls):
        hgdp_tmp = pca_hgdp[pca_hgdp.hgdp_anc==anc]
        plt.plot(hgdp_tmp[f'C{pcs[0]}'],hgdp_tmp[f'C{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(hgdp_anc_ls)]+
                        [Line2D([0],[0],lw=0,markerfacecolor='k',marker='x',label='SPARK',markeredgecolor='k'),
                         Line2D([0],[0],lw=0,color='k',marker='o',label='ref')])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)
    minPCx = pca_hgdp[f'C{pcs[0]}'].min()
    maxPCx = pca_hgdp[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_hgdp[f'C{pcs[1]}'].min()
    maxPCy = pca_hgdp[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark_tmp)} SPARK + '
    title_str += f'{hgdp_n_reported} HGDP, {len(spark_tmp)+hgdp_n_reported} total)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark_preimp7.founders.imus.hgdp{version}_anc.pc{pcs[0]}_pc{pcs[1]}.png',dpi=600)
    

    
    
# try to generalize ancestry in order to plot reference and target ancestries together
generalized_anc_ls = ['africa','america','asia','europe','oceania','other']
#generalized_anc_dict = {'central_south_asia':'asia','east_asia':'asia','middle_east':'other',
#                        'white':'europe','hispanic':'other','native_amer':'america',
#                        'native_hawaiian':'oceania'}
alt_anc_dict = {'africa':['african_amer'],
                'america':['native_amer'],
                'asia':['central_south_asia','east_asia','asian'],
                'europe':['white'],
                'oceania':['native_hawaiian'],
                'other':['hispanic','middle_east']}

# filter to single ancestry
spark_anc = 'white'
hgdp_region = 'middle_east'

anc_df = pca_hgdp[pca_hgdp.FID.str.lower().str.contains(hgdp_region.lower().replace('_',''))|(pca_hgdp.spark_anc==spark_anc)] #african ancestry
anc_subpops = sorted(set(hgdp_labels[hgdp_labels.region==hgdp_region.upper()].population))

for pcs in [[x,y] for x in range(1,6) for y in range(1,6) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    anc_df_spark_tmp = anc_df[~anc_df.spark_anc.isna()]
    plt.plot(anc_df_spark_tmp[f'C{pcs[0]}'],anc_df_spark_tmp[f'C{pcs[1]}'],'x',alpha=0.1,markeredgecolor='k')
    for subpop_idx, subpop in enumerate(anc_subpops):
        subpop_tmp = anc_df[anc_df.FID.str.contains(subpop)]
        plt.plot(subpop_tmp[f'C{pcs[0]}'],subpop_tmp[f'C{pcs[1]}'], 
                 marker='o' if subpop_idx<10 else '*',
                 linestyle='',
                 c=colors[subpop_idx%10],
                 alpha=0.5,
                 markeredgecolor='None',
                 ms=5 if subpop_idx<10 else 10)
    legend_elements = ([Patch(facecolor=colors[subpop_idx%10],label=subpop) for subpop_idx, subpop in enumerate(anc_subpops)]+
                        [Line2D([0],[0],lw=0,markerfacecolor='k',marker='x',label='SPARK',markeredgecolor='k'),
                         Line2D([0],[0],lw=0,color='k',marker='o',label='ref')])
    ax.legend(handles =legend_elements,
              prop={'size':8})
#    minPCx = pca_hgdp[f'C{pcs[0]}'].min()
#    maxPCx = pca_hgdp[f'C{pcs[0]}'].max()
#    rangePCx = maxPCx-minPCx
#    minPCy = pca_hgdp[f'C{pcs[1]}'].min()
#    maxPCy = pca_hgdp[f'C{pcs[1]}'].max()
#    rangePCy = maxPCy-minPCy
#    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
#    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({spark_anc.capitalize()} ancestry subset: {len(anc_df_spark_tmp)} SPARK + '
    title_str += f'{len(anc_df[anc_df.spark_anc.isna()])} HGDP, {len(anc_df)} total)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark_preimp7.founders.imus.hgdp{version}.{spark_anc}.{hgdp_region}.pc{pcs[0]}_pc{pcs[1]}.png',dpi=600)




# get white individuals in SPARK dataset
pca_hgdp[pca_hgdp.spark_anc=='white'][['FID','IID']].to_csv(
        preimp7_wd+f'spark_preimp7.founders.imus.postpca.white.txt',
        sep='\t',
        header=None,
        index=False)

# get European samples in HGDP
pca_hgdp[pca_hgdp.hgdp_anc=='europe'][['FID','IID']].to_csv(
        preimp7_wd+f'hgdp.postpca.european.txt',
        sep='\t',
        header=None,
        index=False)

# check allele frequencies for white/european samples
spark_af = pd.read_csv(preimp7_wd+'spark_preimp7.founders.imus.postpca.white.frq',
                       delim_whitespace=True)

hgdp_af = pd.read_csv(preimp7_wd+'spark_preimp7.founders.imus.postpca.european.frq',
                      delim_whitespace=True)

merge_af = hgdp_af.merge(spark_af,on=['SNP','CHR'],suffixes=('_hgdp','_spark'))
# align SPARK's alleles to hgdp minor allele ("A1")
merge_af.loc[merge_af.A1_hgdp==merge_af.A1_spark,'hgdpA1_sparkAF']=merge_af['MAF_spark']
merge_af.loc[merge_af.A1_hgdp!=merge_af.A1_spark,'hgdpA1_sparkAF']=1-merge_af['MAF_spark']

plt.plot(merge_af.MAF_hgdp,merge_af.hgdpA1_sparkAF,'.',ms=5,alpha=0.5)

# randomly choose A1:
flip_hgdp_A1 = np.random.randint(low=0,high=2,size=len(merge_af))==1
merge_af.loc[~flip_hgdp_A1,'A1_random'] = merge_af.loc[~flip_hgdp_A1,'A1_hgdp']
merge_af.loc[flip_hgdp_A1,'A1_random'] = merge_af.loc[flip_hgdp_A1,'A2_hgdp']

merge_af.loc[merge_af.A1_hgdp==merge_af.A1_random,'A1F_random_hgdp'] = merge_af.loc[merge_af.A1_hgdp==merge_af.A1_random,'MAF_hgdp']
merge_af.loc[merge_af.A1_hgdp!=merge_af.A1_random,'A1F_random_hgdp'] = 1-merge_af.loc[merge_af.A1_hgdp!=merge_af.A1_random,'MAF_hgdp']

merge_af.loc[merge_af.A1_spark==merge_af.A1_random,'A1F_random_spark'] = merge_af.loc[merge_af.A1_spark==merge_af.A1_random,'MAF_spark']
merge_af.loc[merge_af.A1_spark!=merge_af.A1_random,'A1F_random_spark'] = 1-merge_af.loc[merge_af.A1_spark!=merge_af.A1_random,'MAF_spark']

plt.plot(merge_af.A1F_random_hgdp,merge_af.A1F_random_spark,'.',ms=5,alpha=0.5)
#merge_af_tmp = merge_af[abs(merge_af.A1F_random_hgdp-merge_af.A1F_random_spark)>0.2]
#plt.plot(merge_af_tmp.A1F_random_hgdp,merge_af_tmp.A1F_random_spark,'r.',ms=5,alpha=1)
plt.xlabel('HGDP European AF')
plt.ylabel('SPARK white AF')
plt.plot([0,1],[0,1],'k--',alpha=0.5)
plt.plot([0.2,1],[0,0.8],'k:',alpha=0.5)
plt.plot([0,0.8],[0.2,1],'k:',alpha=0.5)
plt.title('Allele frequencies for white/European samples')
plt.savefig(preimp7_wd+'spark_white.hgdp_euoropean.af.png',dpi=300)

# get outlier SNPs
outliers = merge_af[abs(merge_af.A1F_random_hgdp-merge_af.A1F_random_spark)>0.2]['SNP']
outliers.to_csv(preimp7_wd+'spark_hgdp.AFdiff.snps',sep='\t',header=None,index=False)

# check if outlier SNPs are loading on certain PCs
for c in range(1,11):
    qa = pd.read_csv(preimp7_wd+f'preimp7.founders.imus.hgdp_v2.menv.assomds.C{c}.qassoc',delim_whitespace=True)
    field = 'P'
    nonoutlier_metric = qa.loc[~qa.SNP.isin(outliers),field].abs().mean()
    outlier_metric = qa.loc[qa.SNP.isin(outliers),field].abs().mean()
    print(f'PC{c}\nnon-outlier {field} mean: {round(nonoutlier_metric,4)}\noutlier {field} mean: {round(outlier_metric,4)}\n')
    fig,ax=plt.subplots(figsize=(6,4))
    ax.plot(-np.log10(qa.loc[~qa.SNP.isin(outliers)].P),'.',ms=2)
    ax.plot(-np.log10(qa.loc[qa.SNP.isin(outliers)].P),'r.',ms=5)
    ax.plot([0,len(qa.index)],[-np.log10(5e-8)]*2,'k--',alpha=0.9)
    plt.title(f'PC {c}')
    plt.ylabel('-log10(P)')
    plt.savefig(preimp7_wd+f'pc{c}_gwas.png',dpi=300)
    
# get outlier individuals
def add_spark_reported_ancestries(df, fam):
    '''
    NOTE: This is different from add_spark_reported_ancestries() because it adds
          ALL ancestries reported by SPARK, allowing for a single person to have
          multple ancestries.
    df  : The dataframe to which you wish to append the reported ancestry labels
    fam : The fam file containing all individuals 
    '''
    
    df_child = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_child.csv',sep=',')
    df_child = df_child.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
    df_adult = pd.read_csv(spark_wd+'SPARK_Collection_Version2/bghx_adult.csv',sep=',')
    df_adult = df_adult.rename(columns={'subject_sp_id':'IID','family_id':'FID'})

    anc_ls = ['asian','african_amer','native_amer','native_hawaiian','white','hispanic','other']
    anc_dict = dict(zip(anc_ls,[None]*len(anc_ls)))
    
    assert len(df.IID)==len(set(df.IID)), 'there are individuals with duplicate IIDs, be careful!' # needs to test if we can use IIDs as unique identifiers of individuals
    
    #label individuals by ancestry
    for ancestry in anc_ls:
#        df[ancestry] = False
        if ancestry != 'hispanic':
            children_w_ancestry_mother = df_child[(df_child[f'mother_race_{ancestry}']==1)&(df_child.mother_hispanic!=1)][['IID','FID']]
            children_w_ancestry_father = df_child[(df_child[f'father_race_{ancestry}']==1)&(df_child.father_hispanic!=1)][['IID','FID']]
            adults_w_ancestry = df_adult[(df_adult[f'race_{ancestry}']==1)&(df_adult.hispanic!=1)]['IID'].values
        elif ancestry=='hispanic':
            children_w_ancestry_mother = df_child[(df_child[f'mother_{ancestry}']==1)][['IID','FID']]
            children_w_ancestry_father = df_child[(df_child[f'father_{ancestry}']==1)][['IID','FID']]
            adults_w_ancestry = df_adult[(df_adult.hispanic==1)]['IID'].values
        ancestry_mother_IIDs = set(fam[fam.IID.isin(children_w_ancestry_mother.IID)].MAT.values)
        ancestry_father_IIDs = set(fam[fam.IID.isin(children_w_ancestry_father.IID)].PAT.values)
        if '0' in ancestry_mother_IIDs:
            ancestry_mother_IIDs.remove('0')
        if '0' in ancestry_father_IIDs:
            ancestry_father_IIDs.remove('0')
        ancestry_parent_IIDs = ancestry_mother_IIDs.union(ancestry_father_IIDs)
        ancestry_parents_adults_IIDs = ancestry_parent_IIDs.union(adults_w_ancestry)
        anc_dict[ancestry] = ancestry_parents_adults_IIDs
        df.loc[df.IID.isin(ancestry_parents_adults_IIDs), ancestry]  = True #annotate df of SPARK parents + 1kg ref (see top section of code)
    return df

spark_tmp0 = pca_hgdp[~pca_hgdp.spark_anc.isna()].copy()
hgdp_tmp = pca_hgdp[~pca_hgdp.hgdp_anc.isna()].copy()

spark_tmp = add_spark_reported_ancestries(spark_tmp0, preimp3_fam)

anc_ls = ['asian','african_amer','native_amer','native_hawaiian','white','hispanic','other']

spark_tmp0['n_ancestries'] = spark_tmp0[anc_ls].sum(axis=1)


pc = 8

min_hgdp = hgdp_tmp[f'C{pc}'].min()
max_hgdp = hgdp_tmp[f'C{pc}'].max()
print(f'\nNumber of outliers with PC{pc} < {min_hgdp}: {len(spark_tmp0[(spark_tmp0[f"C{pc}"]<min_hgdp)])}')
print(f'Number of outliers with PC{pc} > {max_hgdp}: {len(spark_tmp0[(spark_tmp0[f"C{pc}"]>max_hgdp)])}')
print('NOTE: Includes all SPARK individuals in founder IMUS, not just those with reported ancestry')

spark_tmp = spark_tmp0[spark_tmp['n_ancestries']>0]
spark_outliers = spark_tmp[(spark_tmp[f'C{pc}']<min_hgdp)|(spark_tmp[f'C{pc}']>max_hgdp)].copy()
print(f'Number of outliers with PC{pc} < {min_hgdp}: {len(spark_tmp[(spark_tmp[f"C{pc}"]<min_hgdp)])}')
print(f'Number of outliers with PC{pc} > {max_hgdp}: {len(spark_tmp[(spark_tmp[f"C{pc}"]>max_hgdp)])}')
print(f'Total number of outliers with PC{pc} < {min_hgdp} or PC{pc} > {max_hgdp}: {len(spark_outliers)}')
print('NOTE: Only using individuals with reported ancestry')


for ancestry in anc_ls:
    spark_outliers[ancestry+'_frac'] = spark_outliers[ancestry]/spark_outliers['n_ancestries']

spark_outliers[['IID']+[ancestry+'_frac' for ancestry in anc_ls]]



# use PC range of certain HGDP ancestry
hgdp_anc= 'east_asia'
if pc==8:
    min_hgdp_anc = hgdp_tmp[hgdp_tmp.hgdp_anc==hgdp_anc][f'C{pc}'].min()
    max_hgdp_anc = hgdp_tmp[hgdp_tmp.hgdp_anc==hgdp_anc][f'C{pc}'].max()    
    spark_outliers = spark_tmp[(spark_tmp[f'C{pc}']<min_hgdp_anc)|(spark_tmp[f'C{pc}']>max_hgdp_anc)].copy()
    print(f'Number of outliers with PC{pc} < {min_hgdp}: {len(spark_tmp[(spark_tmp[f"C{pc}"]<min_hgdp)])}')
    print(f'Number of outliers with PC{pc} > {max_hgdp}: {len(spark_tmp[(spark_tmp[f"C{pc}"]>max_hgdp)])}')
    print(f'Total number of outliers with PC{pc} < {min_hgdp} or PC{pc} > {max_hgdp}: {len(spark_outliers)}')
    print(f'NOTE: Only using individuals with reported ancestry.\n      Using PC range defined by HGDP ancestry "{hgdp_anc}"')
        






















##  Create plots for 1kg

pca_1kg = pd.read_csv(preimp7_wd+'preimp7.founders.imus.1kg.menv.mds', delim_whitespace=True)


def add_1kg_ancestry(df):
    df['1kg_anc'] = df['FID'].str.split('_',expand=True)[3]
    df.loc[df['1kg_anc']=='mix','1kg_anc'] = float('NaN')
    return df

pca_1kg = add_spark_reported_ancestry(pca_1kg,preimp3_fam)
pca_1kg = add_1kg_ancestry(pca_1kg)

spark_anc_ls = ['unreported','white','hispanic','african_amer','asian','other','native_amer','native_hawaiian']

for anc in spark_anc_ls :
    print(f'   spark {anc}: {len(pca_1kg[pca_1kg.spark_anc==anc])}')
print(f'** spark total: {len(pca_1kg[~pca_1kg.spark_anc.isna()])}')
spark_n_reported = len(pca_1kg[(~pca_1kg.spark_anc.isna())&(pca_1kg.spark_anc!="unreported")])
print(f'** spark total reported: {spark_n_reported}\n')

kg_anc_ls = sorted(set(pca_1kg.loc[~pca_1kg['1kg_anc'].isna(),'1kg_anc'].values))

for region in kg_anc_ls:
    print(f'   1kg {region}: {len(pca_1kg[pca_1kg["1kg_anc"]==region])}')
print(f'** 1kg total: {len(pca_1kg[~pca_1kg["1kg_anc"].isna()])}')
kg_n_reported = len(pca_1kg[(~pca_1kg["1kg_anc"].isna())&(pca_1kg["1kg_anc"]!="unreported")])
print(f'** 1kg total reported: {kg_n_reported}')


reported_spark_anc_ls = [x for x in spark_anc_ls if x is not 'unreported']
for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    tmp_1kg= pca_1kg[~pca_1kg["1kg_anc"].isna()]
    plt.plot(tmp_1kg[f'C{pcs[0]}'],tmp_1kg[f'C{pcs[1]}'],'x',c='k',ms=5,alpha=0.2)
    for anc_idx, anc in enumerate(reported_spark_anc_ls):
        spark_tmp = pca_1kg[pca_1kg.spark_anc==anc]
        plt.plot(spark_tmp[f'C{pcs[0]}'],spark_tmp[f'C{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(reported_spark_anc_ls)]+
                        [Line2D([0],[0],lw=0,markerfacecolor='k',marker='o',label='SPARK',markeredgecolor='w'),
                         Line2D([0],[0],lw=0,markeredgecolor='k',marker='x',label='ref')])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)    
    minPCx = pca_1kg[f'C{pcs[0]}'].min()
    maxPCx = pca_1kg[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_1kg[f'C{pcs[1]}'].min()
    maxPCy = pca_1kg[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({spark_n_reported} SPARK + '
    title_str += f'{kg_n_reported} 1KG, {spark_n_reported+kg_n_reported} total)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark_preimp7.founders.imus.1kg.pc{pcs[0]}_pc{pcs[1]}.png',dpi=600)
    
# color by spark vs. 1kg
for pcs in [[x,y] for x in range(1,7) for y in range(1,7) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    spark_tmp = pca_1kg[~pca_1kg.spark_anc.isna()]
    plt.plot(spark_tmp[f'C{pcs[0]}'],spark_tmp[f'C{pcs[1]}'],'o',alpha=0.5,markeredgecolor='None')
    tmp_1kg= pca_1kg[~pca_1kg["1kg_anc"].isna()]
    plt.plot(tmp_1kg[f'C{pcs[0]}'],tmp_1kg[f'C{pcs[1]}'],'o',alpha=0.5,markeredgecolor='None')
    ax.legend(['SPARK','1KG'])    
    minPCx = pca_1kg[f'C{pcs[0]}'].min()
    maxPCx = pca_1kg[f'C{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_1kg[f'C{pcs[1]}'].min()
    maxPCy = pca_1kg[f'C{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({spark_n_reported} SPARK + '
    title_str += f'{kg_n_reported} 1KG, {spark_n_reported+kg_n_reported} total)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark_preimp7.founders.imus.spark_vs_1kg.pc{pcs[0]}_pc{pcs[1]}.png',dpi=600)



# check HGDP qassoc files, extract all significant hits
sig_snps = []
for c in range(1,21):
    qa = pd.read_csv(preimp7_wd+f'preimp7.founders.imus.hgdp_v2.menv.assomds.C{c}.qassoc',delim_whitespace=True)
    pval_thresh = 5e-8    
    sig_snps += qa[qa.P<pval_thresh].SNP.values.tolist() #extract all significant hits
    
np.savetxt(fname=preimp7_wd+'preimp7.founders.imus.hgdp_v2.menv.assomds.sig_snps.qassoc', X=sig_snps,fmt='%s')





# Check PCA results on all SPARK samples
#pca_spark = pd.read_csv(preimp7_wd+'spark.all.admixture_tmp1.scores.tsv.bgz',sep='\t',compression='gzip')
pca_spark = pd.read_csv(preimp7_wd+'spark.all.admixture_tmp1.scores.v2.tsv.bgz',sep='\t',compression='gzip')
pca_spark = pca_spark.rename(columns={'s':'IID','fam_id':'FID','pat_id':'PAT','mat_id':'MAT'})

pca_spark, pca_spark_anc_dict = add_spark_reported_ancestry(pca_spark, preimp3_fam)

# color by SPARK self-reported ancestry
spark_anc_ls = ['unreported','white','hispanic','african_amer','asian','other','native_amer','native_hawaiian']
reported_spark_anc_ls = [x for x in spark_anc_ls if x is not 'unreported']
for pcs in [[x,y] for x in range(1,10) for y in range(1,10) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    spark_tmp = pca_spark[pca_spark.spark_anc=='unreported']
    ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c='grey',alpha=0.5,markeredgecolor='None')
    for anc_idx, anc in enumerate(reported_spark_anc_ls):
        spark_tmp = pca_spark[pca_spark.spark_anc==anc]
        ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
    legend_elements = ([Patch(facecolor='k',label='unreported',alpha=0.5)]+
                       [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(reported_spark_anc_ls)])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)    
    minPCx = pca_spark[f'pc{pcs[0]}'].min()
    maxPCx = pca_spark[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_spark[f'pc{pcs[1]}'].min()
    maxPCy = pca_spark[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n(SPARK w/ reported ancestry: {len(pca_spark[pca_spark.spark_anc!="unreported"])}, SPARK total: {len(pca_spark)})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.all.projected.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)
    
print('\n'.join([f' {ancestry}: {len(pca_spark[pca_spark.spark_anc==ancestry])}' for ancestry in spark_anc_ls]))
print(f'*total: {len(pca_spark)}*')


# plot only nonfounders with reported ancestry
spark_anc_ls = ['unreported','white','hispanic','african_amer','asian','other','native_amer','native_hawaiian']
reported_spark_anc_ls = [x for x in spark_anc_ls if x is not 'unreported']
nonfounders = pca_spark[(~pca_spark.PAT.isna())|(~pca_spark.MAT.isna())]
pca_spark_tmp = nonfounders
for pcs in [[x,y] for x in range(1,6) for y in range(1,6) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    spark_tmp = pca_spark_tmp[pca_spark_tmp.spark_anc=='unreported']
    ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c='grey',alpha=0.2,markeredgecolor='None')
    for anc_idx, anc in enumerate(reported_spark_anc_ls):
        spark_tmp = pca_spark_tmp[pca_spark_tmp.spark_anc==anc]
        ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
#    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(spark_anc_ls)])
    legend_elements = ([Patch(facecolor='k',label='unknown',alpha=0.5)]+
                       [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(reported_spark_anc_ls)])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)    
    minPCx = pca_spark_tmp[f'pc{pcs[0]}'].min()
    maxPCx = pca_spark_tmp[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_spark_tmp[f'pc{pcs[1]}'].min()
    maxPCy = pca_spark_tmp[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n(SPARK non-founders: {len(pca_spark_tmp)}, SPARK non-founders w/ reported ancestry: {len(pca_spark_tmp[pca_spark_tmp.spark_anc!="unreported"])})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.nonfounders.projected.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)

    

    
# plot only founders IMUS, reprojected
spark_founders = pca_hgdp[~pca_hgdp.spark_anc.isna()][['FID','IID']]
pca_founders_imus  = pca_spark[pca_spark.IID.isin(spark_founders.IID)]
pca_spark_tmp = pca_founders_imus
for pcs in [[x,y] for x in range(1,5) for y in range(1,5) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    for anc_idx, anc in enumerate(reported_spark_anc_ls):
        spark_tmp = pca_spark_tmp[pca_spark_tmp.spark_anc==anc]
        ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(reported_spark_anc_ls)])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)    
    minPCx = pca_spark_tmp[f'pc{pcs[0]}'].min()
    maxPCx = pca_spark_tmp[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_spark_tmp[f'pc{pcs[1]}'].min()
    maxPCy = pca_spark_tmp[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({len(pca_spark_tmp[(pca_spark_tmp.spark_anc!="unreported")&~(pca_spark_tmp.spark_anc.isna())])} SPARK founders w/ reported ancestry)'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.founders.projected.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)


# color by founders used for pca vs. projected non-founders
spark_founders_imus = pca_hgdp[~pca_hgdp.spark_anc.isna()][['FID','IID']]
pca_founders_imus = pca_spark[pca_spark.IID.isin(spark_founders_imus.IID)]
#non_pca_founders_imus = pca_spark[~pca_spark.IID.isin(spark_founders_imus.IID)] #everyone else in SPARK not included in the founders IMUS used for PCA
nonfounders = pca_spark[(~pca_spark.PAT.isna())|(~pca_spark.MAT.isna())]
for pcs in [[x,y] for x in range(1,5) for y in range(1,5) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    ax.plot(pca_founders_imus[f'pc{pcs[0]}'],pca_founders_imus[f'pc{pcs[1]}'],'o',alpha=0.25,markeredgecolor='None')
    ax.plot(nonfounders[f'pc{pcs[0]}'],nonfounders[f'pc{pcs[1]}'],'o',alpha=0.25,markeredgecolor='None')
    ax.legend(['Founders IMUS used for PCA','Non-founders'])
    minPCx = pca_spark_tmp[f'pc{pcs[0]}'].min()
    maxPCx = pca_spark_tmp[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_spark_tmp[f'pc{pcs[1]}'].min()
    maxPCy = pca_spark_tmp[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n(PCA founders IMUS: {len(pca_founders_imus)}, Non-founders: {len(nonfounders)})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.pcafoundersimus_vs_nonfounders.projected.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)




# check ADMIXTURE results from founders + HGDP reference
k=5

admix0 = pd.read_csv(preimp7_wd+f'spark.hgdp.admixture_tmp2.{k}.Q',
                   delim_whitespace=True,
                   names=[f'pop{x}' for x in range(k)])
admix_fam = pd.read_csv(preimp7_wd+f'spark.hgdp.admixture_tmp2.fam',
                   delim_whitespace=True,
                   names=['FID','IID','PAT','MAT','SEX','PHEN'])
admix1 = pd.concat([admix_fam,admix0],axis=1)

#remove duplicate IIDs
admix2 = admix1.loc[~admix1.FID.str.contains('preimp7'),:]
#admix2 = admix1.drop_duplicates(subset='IID',keep='last')

#print(f'number of HGDP individuals: {admix2.FID.str.contains('HGDP').sum()}') # number of HGDP individuals


admix3, anc_dict = add_spark_reported_ancestry(df=admix2, 
                                               fam=preimp3_fam)
admix = add_hgdp_ancestry(df=admix3)

# check which pop corresponds to which population
#admix[admix.hgdp_anc=='africa'][[f'pop{x}' for x in range(k)]]
spark_anc_ls = ['unreported','white','hispanic','african_amer','asian','other','native_amer','native_hawaiian']
for spark_anc in spark_anc_ls:
    print(f'\n** {spark_anc} **\n{admix[admix.spark_anc==spark_anc][[f"pop{x}" for x in range(k)]].mean()}\nn = {admix[admix.spark_anc==spark_anc].shape[0]}')
print('\nSPARK\n==================================\nHGDP') # formatting to separate two print statements

hgdp_anc_ls = ['africa', 'america', 'central_south_asia', 'east_asia', 'europe', 'middle_east', 'oceania']
for hgdp_anc in hgdp_anc_ls:
    print(f'\n** {hgdp_anc} **\n{admix[admix.hgdp_anc==hgdp_anc][[f"pop{x}" for x in range(k)]].mean()}\nn = {admix[admix.hgdp_anc==hgdp_anc].shape[0]}')

spark_pop_assign_dict = dict(zip([f'pop{i}' for i in range(k)],[[]]*k))
hgdp_pop_assign_dict = dict(zip([f'pop{i}' for i in range(k)],[[]]*k))

for spark_anc in spark_anc_ls:
    pop = admix[admix.spark_anc==spark_anc][[f"pop{x}" for x in range(k)]].mean().idxmax()
    spark_pop_assign_dict[pop] = spark_pop_assign_dict[pop]+[spark_anc]
    
for hgdp_anc in hgdp_anc_ls:
    pop = admix[admix.hgdp_anc==hgdp_anc][[f"pop{x}" for x in range(k)]].mean().idxmax()
    hgdp_pop_assign_dict[pop] = hgdp_pop_assign_dict[pop]+[hgdp_anc]

inv_hgdp_pop_assign_dict = {}
for key, val in hgdp_pop_assign_dict.items():
    val=val[0]
    inv_hgdp_pop_assign_dict[val] = inv_hgdp_pop_assign_dict.get(val, [])
    inv_hgdp_pop_assign_dict[val].append(key)

    
print('\n'+'\n'.join([f'{k}:{v}' for k,v in spark_pop_assign_dict.items()]))
print('\nSPARK\n==================================\nHGDP') # formatting to separate two print statements
print('\n'.join([f'{k}:{v}' for k,v in hgdp_pop_assign_dict.items()]))

# make stacked bar graph (old version)
admix_reported = admix[(admix.spark_anc!='unreported')&(~admix.spark_anc.isna())].copy()

##rearrange to have european as pop0
#admix_reported['pop0_tmp'] = admix_reported[inv_hgdp_pop_assign_dict['europe']].copy()
#admix_reported[inv_hgdp_pop_assign_dict['europe'][0]] = admix_reported['pop0'].copy()
#admix_reported['pop0'] = admix_reported['pop0_tmp'].copy()
#admix_reported = admix_reported.drop('pop0_tmp',axis=1)
#hgdp_pop_assign_dict[inv_hgdp_pop_assign_dict['europe'][0]] = hgdp_pop_assign_dict['pop0']
#hgdp_pop_assign_dict['pop0'] = ['europe']
#inv_hgdp_pop_assign_dict = {}
#for key, val in hgdp_pop_assign_dict.items():
#    val=val[0]
#    inv_hgdp_pop_assign_dict[val] = inv_hgdp_pop_assign_dict.get(val, [])
#    inv_hgdp_pop_assign_dict[val].append(key)




## plot pca colored by max ADMIXTURE population proportion

max_anc_prop = admix[[f'pop{pop}' for pop in  range(k)]].values.argmax(axis=1) # admixture results with ancestry guess based on maximum ancestry proportion
admix_w_anc_guess = admix.copy()
hgdp_pop_assign_dict = dict(zip([f'pop{i}' for i in range(k)],[[]]*k))
for hgdp_anc in hgdp_anc_ls:
    pop = admix[admix.hgdp_anc==hgdp_anc][[f"pop{x}" for x in range(k)]].mean().idxmax()
    hgdp_pop_assign_dict[pop] = hgdp_pop_assign_dict[pop]+[hgdp_anc]
inv_hgdp_pop_assign_dict = {}
for key, val in hgdp_pop_assign_dict.items():
    val=val[0]
    inv_hgdp_pop_assign_dict[val] = inv_hgdp_pop_assign_dict.get(val, [])
    inv_hgdp_pop_assign_dict[val].append(key)
admix_w_anc_guess['admixture_anc'] = [hgdp_pop_assign_dict[f'pop{i}'][0] for i in max_anc_prop]

spark_tmp = admix_w_anc_guess.merge(pca_spark,on=['IID'])

alt_hgdp_anc_ls = ['europe','east_asia','middle_east','america','central_south_asia','africa','oceania']

for pcs in [[x,y] for x in range(1,9) for y in range(1,9) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    for anc_idx, anc in enumerate(alt_hgdp_anc_ls):
        spark_anc_tmp = spark_tmp[spark_tmp.admixture_anc==anc]
        plt.plot(spark_anc_tmp[f'pc{pcs[0]}'],spark_anc_tmp[f'pc{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(alt_hgdp_anc_ls)])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)
    minPCx = spark_tmp[f'pc{pcs[0]}'].min()
    maxPCx = spark_tmp[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = spark_tmp[f'pc{pcs[1]}'].min()
    maxPCy = spark_tmp[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({len(spark_tmp)} SPARK, k= {k})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.admixture_anc_guess.k_{k}.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)

print(f'Population count (k={k}):')
for hgdp_anc in alt_hgdp_anc_ls:
    print(f'\r{hgdp_anc}: {admix_w_anc_guess[admix_w_anc_guess.admixture_anc==hgdp_anc].shape[0]}')

admix_reported_sub = admix_w_anc_guess.loc[~admix_w_anc_guess.spark_anc.isna()]
admix_reported_sub = admix_reported_sub.sort_values(by=['admixture_anc'],ascending=True)
admix_reported_sub = admix_reported_sub[~admix_reported_sub.pop0.isna()]

fig,ax=plt.subplots(figsize=(6*2,3*2))
prev=0
for pop in range(k):
#    print(pop)
    bar = ax.bar(admix_reported_sub.index, admix_reported_sub[f'pop{pop}'],
                  bottom=(0*admix_reported_sub.index if pop==0 else prev),
                  width=1)
    prev = admix_reported_sub[f'pop{pop}']+prev
plt.legend([v[0] for k,v in hgdp_pop_assign_dict.items()])
plt.xlim([-0.5,len(admix_reported_sub)-0.5])
plt.ylim([0,1])
plt.title(f'ADMIXTURE results (k={k})')
plt.savefig(preimp7_wd+f'admixture_proportions.k_{k}.v2.png',dpi=300)


# make stacked bar graph (updated version)
k=5

preimp7_wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/preimp7_imus/'

admix0 = pd.read_csv(preimp7_wd+f'spark.hgdp.admixture_tmp2.{k}.Q',
                   delim_whitespace=True,
                   names=[f'pop{x}' for x in range(k)])
admix_fam = pd.read_csv(preimp7_wd+f'spark.hgdp.admixture_tmp2.fam',
                   delim_whitespace=True,
                   names=['FID','IID','PAT','MAT','SEX','PHEN'])
admix1 = pd.concat([admix_fam,admix0],axis=1)
#remove duplicate IIDs
admix2 = admix1.loc[~admix1.FID.str.contains('preimp7'),:] #remove duplicates (individuals in SPARK founders IMUS)
#admix2 = admix1.drop_duplicates(subset='IID',keep='last')

print(f'number of SPARK individuals: {admix2.FID.str.contains("SF").sum()}') # number of SPARK individuals
print(f'number of HGDP individuals: {admix2.FID.str.contains("HGDP").sum()}') # number of HGDP individuals
admix3, anc_dict = add_spark_reported_ancestry(df=admix2, 
                                               fam=preimp3_fam)
anc_sum = 0
for key,val in anc_dict.items():
    print(f'{key}: {val[1]}')
    print(key, sum(admix3.spark_anc==key))
    anc_sum += val[1]
admix = add_hgdp_ancestry(df=admix3)
admix = admix.reset_index()

# check which pop corresponds to which population
#admix[admix.hgdp_anc=='africa'][[f'pop{x}' for x in range(k)]]
spark_anc_ls = ['unreported','white','hispanic','african_amer','asian','other','native_amer','native_hawaiian']
for spark_anc in spark_anc_ls:
    print(f'\n** {spark_anc} **\n{admix[admix.spark_anc==spark_anc][[f"pop{x}" for x in range(k)]].mean()}\nn = {admix[admix.spark_anc==spark_anc].shape[0]}')
print('\nSPARK\n==================================\nHGDP') # formatting to separate two print statements

hgdp_anc_ls = ['africa', 'america', 'central_south_asia', 'east_asia', 'europe', 'middle_east', 'oceania']
for hgdp_anc in hgdp_anc_ls:
    print(f'\n** {hgdp_anc} **\n{admix[admix.hgdp_anc==hgdp_anc][[f"pop{x}" for x in range(k)]].mean()}\nn = {admix[admix.hgdp_anc==hgdp_anc].shape[0]}')

#k_eurpop_dict = {6:'pop0', 7:'', 8:''} #dictionary mapping k to which population is matched with white/European ancestry

np.argmax(admix[[f'pop{i}' for i in range(k)]].values,axis=1) #get maximum proportion
for i in range(k):
    print(f'\n** pop{i} **\n{np.sum(np.argmax(admix[[f"pop{i}" for i in range(k)]].values,axis=1)==i)}') #print number of individuals with maximum proportion at pop{i}
    
pop_ct = [np.sum(np.argmax(admix[[f"pop{i}" for i in range(k)]].values,axis=1)==i) for i in range(k)]
sorted_pop_idx = np.argsort(pop_ct)[::-1].tolist() #population indices by count, in decreasing order


print('... subsetting to SPARK datasets! ...')
admix = admix[admix.hgdp_anc.isna()]



for prop_threshold in [0.85]:
    #assert prop_threshold>=0.5,'Threshold should be greater than or equal to 0.5, otherwise multiple population assignments are possible'
    
    for idx in admix.index:
        admix.loc[admix.index==idx,'admix_pop_guess'] = np.argmax(
                (admix.loc[admix.index==idx][[f"pop{i}" for i in range(k)]]>prop_threshold).values).astype(int) if (
                        admix.loc[admix.index==idx][[f"pop{i}" for i in range(k)]].max(axis=1)>prop_threshold).values[0] else 'unknown'
    
    print('UNSORTED admixture population guess')
    for i in sorted_pop_idx:
        print(f'\n** pop {i} **\n{np.sum(admix.admix_pop_guess==i)}') #print number of individuals in pop{i}
    
    
    assigned_pop_ct = [np.sum(admix.admix_pop_guess==i) for i in sorted_pop_idx]
    sorted_assigned_pop_idx = np.argsort(assigned_pop_ct)[::-1].tolist() #population indices by count, in decreasing order
    sorted_assigned_pop_idx_dict = dict(zip(sorted_pop_idx, sorted_assigned_pop_idx))
    inv_sorted_assigned_pop_idx_dict = dict(zip(sorted_assigned_pop_idx,sorted_pop_idx)) #inverse map
    
    print('SORTED admixture population guess')
    admix['sorted_admix_pop_guess'] = admix.admix_pop_guess.map(lambda x: sorted_assigned_pop_idx_dict[x] if x is not 'unknown' else '999').astype('int')
    for i in list(range(k))+[999]:
        print(f'\n** pop {i} **\n{np.sum(admix.sorted_admix_pop_guess==i)}') #print number of individuals with maximum proportion at pop{i}
    
    admix_tmp = admix.copy() #.sample(n=500, axis=0)
    #admix_tmp = admix_tmp.sort_values(by=f'pop{sorted_pop_idx[0]}',ascending=False) #sort by proportion of the most prevalent ancestral population
    admix_tmp = admix_tmp.sort_values(by=f'sorted_admix_pop_guess',ascending=True) #sort by population size
    admix_tmp = admix_tmp.reset_index()
    
    
    spark_pop_assign_dict = dict(zip([f'pop{i}' for i in range(k)],[[]]*k))
    hgdp_pop_assign_dict = dict(zip([f'pop{i}' for i in range(k)],[[]]*k))
    
    for spark_anc in spark_anc_ls:
        pop = admix[admix.spark_anc==spark_anc][[f"pop{x}" for x in range(k)]].mean().idxmax()
        spark_pop_assign_dict[pop] = spark_pop_assign_dict[pop]+[spark_anc]
        
    for hgdp_anc in hgdp_anc_ls:
        pop = admix[admix.hgdp_anc==hgdp_anc][[f"pop{x}" for x in range(k)]].mean().idxmax()
        hgdp_pop_assign_dict[pop] = hgdp_pop_assign_dict[pop]+[hgdp_anc]
    
    inv_hgdp_pop_assign_dict = {}
    for key, val in hgdp_pop_assign_dict.items():
        val=val[0]
        inv_hgdp_pop_assign_dict[val] = inv_hgdp_pop_assign_dict.get(val, [])
        inv_hgdp_pop_assign_dict[val].append(key)
    
    for pop in [f'pop{pop}' for pop in range(k)]:
        hgdp_pop_assign_dict[pop] = '/'.join(hgdp_pop_assign_dict[pop])
        
    print('\n'+'\n'.join([f'{k}:{v}' for k,v in spark_pop_assign_dict.items()]))
    print('\nSPARK\n==================================\nHGDP') # formatting to separate two print statements
    print('\n'.join([f'{k}:{v}' for k,v in hgdp_pop_assign_dict.items()]))
    
    for pop_idx in range(k):
        pop = f'pop{pop_idx}'
        print(f'{hgdp_pop_assign_dict[pop]} ({pop}): {admix[admix.admix_pop_guess==pop_idx].shape[0]} ')
    
    
    
    # stacked bar plot, sorted by population with highest prevalence

    start = dt.now()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig,ax=plt.subplots(figsize=(8*2,3*2))
    for i in list(range(k))+[999]:
        admix_pop_tmp = admix_tmp[admix_tmp.sorted_admix_pop_guess==i]
        min_idx = admix_pop_tmp.index.min()
        if i != 999:
            original_pop = inv_sorted_assigned_pop_idx_dict[i]
            admix_pop_tmp = admix_pop_tmp.sort_values(by=f'pop{original_pop}',ascending=False) #sort by proportion of population that individual belongs to
        else:
            admix_pop_tmp = admix_pop_tmp.sort_values(by=f'pop{sorted_pop_idx[0]}',ascending=False) #sort unknowns by proportion of the most prevalent ancestral population
        admix_pop_tmp = admix_pop_tmp.reset_index(drop=True)
    #    print(pop)
        prev=0
        for j in sorted_pop_idx:
            print(f'{i}.{j}')
            bar = ax.bar(min_idx+admix_pop_tmp.index, admix_pop_tmp[f'pop{j}'],
                         bottom=(0*admix_pop_tmp.index if j==sorted_pop_idx[0] else prev),
                         width=1, color=colors[j])
            prev = admix_pop_tmp[f'pop{j}']+prev
    #plt.legend([v[0] for k,v in hgdp_pop_assign_dict.items()])
#    plt.legend([f'pop{i}' for i in sorted_pop_idx])
    legend_elements = [Patch(facecolor=colors[inv_sorted_assigned_pop_idx_dict[pop_idx]],label=['Europe/Middle East','East Asia/Oceania','Africa','America','Central South Asia'][pop_idx]) for pop_idx in range(k)]
    ax.legend(handles =legend_elements)
    ax.plot([0,len(admix_tmp)],[prop_threshold]*2,'k--',)
    plt.xlim([-0.5,len(admix_tmp)-0.5])
    plt.ylim([0,1])
    plt.title(f'ADMIXTURE results (k={k})')
    plt.savefig(preimp7_wd+f'admixture_proportions.k_{k}.threshold_{prop_threshold}.v4.png',dpi=600)    
    elapsed = dt.now()-start
    print(f'elapsed: {round(elapsed.total_seconds()/60,2)} min')
    
    # PCA plot colored by assigned ancestry according to ADMIXTURE result
    pca_spark = pd.read_csv(preimp7_wd+'spark.all.admixture_tmp1.scores.v2.tsv.bgz',sep='\t',compression='gzip')
    pca_spark = pca_spark.rename(columns={'s':'IID','fam_id':'FID','pat_id':'PAT','mat_id':'MAT'})
    pca_spark = pca_spark.merge(admix, on='IID')
#    
#    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#    colors.remove(colors[7])
#    for pcs in [[x,y] for x in range(1,11) for y in range(1,11) if (x-y<0 and x-y>=-1)]:
#        fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
#        pca_spark_tmp = pca_spark[pca_spark.sorted_admix_pop_guess==999]
#        plt.plot(pca_spark_tmp[f'pc{pcs[0]}'],pca_spark_tmp[f'pc{pcs[1]}'],'o',c='grey',alpha=0.1,markeredgecolor='None')
#        for pop_idx, pop in enumerate(list(range(k))):
#            pca_spark_tmp = pca_spark[pca_spark.sorted_admix_pop_guess==pop]
#            plt.plot(pca_spark_tmp[f'pc{pcs[0]}'],pca_spark_tmp[f'pc{pcs[1]}'],'o',c=colors[inv_sorted_assigned_pop_idx_dict[pop_idx]],alpha=0.5,markeredgecolor='None')
#        legend_elements = ([Patch(facecolor='grey',label='unknown')]+
#                           [Patch(facecolor=colors[inv_sorted_assigned_pop_idx_dict[pop_idx]],label=hgdp_pop_assign_dict[f'pop{inv_sorted_assigned_pop_idx_dict[pop_idx]}']) for pop_idx,pop in enumerate(list(range(k)))])
#        #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
#        ax.legend(handles =legend_elements)
#        minPCx = pca_spark[f'pc{pcs[0]}'].min()
#        maxPCx = pca_spark[f'pc{pcs[0]}'].max()
#        rangePCx = maxPCx-minPCx
#        minPCy = pca_spark[f'pc{pcs[1]}'].min()
#        maxPCy = pca_spark[f'pc{pcs[1]}'].max()
#        rangePCy = maxPCy-minPCy
#        plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
#        plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
#        title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({len(pca_spark)} SPARK, k= {k}, minimum pop. proportion threshold = {prop_threshold})'
#        plt.title(title_str)
#        plt.xlabel(f'PC{pcs[0]}')
#        plt.ylabel(f'PC{pcs[1]}')
#        plt.savefig(preimp7_wd+f'spark.admixture_anc_guess.k_{k}.threshold_{prop_threshold}.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)
        
    
        
    anc_sum = 0
    print(f'... prop_threshold: {prop_threshold} ...')
    for pop_idx in range(k):
        pop = f'pop{pop_idx}'
        pop_ct = pca_spark[pca_spark.admix_pop_guess==pop_idx].shape[0]
    #    print(f'{hgdp_pop_assign_dict[pop]} ({pop}): {pop_ct} ')
        print(f'{hgdp_pop_assign_dict[pop]}: {pop_ct} ')
        print(f'number of cases: {round((pca_spark[pca_spark.admix_pop_guess==pop_idx].PHEN==2).sum(), 3)}')
        if pop_idx==0:
            print(f'number of self-reported cases: {(pca_spark[(pca_spark.admix_pop_guess==pop_idx)&(pca_spark.spark_anc=="white")].PHEN==2).sum()}')
            print(f'number of dropped self-reported cases: {(pca_spark[(pca_spark.admix_pop_guess!=pop_idx)&(pca_spark.spark_anc=="white")].PHEN==2).sum()}')
        if pop_idx==4:
            print(f'number of self-reported cases: {(pca_spark[(pca_spark.admix_pop_guess==pop_idx)&(pca_spark.spark_anc=="african_amer")].PHEN==2).sum()}')
            print(f'number of dropped self-reported cases: {(pca_spark[(pca_spark.admix_pop_guess!=pop_idx)&(pca_spark.spark_anc=="african_amer")].PHEN==2).sum()}')
        print('\n')
        anc_sum+=pop_ct
    unknown_ct = pca_spark[pca_spark.admix_pop_guess=="unknown"].shape[0]
    print(f'unknown: {unknown_ct}')
    print(f'number of cases: {round((pca_spark[pca_spark.admix_pop_guess=="unknown"].PHEN==2).sum(), 3)}')
    print(f'total: {pca_spark[~pca_spark.admix_pop_guess.isna()].shape[0]} (exp: {anc_sum+unknown_ct})\n')
round((pca_spark.PHEN==2).sum(), 3)
    


# plot PCs of PCA run only on founders IMUS
import ast

founders = pd.read_csv(preimp7_wd+'preimp7.founders.imus.hgdp_v3.menv.scores.founders.tsv',
                       sep='\t')
loadings_ls = [ast.literal_eval(x) for x in founders['scores']]
founders['scores'] = loadings_ls
for i in range(20):
    founders[f'pc{i+1}'] = np.asarray(founders.scores.values.tolist())[:,i]
founders = founders.rename({'s':'IID'},axis=1)
founders, anc_dict = add_spark_reported_ancestry(founders, preimp3_fam)    

pca_spark_tmp = founders
for pcs in [[x,y] for x in range(1,12) for y in range(1,12) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    spark_tmp = pca_spark_tmp[pca_spark_tmp.spark_anc=='unreported']
    ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c='grey',alpha=0.2,markeredgecolor='None')
    for anc_idx, anc in enumerate(reported_spark_anc_ls):
        spark_tmp = pca_spark_tmp[pca_spark_tmp.spark_anc==anc]
        ax.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',c=colors[anc_idx],alpha=0.5,markeredgecolor='None')
#    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(spark_anc_ls)])
    legend_elements = ([Patch(facecolor='k',label='unknown',alpha=0.5)]+
                       [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(reported_spark_anc_ls)])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)    
    minPCx = pca_spark_tmp[f'pc{pcs[0]}'].min()
    maxPCx = pca_spark_tmp[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = pca_spark_tmp[f'pc{pcs[1]}'].min()
    maxPCy = pca_spark_tmp[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n(SPARK founders: {len(pca_spark_tmp)}, SPARK non-founders w/ reported ancestry: {len(pca_spark_tmp[pca_spark_tmp.spark_anc!="unreported"])})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.onlyfounderspca.projected.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)

    
    
    
    
## Plot CV error
## use this command: 
## gsutil cat gs://qc-nbaya/spark/array_May2019/preimputation/preimp7_imus/admixture_hgdp/gcp_results/cv/log*.out | grep -h CV
k_ls = range(1,12)
cv_ls = [0.56031, 0.54731, 0.54045, 0.53801, 0.53717, 0.53664, 0.53631, 0.53622, 
         0.53615, 0.53602, 0.53599]
plt.plot(k_ls, cv_ls, '.-')
plt.xlabel('K (number of populations)')
plt.ylabel('CV error')
plt.xticks(k_ls)
plt.title('K vs. CV error')
plt.savefig(preimp7_wd+f'sparkfoundersimus_w_hgdp.cv_error.png',dpi=300)



## cases in bghx_{child,adult} but not PLINK
bghx_not_plink_cases = set(df_child.IID).union(df_adult.IID).difference(df.IID)
df1 = df_child[df_child.IID.isin(bghx_not_plink_cases)][['FID','IID']]
df1 = df1.append(df_adult[df_adult.IID.isin(bghx_not_plink_cases)][['FID','IID']])
df1.to_csv(preimp7_wd+'cases.bghx_not_plink.tsv',sep='\t',index=False)

## cases in PLINK but not bghx_{child,adult} 
plink_not_bghx_cases = set(df[df.PHEN==2].IID).difference(set(df_child.IID).union(df_adult.IID))
df2 = df[df.IID.isin(plink_not_bghx_cases)][['FID_x','IID']]
df2 = df2.rename({'FID_x':'FID'},axis=1)
df2.to_csv(preimp7_wd+'cases.plink_not_bghx.tsv',sep='\t',index=False)


## cases in PLINK and bghx_{child,adult} 
plink_and_bghx_cases = set(df[df.PHEN==2].IID).intersection(set(df_child.IID).union(df_adult.IID))
df3 = df[df.IID.isin(plink_and_bghx_cases)][['FID_x','IID']]
df3 = df3.rename({'FID_x':'FID'},axis=1)
df3.to_csv(preimp7_wd+'cases.plink_and_bghx.tsv',sep='\t',index=False)


## Plot PCA with ADMIXTURE results, highlighting missing samples and EUR parents
missing = pd.read_csv('/Users/nbaya/Downloads/missing_eur_offspring.txt',
                      delim_whitespace=True)
missing = missing.rename(columns={'offspring':'IID', 'father':'PAT','mother':'MAT'})

preimp7_wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/preimp7_imus/'

pca_spark = pd.read_csv(preimp7_wd+'spark.all.admixture_tmp1.scores.v2.tsv.bgz',sep='\t',compression='gzip')
pca_spark = pca_spark.rename(columns={'s':'IID','fam_id':'FID'})
pca_spark = pca_spark.drop(['pat_id','mat_id'], axis=1)
pca_spark, pca_spark_anc_dict = add_spark_reported_ancestry(pca_spark,preimp3_fam)

k=5

admix0 = pd.read_csv(preimp7_wd+f'spark.hgdp.admixture_tmp2.{k}.Q',
                   delim_whitespace=True,
                   names=[f'pop{x}' for x in range(k)])
admix_fam = pd.read_csv(preimp7_wd+f'spark.hgdp.admixture_tmp2.fam',
                   delim_whitespace=True,
                   names=['FID','IID','PAT','MAT','SEX','PHEN'])
admix1 = pd.concat([admix_fam,admix0],axis=1)
#remove duplicate IIDs
admix2 = admix1.loc[~admix1.FID.str.contains('preimp7'),:] #remove duplicates (individuals in SPARK founders IMUS)

print(f'Number of missing offspring that made it to ADMIXTURE: {admix2.merge(missing, on="IID").shape[0]}')

#admix2 = admix1.drop_duplicates(subset='IID',keep='last')

print(f'number of SPARK individuals in ADMIXTURE results: {admix2.FID.str.contains("SF").sum()}') # number of SPARK individuals
print(f'number of HGDP individuals in ADMIXTURE results: {admix2.FID.str.contains("HGDP").sum()}') # number of HGDP individuals


admix_pca = admix2.merge(pca_spark, on=['FID','IID'])

post_admix_eur_fam = pd.read_csv(preimp7_wd+'spark.eur.20191010.fam',
                                 delim_whitespace=True)
post_admix_eur_fam = post_admix_eur_fam.drop('PHEN', axis=1)
                                 

preimp_final = pd.read_csv('/Users/nbaya/Downloads/spark.eur.preimp_final.fam', delim_whitespace=True,
                           names=['FID','IID','PAT','MAT','SEX','PHEN'])
missing['parent_in_preimp_final_ct'] = missing.PAT.isin(preimp_final.IID).astype(int) + missing.MAT.isin(preimp_final.IID).astype(int)
missing['parent_in_preimp_final_ct'].value_counts()

admix_pca['parent_in_preimp_final_ct'] = admix_pca.PAT.isin(preimp_final.IID).astype(int) + admix_pca.MAT.isin(preimp_final.IID).astype(int)

admix_pca['missing_offspring'] = admix_pca.IID.isin(missing.IID)

admix_pca['eur_parent_of_missing_offspring'] = (admix_pca.pop0>0.85)&(admix_pca.IID.isin(missing.PAT)|admix_pca.IID.isin(missing.MAT))
admix_pca['post_admix_eur_parent_of_missing_offspring'] = (admix_pca.IID.isin(post_admix_eur_fam.IID))&(admix_pca.IID.isin(missing.PAT)|admix_pca.IID.isin(missing.MAT))
admix_pca['imputed_parent_of_missing_offspring'] = (admix_pca.IID.isin(preimp_final.IID))&(admix_pca.IID.isin(missing.PAT)|admix_pca.IID.isin(missing.MAT))

admix_pca['noneur_parent_of_missing_offspring'] = (admix_pca.pop0<=0.85)&(admix_pca.IID.isin(missing.PAT)|admix_pca.IID.isin(missing.MAT))

admix_pca[admix_pca.missing_offspring].parent_in_preimp_final_ct.value_counts()

admix_pca[admix_pca.missing_offspring&(admix_pca.parent_in_preimp_final_ct==2)]

original_fam = pd.read_csv(f'{spark_wd}/SPARK.27K.genotype.20190501.fam',
                           delim_whitespace=True, names=['FID','IID','PAT','MAT','SEX','PHEN'])

missing[(missing.parent_in_preimp_final_ct==1)&(missing.IID.isin(original_fam.IID))&(missing.IID.isin(admix1.IID))].shape

admix_pca[admix_pca['noneur_parent_of_missing_offspring']].IID.isin(missing.PAT).sum()
admix_pca[admix_pca['noneur_parent_of_missing_offspring']].IID.isin(missing.MAT).sum()

missing[(missing.PAT.isin(admix_pca[admix_pca['noneur_parent_of_missing_offspring']].IID))|
        (missing.MAT.isin(admix_pca[admix_pca['noneur_parent_of_missing_offspring']].IID))]

missing[missing.PAT.isin(preimp_final.IID)&(missing.MAT.isin(preimp_final.IID))].shape
missing[missing.PAT.isin(preimp_final.IID)&(missing.MAT.isin(preimp_final.IID))&(missing.IID.isin(admix_pca.IID))].shape

admix_pca['report_eur_parent_of_missing_offspring'] = (admix_pca.spark_anc=='white')&(admix_pca.IID.isin(missing.PAT)|admix_pca.IID.isin(missing.MAT))

admix_pca['report_noneur_parent_of_missing_offspring'] = (admix_pca.pop0!='white')&(admix_pca.IID.isin(missing.PAT)|admix_pca.IID.isin(missing.MAT))

admix_pca['missing_offspring_with_pop0>0.85'] = (admix_pca.pop0>0.85)&(admix_pca.missing_offspring)
admix_pca['missing_offspring_in_post_admix_eur_set'] = (admix_pca.IID.isin(post_admix_eur_fam.IID))&(admix_pca.missing_offspring)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
sample_types = ['other','imputed_parent_of_missing_offspring','missing_offspring']
#sample_types = ['other','imputed_parent_of_missing_offspring','missing_offspring_with_pop0>0.85']
sample_types = ['other','noneur_parent_of_missing_offspring','missing_offspring']
for pcs in [(2*x+1, 2*(x+1)) for x in range(0,5)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    size=5
    for sample_type in sample_types:
        if sample_type=='other':
            spark_tmp = admix_pca[(~admix_pca[sample_types[1]])&(~admix_pca[sample_types[2]])]
            markeredgecolor=None
            color='grey'
            alpha=0.05
        elif sample_type == sample_types[1]:
            spark_tmp = admix_pca[admix_pca[sample_type]]
            color=colors[0]
            alpha=0.2
        elif sample_type is sample_types[2]:
            spark_tmp = admix_pca[admix_pca[sample_type]]
            color= colors[1]
            alpha=0.2
        print(sample_type, ':', spark_tmp.shape[0])
        plt.plot(spark_tmp[f'pc{pcs[0]}'],spark_tmp[f'pc{pcs[1]}'],'o',
                 c=color,alpha=alpha,markeredgecolor=markeredgecolor, ms=size)
    legend_elements = ([Patch(facecolor=c,label=sample_type.replace('_',' ')) for sample_type, c in list(zip(sample_types, ['grey',colors[0], colors[1]]))[::-1]])
#    legend_elements = ([Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(alt_hgdp_anc_ls)])
    #legend_elements = [Patch(facecolor=colors[anc_idx],label=anc) for anc_idx,anc in enumerate(ancestry_ls)]
    ax.legend(handles =legend_elements)
    minPCx = admix_pca[f'pc{pcs[0]}'].min()
    maxPCx = admix_pca[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = admix_pca[f'pc{pcs[1]}'].min()
    maxPCy = admix_pca[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n({len(admix_pca)} SPARK, k= {k})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
#    plt.savefig(f'/Users/nbaya/Downloads/missing_eur_offspring.v1.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)
#    plt.savefig(f'/Users/nbaya/Downloads/missing_eur_offspring.v2.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)
    plt.savefig(f'/Users/nbaya/Downloads/missing_eur_offspring.v3.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)


## check samples removed by ListSPID27270

new = pd.read_csv(preimp7_wd+'ListSPID27270.txt',delim_whitespace=True)

pca_spark['keep'] = pca_spark.IID.isin(new.spid)

master_fam = pd.read_csv('/Users/nbaya/Documents/lab/genotype-qc/spark/SPARK.30K.array_genotype.20190423.fam',
                         delim_whitespace=True,
                         names=['FID','IID','PAT','MAT','SEX','PHEN'])
master_fam['keep'] = master_fam.IID.isin(new.spid)

bad_array= pd.read_csv(spark_wd+'SPARK.27K.mastertable.20190501.with-bad-array-marker.tsv',
                                sep='\t')

not_included = set(new.spid).difference(master_fam.IID) #individuals who were not included (for genotyping?)
plink_not_spark = set(master_fam.IID).difference(new.spid)
master_not_pca = set(master_fam.IID).difference(pca_spark.IID)

not_inc_bad_array = bad_array[bad_array.spid.isin(not_included)]

(not_inc_bad_array['call_rate_below_0.9']=='X').sum()
not_inc_bad_array.array_idat.isna().mean()

(bad_array[bad_array.spid.isin(master_not_pca)]['call_rate_below_0.9']=='X').sum()


len(set(pca_spark.IID).difference(new.spid))

len(not_included)

new[new.spid.isin(not_included)].to_csv(spark_wd+'ListSPID27270_not_in_SPARK.30K.array_genotype.20190423.tsv',sep='\t',index=False)

# remove individuals NOT in ListSPID27270
pca_spark = pca_spark[~pca_spark.IID.isin(plink_not_spark)]

# select EUR individuals
spark_eur = pca_spark[pca_spark.pop0>=0.85]
spark_eur_fam = master_fam[master_fam.IID.isin(spark_eur.IID)]
spark_eur_fam[['FID','IID','PAT', 'MAT', 'SEX', 'PHEN']].to_csv(preimp7_wd+'spark.eur.20191010.fam',sep='\t',index=False)



# plot PCs for EUR samples
eur_pca2_wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/eur_pca2/'
eur_pca = pd.read_csv(eur_pca2_wd+'spark.eur.preimp3.qc_snp.scores.tsv.bgz',
                      compression='gzip',sep='\t')
eur_pca = eur_pca.rename(columns={'s':'IID','fam_id':'FID','pat_id':'PAT',
                                  'mat_id':'MAT'})
eur_pca_founders = pd.read_csv(eur_pca2_wd+'spark.eur.founders.menv.scores.founders.tsv', # post pcaer founders (i.e. filtered by pcaer)
                               sep='\t')
eur_pca_founders = eur_pca_founders.rename(columns={'s':'IID','fam_id':'FID','pat_id':'PAT',
                                                    'mat_id':'MAT'})
eur_pca, anc_dict = add_spark_reported_ancestry(df=eur_pca, 
                                               fam=preimp3_fam)        
#eur_founders_pca = eur_pca[(eur_pca.PAT.isna())&(eur_pca.MAT.isna())]
eur_filtered_founders_pca = eur_pca[eur_pca.IID.isin(eur_pca_founders.IID)]
eur_everyone_else_pca = eur_pca[~(eur_pca.IID.isin(eur_filtered_founders_pca.IID))]
for pcs in [[x,y] for x in range(1,11) for y in range(1,11) if (x-y<0 and x-y>=-1)]:
    fig,ax=plt.subplots(figsize=(1.5*6,1.5*4))
    ax.plot(eur_filtered_founders_pca[f'pc{pcs[0]}'],eur_filtered_founders_pca[f'pc{pcs[1]}'],'o',alpha=0.25,markeredgecolor='None')
    ax.plot(eur_everyone_else_pca[f'pc{pcs[0]}'],eur_everyone_else_pca[f'pc{pcs[1]}'],'o',alpha=0.25,markeredgecolor='None')
    ax.legend(['Founders used for PCA','Non-founders + founders not in PCA'])
    minPCx = eur_pca[f'pc{pcs[0]}'].min()
    maxPCx = eur_pca[f'pc{pcs[0]}'].max()
    rangePCx = maxPCx-minPCx
    minPCy = eur_pca[f'pc{pcs[1]}'].min()
    maxPCy = eur_pca[f'pc{pcs[1]}'].max()
    rangePCy = maxPCy-minPCy
    plt.xlim([minPCx-rangePCx*0.05, maxPCx+rangePCx*0.05])
    plt.ylim([minPCy-rangePCy*0.05, maxPCy+rangePCy*0.05])
    title_str =  f'PC{pcs[0]} vs PC{pcs[1]}\n(PCA founders used for PCA: {len(eur_filtered_founders_pca)}, Non-founders + founders not in PCA: {len(eur_everyone_else_pca)}, total: {len(eur_pca)})'
    plt.title(title_str)
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.savefig(preimp7_wd+f'spark.eur.projected.pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)






eigval = pd.read_csv('/Users/nbaya/Downloads/AFR_22_postqc.eigenval',delim_whitespace=True, names=['eigval'])
plt.plot(range(1,11), eigval.eigval.values/eigval.eigval.values.sum(),'.-')
plt.xlabel('PC')
plt.ylabel('Percent variance explained')
plt.xticks(np.arange(1, 11, step=1))
plt.xlim([0, 11])
plt.savefig('/Users/nbaya/Downloads/eigvals.png',dpi=300)

df_pcs = pd.read_csv('/Users/nbaya/Downloads/AFR_22_postqc.eigenvec',delim_whitespace=True, 
                  names=['FID','IID']+[f'PC{pc}' for pc in range(1,11)])
for pcs in [[x,y] for x in range(7,11) for y in range(7,11) if (x-y<0 and x-y>=-1)]:
    fig, ax = plt.subplots(figsize=(6,4))
    plt.plot(df_pcs[f'PC{pcs[0]}'], df_pcs[f'PC{pcs[1]}'],'.')
    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}')
    plt.tight_layout()
    plt.savefig(f'/Users/nbaya/Downloads/pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)



df_pcs = pd.read_csv('/Users/nbaya/Downloads/NeuroGAP_pilotData_clean.eigenvec',delim_whitespace=True, 
                  names=['FID','IID']+[f'PC{pc}' for pc in range(1,15)])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
pops = sorted(list(set(df_pcs.index)))
for pcs in [[x,y] for x in range(1,4) for y in range(1,4) if (x-y<0 and x-y>=-1)]:
    fig, ax = plt.subplots(figsize=(6,4))
    for pop in pops:
        plt.plot(df_pcs.loc[df_pcs.index==pop, f'PC{pcs[0]}'], 
                 df_pcs.loc[df_pcs.index==pop,f'PC{pcs[1]}'],'.')

    plt.xlabel(f'PC{pcs[0]}')
    plt.ylabel(f'PC{pcs[1]}')
    plt.title(f'PC{pcs[0]} vs PC{pcs[1]}')
    plt.tight_layout()
#    plt.savefig(f'/Users/nbaya/Downloads/pc{pcs[0]}_pc{pcs[1]}.png',dpi=300)
plt.legend(pops)





# plot GWAS from post-imputation SPARK unrelateds
# on Broad server: /stanley/genetics/users/nbaya/spark/array_May2019/spark_ricopili/gwas_check/plink_gwas/spark_plink_gwas_check.assoc.pval_lt_0.01.tsv.gz
spark = pd.read_csv('/Users/nbaya/Downloads/spark_plink_gwas_check.assoc.tsv.gz', # NOTE: OR is with respect to A1 allele
                    compression='gzip', delim_whitespace=True,
                    dtype = {'CHR': np.int, 'BP':np.int})

# iPSYCH meta-analyzed with PGC meta
ipsych = pd.read_csv('/Users/nbaya/Downloads/iPSYCH-PGC_ASD_Nov2017.gz', # NOTE: OR is the effect of the A1 allele
                     compression='gzip', delim_whitespace=True, 
                     dtype = {'CHR': np.int, 'BP':np.int})


merged = spark.merge(ipsych, on=['CHR','SNP','BP'], suffixes=('_spark','_ipsych'))


merged.loc[(merged.A1_spark==merged.A1_ipsych)&(merged.A2_spark==merged.A2_ipsych),'OR_ipsych_flipped'] = merged.OR_ipsych
merged.loc[(merged.A1_spark==merged.A2_ipsych)&(merged.A2_spark==merged.A1_ipsych),'OR_ipsych_flipped'] = 1/merged.OR_ipsych

# get SNPs with A1/A2 (in whatever order) in common between SPARK and the iPSYCH + PGC ASD meta-analysis
merged = merged.loc[((merged.A1_spark==merged.A2_ipsych)&(merged.A2_spark==merged.A2_ipsych))|
                    ((merged.A2_spark==merged.A1_ipsych)&(merged.A1_spark==merged.A2_ipsych))]

merged[['CHR','BP','SNP','P_ipsych']].to_csv('/Users/nbaya/Downloads/spark.meta_ipsych.merged.tsv.gz',
                                              sep='\t', compression='gzip', index=None)





#merged.loc[merged.OR_flipped.isna()]
#
#merged.loc['matching'] = 0
#merged.loc[((merged.OR_spark<1)&(merged.OR_ipsych<1))|((merged.OR_spark>1)&(merged.OR_ipsych>1)),'matching'] = 1
merged.loc[merged.OR_ipsych_flipped.isna(),'matching_flipped'] = np.nan
merged.loc[~(((merged.OR_spark<1)&(merged.OR_ipsych_flipped>1))|((merged.OR_spark>1)&(merged.OR_ipsych_flipped>1)))&(~merged.OR_ipsych_flipped.isna()),'matching_flipped'] = 0
merged.loc[((merged.OR_spark<1)&(merged.OR_ipsych_flipped<1))|((merged.OR_spark>1)&(merged.OR_ipsych_flipped>1)),'matching_flipped'] = 1
#
#pruned = pd.read_csv('/Users/nbaya/Downloads/tmp4_spark.eur.cobg.filtered.uniq_fid_finalpruned.bim',
#                     names=['CHR','SNP','CM','BP','REF','ALT'], sep='\t')
pruned = pd.read_csv('/Users/nbaya/Downloads/spark.meta_ipsych.merged.clumped2.tsv.gz',
                     compression='gzip', delim_whitespace=True)

merged_pruned = merged.merge(pruned[['CHR','SNP','BP']], on=['CHR','SNP','BP'])
#merged_pruned = merged.merge(pruned[['CHR','SNP','BP']], on=['CHR','SNP','BP'])

#removed = merged.loc[(merged.P_ipsych<1e-5)&(merged.P_ipsych>8e-6)][['CHR','BP','SNP','OR_spark','P_spark','OR_ipsych_flipped','P_ipsych']]
#removed.to_csv('/Users/nbaya/Downloads/spark.removed_snps.tsv',sep='\t', index=None)

pvals = np.logspace(-1, -7, 56)
#
#match_rate = [(merged.loc[merged.P_ipsych<pval,'matching'].shape[0], merged.loc[merged.P_ipsych<pval,'matching'].mean()) for pval in pvals]
#match_rate = [(merged.loc[merged.P_ipsych<pval,'matching_flipped'].shape[0], merged.loc[merged.P_ipsych<pval,'matching_flipped'].mean()) for pval in pvals]
#match_rate = [(merged.loc[merged.P_spark<pval,'matching_flipped'].shape[0], merged.loc[merged.P_spark<pval,'matching_flipped'].mean()) for pval in pvals]
match_rate = [(merged_pruned.loc[merged_pruned.P_ipsych<pval,'matching_flipped'].shape[0], merged_pruned.loc[merged_pruned.P_ipsych<pval,'matching_flipped'].mean()) for pval in pvals]



fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.plot(-np.log10(
        pvals), [rate for rows, rate in match_rate],'.-', c=color)
ax1.set_ylabel('proportion matching', color=color)
ax1.set_xlabel(r'$-\log_{10}(p)$')
plt.ylim([0,1])
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:red'
plt.plot(-np.log10(pvals), [rows for rows, rate in match_rate],'.-', c=color)
ax2.set_ylabel('SNPs in common', color=color)
#plt.ylim([0,max([rows for rows, rate in match_rate])*1.01])
plt.yscale('log')
ax2.tick_params(axis='y', labelcolor=color)

merged.loc[(merged.matching==1)&(merged.matching_flipped==0)]

df = merged
max_pos = df.BP.max()
df['chr_pos'] = (df.CHR.astype(str)+ '.'+ df.BP.astype(str).str.zfill(len(str(max_pos)))).astype(float)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, ax = plt.subplots(figsize=(9,6))
for chrom in range(1,23):
    ax.plot(df[df.CHR==chrom].chr_pos, -np.log10(df[df.CHR==chrom].P_ipsych),'.',
             color = colors[chrom % 2])
ax.axhline(-np.log10(5e-8),color='k', ls='--')
plt.ylabel(r'$-\log_{10}(p)$')
plt.xticks(range(1,23))


## Check allele freqs against reference

# SPARK Europeans (with MAF (A1 freq) and INFO score, only SNPs with genotyping rate >0.995 retained)
spark_eur = pd.read_csv(f'{spark_wd}/spark.eur.cobg.filtered.geno_0.005.frq.info',
                         sep='\t', names = ['CHR','SNP','BP','A1','A2','FRQ_A1','INFO','EXP_A1','TYPED'])
info_thresh = 0.9
spark_eur = spark_eur[spark_eur.INFO>info_thresh]

# randomly choose an effect allele in SPARK (because A1 is always minor allele)
spark_eur['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=spark_eur.shape[0]) # randomly choose an effect allele
spark_eur['EAF'] = spark_eur.A1_is_effect_allele*(spark_eur.FRQ_A1)+(spark_eur.A1_is_effect_allele==0)*(1-spark_eur.FRQ_A1)


## Check allele freqs against CEU meta-analysis from PGC
spark_frq_chr23 = pd.read_csv('/Users/nbaya/Downloads/spark.eur.chr23.cobg.filtered.geno_0.005.frq.info', # filtered to SNPs with genotyping rate > 0.005
                              delim_whitespace=True, names=['CHR','SNP','BP','A1','A2','MAF','INFO','f1','f2'])
plt.hist(spark_frq_chr23.MAF,np.linspace(0,0.2,41))
plt.xlabel('MAF')
plt.ylabel('density')
plt.xlim([0,0.2])
plt.title(f'All SPARK chr X SNPs w/ genotyping rate > 0.005\n({spark_frq_chr23.shape[0]} total)')
plt.savefig('/Users/nbaya/Downloads/chr23_maf.geno_0.005.png',dpi=300)

spark_frq_chr23 = pd.read_csv('/Users/nbaya/Downloads/spark.eur.chr23.cobg.filtered.frq', # not filtered by genotyping rate
                              delim_whitespace=True)
plt.hist(spark_frq_chr23.MAF,np.linspace(0,0.2,41))
plt.xlabel('MAF')
plt.ylabel('density')
plt.xlim([0,0.20])
plt.title(f'All SPARK chr X SNPs\n({spark_frq_chr23.shape[0]} total)')
plt.savefig('/Users/nbaya/Downloads/chr23_maf.png',dpi=300)

pgc = pd.read_csv('/Users/nbaya/Downloads/daner_AUT_meta14.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017/daner_AUT_meta14_CEU_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv.gz', # NOTE: OR is the effect of the A1 allele
                     compression='gzip', delim_whitespace=True, 
                     names = ['CHR','BP','SNP','A1','A2','OR','lb95','ub95',
                              'effect','SE','P','FRQ_A1','INFO','n','direction'], # not sure what the 2nd to last column is, possibly n?
                     dtype = {'CHR': np.int, 'BP':np.int})


spark_pgc_chr23 = spark_frq_chr23.merge(pgc[['CHR','SNP','BP','A1','A2','FRQ_A1']], on=['CHR','SNP','BP'],
                        suffixes = ('_spark','_pgc'))

spark_pgc.loc[(spark_pgc.A1_spark==spark_pgc.A1_pgc),'FRQ_A1_pgc_flipped'] = spark_pgc.FRQ_A1_pgc
spark_pgc.loc[(spark_pgc.A1_spark==spark_pgc.A2_pgc),'FRQ_A1_pgc_flipped'] = 1-spark_pgc.FRQ_A1_pgc
spark_pgc.loc[(spark_pgc.A1_is_effect_allele==1),'EAF_pgc'] = spark_pgc.FRQ_A1_pgc_flipped
spark_pgc.loc[(spark_pgc.A1_is_effect_allele==0),'EAF_pgc'] = 1-spark_pgc.FRQ_A1_pgc_flipped

plt.plot(spark_pgc.EAF, spark_pgc.EAF_pgc, '.', alpha=0.005)
plt.plot([0,1],[0,1],'k--')
plt.plot([0.2,1],[0,0.8],'k:',alpha=0.5)
plt.plot([0,0.8],[0.2,1],'k:',alpha=0.5)
plt.xlabel('SPARK')
plt.ylabel('PGC meta-analysis (no iPSYCH)')
plt.title(f'EAF (case+control combined, # SNPs = {spark_pgc.shape[0]})')

# check allele freq against 1kg europeans
eur_1kg = pd.read_csv('/Users/nbaya/Documents/lab/genotype-qc/eur.frq',
                      delim_whitespace=True)
eur_1kg = eur_1kg.rename(columns={'MAF':'FRQ_A1'})
spark_eur1kg = spark_eur.merge(eur_1kg[['CHR','SNP','A1','A2','FRQ_A1']], on=['CHR','SNP'],
                               suffixes=('_spark','_eur1kg'))

spark_eur1kg.loc[(spark_eur1kg.A1_spark==spark_eur1kg.A1_eur1kg),'FRQ_A1_eur1kg_flipped'] = spark_eur1kg.FRQ_A1_eur1kg
spark_eur1kg.loc[(spark_eur1kg.A1_spark==spark_eur1kg.A2_eur1kg),'FRQ_A1_eur1kg_flipped'] = 1-spark_eur1kg.FRQ_A1_eur1kg
spark_eur1kg.loc[(spark_eur1kg.A1_is_effect_allele==1),'EAF_eur1kg'] = spark_eur1kg.FRQ_A1_eur1kg_flipped
spark_eur1kg.loc[(spark_eur1kg.A1_is_effect_allele==0),'EAF_eur1kg'] = 1-spark_eur1kg.FRQ_A1_eur1kg_flipped


plt.plot(spark_eur1kg.EAF, spark_eur1kg.EAF_eur1kg, '.')
plt.plot([0,1],[0,1],'k--')
plt.plot([0.2,1],[0,0.8],'k:',alpha=0.5)
plt.plot([0,0.8],[0.2,1],'k:',alpha=0.5)
plt.xlabel('SPARK')
plt.ylabel('1KG Europeans')
plt.title(f'EAF (# SNPs = {spark_eur1kg.shape[0]})')         
        
        
        
        
        
        
        
# chr X GWAS sign concordance check

# on Broad server: /stanley/genetics/users/nbaya/spark/array_May2019/spark_ricopili/gwas_check/plink_gwas/spark_plink_gwas_check.assoc.pval_lt_0.01.tsv.gz
spark_chr23 = pd.read_csv('/Users/nbaya/Downloads/spark_plink_gwas_check_chr23.assoc', # NOTE: OR is with respect to A1 allele
                          delim_whitespace=True, dtype = {'CHR': np.int, 'BP':np.int})

# iPSYCH meta-analyzed with PGC meta
pgc_chr23 = pd.read_csv('/Users/nbaya/Downloads/iPSYCH-PGC_ASD_Nov2017_chrX.gz',
                          compression='gzip', delim_whitespace=True,
                          dtype = {'CHR': np.int, 'BP':np.int})
pgc_chr23['FRQ_A1'] = (18313*pgc_chr23.FRQ_A_18313+27901*pgc_chr23.FRQ_U_27901)/(18313+27901)


merged = spark_chr23.merge(pgc_chr23, on=['CHR','SNP','BP'], suffixes=('_spark','_ipsych'))


merged.loc[(merged.A1_spark==merged.A1_ipsych)&(merged.A2_spark==merged.A2_ipsych),'OR_ipsych_flipped'] = merged.OR_ipsych
merged.loc[(merged.A1_spark==merged.A2_ipsych)&(merged.A2_spark==merged.A1_ipsych),'OR_ipsych_flipped'] = 1/merged.OR_ipsych

# get SNPs with A1/A2 (in whatever order) in common between SPARK and the iPSYCH + PGC ASD meta-analysis
merged = merged.loc[((merged.A1_spark==merged.A2_ipsych)&(merged.A2_spark==merged.A2_ipsych))|
                    ((merged.A2_spark==merged.A1_ipsych)&(merged.A1_spark==merged.A2_ipsych))]

merged[['CHR','BP','SNP','P_ipsych']].to_csv('/Users/nbaya/Downloads/spark.meta_ipsych.merged.chr23.tsv.gz',
                                              sep='\t', compression='gzip', index=None)

merged.loc[merged.OR_ipsych_flipped.isna(),'matching_flipped'] = np.nan
merged.loc[~(((merged.OR_spark<1)&(merged.OR_ipsych_flipped>1))|((merged.OR_spark>1)&(merged.OR_ipsych_flipped>1)))&(~merged.OR_ipsych_flipped.isna()),'matching_flipped'] = 0
merged.loc[((merged.OR_spark<1)&(merged.OR_ipsych_flipped<1))|((merged.OR_spark>1)&(merged.OR_ipsych_flipped>1)),'matching_flipped'] = 1


pruned = pd.read_csv('/Users/nbaya/Downloads/spark.meta_ipsych.merged.chr23.v2.clumped',
                     delim_whitespace=True)

merged_pruned = merged.merge(pruned[['CHR','SNP','BP']], on=['CHR','SNP','BP'])

max_log10p = -1
min_log10p = int(round(np.log10(merged_pruned.P_ipsych.min())))
pvals = np.logspace(max_log10p,min_log10p , (max_log10p-min_log10p)*10+1)

match_rate = [(merged.loc[merged.P_ipsych<pval,'matching_flipped'].shape[0], merged.loc[merged.P_ipsych<pval,'matching_flipped'].mean()) for pval in pvals] # get sign concordance w/o pruning
match_rate = [(merged_pruned.loc[merged_pruned.P_ipsych<pval,'matching_flipped'].shape[0], merged_pruned.loc[merged_pruned.P_ipsych<pval,'matching_flipped'].mean()) for pval in pvals]



fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.plot(-np.log10(
        pvals), [rate for rows, rate in match_rate],'.-', c=color)
ax1.set_ylabel('proportion matching', color=color)
ax1.set_xlabel(r'$-log_{10}(p)$')
plt.ylim([0,1])
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:red'
plt.plot(-np.log10(pvals), [rows for rows, rate in match_rate],'.-', c=color)
ax2.set_ylabel('SNPs in common', color=color)
#plt.ylim([0,max([rows for rows, rate in match_rate])*1.01])
plt.yscale('log')
ax2.tick_params(axis='y', labelcolor=color)


# check chr23 allele frequencies
spark_chr23 = pd.read_csv('/Users/nbaya/Downloads/spark.eur.chr23.cobg.filtered.geno_0.005.frq.info',
                          delim_whitespace=True, names = ['CHR','SNP','BP','A1','A2','FRQ_A1','INFO','EXP_A1','TYPED'])

spark_chr23['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=spark_chr23.shape[0]) # randomly choose an effect allele
spark_chr23['EAF'] = spark_chr23.A1_is_effect_allele*(spark_chr23.FRQ_A1)+(spark_chr23.A1_is_effect_allele==0)*(1-spark_chr23.FRQ_A1)


plt.hist(pgc_chr23.FRQ_A1,50)
plt.title('PGC + iPSYCH ASD meta-analysis (chr X)')
plt.xlabel('A1 freq (combined case/control)')
plt.savefig('/Users/nbaya/Downloads/af_spectrum.png',dpi=300)

spark_pgc = spark_chr23.merge(pgc_chr23[['CHR','SNP','BP','A1','A2','FRQ_A1']], on=['CHR','SNP'],
                        suffixes = ('_spark','_pgc'))

spark_pgc.loc[(spark_pgc.A1_spark==spark_pgc.A1_pgc),'FRQ_A1_pgc_flipped'] = spark_pgc.FRQ_A1_pgc
spark_pgc.loc[(spark_pgc.A1_spark==spark_pgc.A2_pgc),'FRQ_A1_pgc_flipped'] = 1-spark_pgc.FRQ_A1_pgc
spark_pgc.loc[(spark_pgc.A1_is_effect_allele==1),'EAF_pgc'] = spark_pgc.FRQ_A1_pgc_flipped
spark_pgc.loc[(spark_pgc.A1_is_effect_allele==0),'EAF_pgc'] = 1-spark_pgc.FRQ_A1_pgc_flipped

plt.plot(spark_pgc.EAF, spark_pgc.EAF_pgc, '.', alpha=0.005)
plt.plot([0,1],[0,1],'k--')
plt.plot([0.2,1],[0,0.8],'k:',alpha=0.5)
plt.plot([0,0.8],[0.2,1],'k:',alpha=0.5)
plt.xlabel('SPARK')
plt.ylabel('PGC meta-analysis')
plt.title(f'EAF (case+control combined, # SNPs = {spark_pgc.shape[0]}), chr23')


# all SNPs in chr23 SPARK imputation
spark_chr23 = pd.read_csv('/Users/nbaya/Downloads/spark.eur.chr23.cobg.filtered.frq',
                          delim_whitespace=True)
spark_chr23['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=spark_chr23.shape[0]) # randomly choose an effect allele
spark_chr23['EAF'] = spark_chr23.A1_is_effect_allele*(spark_chr23.MAF)+(spark_chr23.A1_is_effect_allele==0)*(1-spark_chr23.MAF)

plt.hist(spark_chr23.EAF,50)
plt.title('SPARK (chr X)')
plt.xlabel('"Effect" allele freq\n(effect allele chosen randomly)')
plt.savefig('/Users/nbaya/Downloads/spark.af_spectrum.png',dpi=300)


# males
spark_chr23 = pd.read_csv('/Users/nbaya/Downloads/males.frq',
                          delim_whitespace=True)
spark_chr23['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=spark_chr23.shape[0]) # randomly choose an effect allele
spark_chr23['EAF'] = spark_chr23.A1_is_effect_allele*(spark_chr23.MAF)+(spark_chr23.A1_is_effect_allele==0)*(1-spark_chr23.MAF)
plt.hist(spark_chr23.EAF,50)
plt.title('SPARK males (chr X)')
plt.xlabel('"Effect" allele freq\n(effect allele chosen randomly)')
plt.savefig('/Users/nbaya/Downloads/spark.af_spectrum.males.png',dpi=300)

# females
spark_chr23 = pd.read_csv('/Users/nbaya/Downloads/females.frq',
                          delim_whitespace=True)
spark_chr23['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=spark_chr23.shape[0]) # randomly choose an effect allele
spark_chr23['EAF'] = spark_chr23.A1_is_effect_allele*(spark_chr23.MAF)+(spark_chr23.A1_is_effect_allele==0)*(1-spark_chr23.MAF)
plt.hist(spark_chr23.EAF,50)
plt.title('SPARK females (chr X)')
plt.xlabel('"Effect" allele freq\n(effect allele chosen randomly)')
plt.savefig('/Users/nbaya/Downloads/spark.af_spectrum.females.png',dpi=300)

# preimp
spark_chr23 = pd.read_csv('/Users/nbaya/Downloads/preimp.frq',
                          delim_whitespace=True)
spark_chr23['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=spark_chr23.shape[0]) # randomly choose an effect allele
spark_chr23['EAF'] = spark_chr23.A1_is_effect_allele*(spark_chr23.MAF)+(spark_chr23.A1_is_effect_allele==0)*(1-spark_chr23.MAF)
plt.hist(spark_chr23.EAF,50)
plt.title('Preimputation SPARK (chr X)')
plt.xlabel('"Effect" allele freq\n(effect allele chosen randomly)')
plt.savefig('/Users/nbaya/Downloads/spark.af_spectrum.preimp.png',dpi=300)

plt.hist(spark_chr23.EAF,50)
plt.title('SPARK (chr X)')
plt.xlabel('"Effect" allele freq\n(effect allele chosen randomly)')
plt.savefig('/Users/nbaya/Downloads/spark.af_spectrum.png',dpi=300)

# 1KG EUR (CEU, FIN, GBR, IBS, TSI) chrX frequencies
eur_1kg = pd.read_csv('/Users/nbaya/Downloads/1KG_EUR.chrX.frq_alt',delim_whitespace=True,
                      names = ['CHR','SNP','BP','A1','A2','MAF'])
eur_1kg['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=eur_1kg.shape[0]) # randomly choose an effect allele
eur_1kg['EAF'] = eur_1kg.A1_is_effect_allele*(eur_1kg.MAF)+(eur_1kg.A1_is_effect_allele==0)*(1-eur_1kg.MAF)
eur_1kg = eur_1kg[eur_1kg.MAF>0.05]
plt.hist(eur_1kg.EAF, 50)
plt.title('1KG EUR (chr X, sample size=503, MAF>0.05)')
plt.xlabel('"Effect" allele freq\n(effect allele chosen randomly)')


# Compare 1KG EUR to PGC meta-analysis
pgc_chr23 = pd.read_csv('/Users/nbaya/Downloads/iPSYCH-PGC_ASD_Nov2017_chrX.gz',
                          compression='gzip', delim_whitespace=True,
                          dtype = {'CHR': np.int, 'BP':np.int})
pgc_chr23['FRQ_A1'] = (18313*pgc_chr23.FRQ_A_18313+27901*pgc_chr23.FRQ_U_27901)/(18313+27901)
eur_1kg = pd.read_csv('/Users/nbaya/Downloads/1KG_EUR.chrX.frq_alt',delim_whitespace=True,
                      names = ['CHR','SNP','BP','A1','A2','MAF'])
eur_1kg['A1_is_effect_allele'] = np.random.randint(low=0,high=2,size=eur_1kg.shape[0]).astype('bool') # randomly choose an effect allele
eur_1kg.loc[eur_1kg.A1_is_effect_allele, 'effect_allele_1kg'] = eur_1kg.A1
eur_1kg.loc[~eur_1kg.A1_is_effect_allele, 'effect_allele_1kg'] = eur_1kg.A2
eur_1kg.loc[eur_1kg.A1_is_effect_allele, 'EAF_1kg'] = eur_1kg.MAF
eur_1kg.loc[~eur_1kg.A1_is_effect_allele, 'EAF_1kg'] = 1-eur_1kg.MAF

merged_chr23 = pgc_chr23.merge(eur_1kg, on=['CHR','BP'],
                               suffixes=('_pgc','_1kg'))
merged_chr23 = merged_chr23[((merged_chr23.A1_pgc==merged_chr23.A1_1kg)&
                             (merged_chr23.A2_pgc==merged_chr23.A2_1kg))|
                            ((merged_chr23.A1_pgc==merged_chr23.A2_1kg)&
                             (merged_chr23.A2_pgc==merged_chr23.A1_1kg))]
merged_chr23['EAF_pgc'] = (merged_chr23.A1_pgc==merged_chr23.effect_allele_1kg)*merged_chr23.FRQ_A1 + (merged_chr23.A1_pgc!=merged_chr23.effect_allele_1kg)*(1-merged_chr23.FRQ_A1)

plt.plot(merged_chr23.EAF_1kg, merged_chr23.EAF_pgc, '.', alpha=0.005)
plt.plot([0,1],[0,1],'k--')
plt.plot([0.2,1],[0,0.8],'k:',alpha=0.5)
plt.plot([0,0.8],[0.2,1],'k:',alpha=0.5)
plt.xlabel('1KG Europeans')
plt.ylabel('PGC meta-analysis')
plt.title(f'EAF\n("effect" allele chosen randomly, total SNPs = {merged_chr23.shape[0]})')


