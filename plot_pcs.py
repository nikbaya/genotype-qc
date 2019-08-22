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
#df = pd.read_csv(wd+'preimp3/SPARK.parents.IMUS.menv.mds',delim_whitespace=True) #NOTE: number of samples is 6082, 1128 fewer than ADMIXTURE results (compared to 7173 from the previous non-liftover, autosomes+non-autosomes run)
df = pd.read_csv(wd+'preimp3/SPARK.parents.IMUS_v2.menv.mds',delim_whitespace=True) #NOTE: Includes MAF filter >0.5%, number of samples is 6081, 1129 fewer than ADMIXTURE results (compared to 7173 from the previous non-liftover, autosomes+non-autosomes run)

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


fam = pd.read_csv(wd+'preimp3/SPARK.parents.IMUS_v2.menv.fam',delim_whitespace=True,header=None,names=['FID','IID','PAT','MAT','SEX','PHEN'])

df = df.merge(fam, on=['FID','IID'])


## View results from ADMIXTURE
admix_Q = pd.read_csv(wd+'merged_v2.4.Q',delim_whitespace=True,header=None,names=[f'Q{i}' for i in range(4)]) # ancestry fractions
admix_P = pd.read_csv(wd+'merged_v2.4.P',delim_whitespace=True,header=None,names=[f'Q{i}' for i in range(4)]) # allele frequencies
'''
Pop0 = African (afr)
Pop1 = Native American (amr)
Pop2 = Asian (asn)
Pop3 = European (eur)
'''
#fam = pd.read_csv(wd+'merged_v2.fam',delim_whitespace=True,header=None,names=['FID','IID','PAT','MAT','SEX','PHEN'])
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
#    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.minQ_{min_Q}.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)
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
    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.male_female.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)


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
    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.spark_vs_1kg.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)


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
    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.minQ_{min_Q}.PC{pcs[0]}PC{pcs[1]}.ancestry_{ancestry}.png',dpi=600)
    plt.close()
    



#plot PC loadings
PC=4
for PC in range(1,11):
    fig,ax = plt.subplots(figsize=(6*1.5,4*1.5))
    plt.plot(df[f'C{PC}']**2,'.',ms=2)
    plt.xlabel('variant index')
    plt.ylabel(f'PC{PC}^2')
    plt.title(f'PC loadings for PC{PC}')
    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.pc_loadings.PC{PC}.png',dpi=600)
    
    
    
#plot mendelian errors across variants
me_l = pd.read_csv(wd+'preimp3/SPARK.27K.genotype.20190501.hg19_preimp3.lmendel',delim_whitespace=True)
len(me_l)
plt.plot(me_l.N, '.',alpha=1)
plt.xlabel('variant index')
plt.ylabel('# of mendelian errors')

me_l.N.mean()
plt.hist(me_l.N,100)

#plot mendelian errors across FIDs
me_i = pd.read_csv(wd+'preimp3/SPARK.27K.genotype.20190501.hg19_preimp3.imendel',delim_whitespace=True)
me_i = me_i.sort_values(by='FID')

plt.plot(me_i.sort_values(by='N').N.values,'.')
plt.xlabel('individual index')
plt.ylabel('# of mendelian errors')

mean_N = me_i.groupby('FID')['N'].mean()
count_FID = me_i.groupby('FID')['N'].count()

# get FIDS of families with mean mendelian errors > 300
df = pd.DataFrame(mean_N[mean_N>300].index)
df.to_csv(wd+'preimp3/FID_avg_mendelian_gt_300.tsv',sep='\t',index=False,header=False)


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

df_child = pd.read_csv(wd+'SPARK_Collection_Version2/bghx_child.csv',sep=',')
df_child = df_child.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
df_adult = pd.read_csv(wd+'SPARK_Collection_Version2/bghx_adult.csv',sep=',')
df_sibling = pd.read_csv(wd+'SPARK_Collection_Version2/bghx_sibling.csv',sep=',')

#df_adult = df_adult.rename(columns={'subject_sp_id':'IID','family_id':'FID'})
df_all = df_child.append(df_adult)
df_all = df_all.append(df_sibling)
df_all = df_all.rename(columns={'subject_sp_id':'IID','family_id':'FID'})

eur = df_all[(df_all.race_white==1)&(df_all.race_more_than_one_calc!=1)][['IID','FID']]

preimp3_fam = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])

## get white parents ##
children_w_white_mother = df_child[(df_child.mother_race_white==1)&(df_child.mother_race_more_than_one_calc!=1)][['IID','FID']]
children_w_white_father = df_child[(df_child.father_race_white==1)&(df_child.father_race_more_than_one_calc!=1)][['IID','FID']]
white_mother_IIDs = set(preimp3_fam[preimp3_fam.IID.isin(children_w_white_mother.IID)].MAT.values)
white_father_IIDs = set(preimp3_fam[preimp3_fam.IID.isin(children_w_white_father.IID)].PAT.values)
white_mother_IIDs.remove('0')
white_father_IIDs.remove('0')
white_parent_IIDs = white_mother_IIDs.union(white_father_IIDs)
white_parents = preimp3_fam[preimp3_fam.IID.isin(white_parent_IIDs)]

# get all individuals who identify as only being white
white_individuals = preimp3_fam[preimp3_fam.IID.isin(white_parent_IIDs)|(preimp3_fam.MAT.isin(white_parent_IIDs)&preimp3_fam.PAT.isin(white_parent_IIDs))]

# get parents labeled by ancestry
anc_ls = ['asian','african_amer','native_amer','native_hawaiian','white','hispanic']
anc_dict = dict(zip(anc_ls,[None]*len(anc_ls)))
assert len(preimp3_fam.IID)==len(set(preimp3_fam.IID)), 'there are individuals with duplicate IIDs, be careful!'

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
    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.reported_ancestry.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)

for anc, iids in anc_dict.items():
#    print(f'number of {anc} parents in fam file: {len(iids)}')
    print(f'number of {anc} parents in PC results: {spark_reported[spark_reported.anc_reported==anc].shape[0]}')

    


preimp3_fam[preimp3_fam.IID.isin(anc_dict['white'])|((preimp3_fam.PAT.isin(anc_dict['white']))&(preimp3_fam.MAT.isin(anc_dict['white'])))].to_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.eur.fam', sep='\t',index=False)

#merged1 = mastertable_20190501.merge(eur, on=['IID'])

preimp3_eur = preimp3_fam[preimp3_fam.FID.isin(eur.FID)]

#df_all.subject_sp_id.shape

all_fam = pd.read_csv(wd+'preimp3/SPARK.27K.genotype.20190501.liftover.v2.autosome.fam',
                      names=['FID','IID','PAT','MAT','SEX','PHEN'],
                      delim_whitespace=True)


len(set(df_all.IID).intersection(all_fam.IID))

parents_fam = pd.read_csv(wd+'preimp3/SPARK.27K.genotype.20190501.hg19_preimp3.parents.fam',
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
    plt.savefig(wd+f'plots/spark_preimp3parentsIMUS.spark_vs_1kg.reported_ancestry.PC{pcs[0]}PC{pcs[1]}.png',dpi=600)

    

# get set of 119 individuals included in 27k mastertable from 05/01 but not in SPARK.30K.array_genotype.20190423.fam
# to check if there is overlap
mastertable_20190501 = pd.read_csv(wd+'SPARK.27K.mastertable.20190501.fam',
                                   delim_whitespace=True,
                                   header=None,
                                   names=['FID','IID','PAT','MAT','SEX','PHEN'])
fam_20190423 = pd.read_csv(wd+'SPARK.30K.array_genotype.20190423.fam',
                           delim_whitespace=True,
                           header=None,
                           names=['FID','IID','PAT','MAT','SEX','PHEN'])
len(mastertable_20190501 )
len(fam_20190423)

removed = mastertable_20190501[~mastertable_20190501.IID.isin(fam_20190423.IID)]

bad_array_markers = pd.read_csv(wd+'SPARK.27K.mastertable.20190501.with-bad-array-marker.tsv',
                                delim_whitespace=True)

removed_merged_bad_array = removed.merge(bad_array_markers, on=['FID','IID'])

master_merged_bad_array = mastertable_20190501.merge(bad_array_markers, on=['FID','IID','PAT','MAT','SEX','PHEN'])

master_merged_bad_array = master_merged_bad_array.drop('call_rate_below_0.9',axis=1)

master_merged_bad_array.to_csv(wd+'SPARK.27K.mastertable.20190501.bad_array.fam',
                               sep='\t',
                               index=False)


preimp3_fam = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])

preimp3_fam[preimp3_fam.FID.isin(df.FID)]





# change PAT/MAT for individuals identified as having wrong parents

preimp3_fam = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.fam',
                          delim_whitespace=True,
                          header=None,
                          names=['FID','IID','PAT','MAT','SEX','PHEN'])

wrong_parents = pd.read_csv(wd+'me_gt_300.wrongparents.txt',
                            delim_whitespace=True)

# remove PAT/MAT for selected individuals
preimp3_fam.loc[preimp3_fam.IID.isin(wrong_parents.IID1)|preimp3_fam.IID.isin(wrong_parents.IID2),['PAT','MAT']] = 0

preimp3_fam[['FID','IID','PAT','MAT']].to_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.parents_fixed.fam',
                                              sep='\t',
                                              index=None)


# plot mendelian errors across individuals, AFTER fixing fam file for mendelian errors > 300, and filtering for mendelian errors

me_f = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp6.fmendel',delim_whitespace=True)
me_i = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp6.imendel',delim_whitespace=True) #mendelian errors per individual across variants
me_l = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp6.lmendel',delim_whitespace=True) #mendelian errors per variant across individuals

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
#df.to_csv(wd+'preimp3/FID_avg_mendelian_gt_300.tsv',sep='\t',index=False,header=False)
mean_N.mean()

plt.plot([0,len(mean_N)],[300, 300],'k--',alpha=0.5)
plt.plot(mean_N.values,'.')
plt.xlabel('Family index')
plt.ylabel('Mean number of mendelian errors')



#compare EUR maf between SPARK and 1kg

maf = pd.read_csv(wd+'SPARK_vs_1kg.eur_maf.frq',delim_whitespace=True) #first set is from SPARK, second is 1kg eur
bim = pd.read_csv(wd+'SPARK.27K.genotype.20190501.cluster_hg19_preimp3.bim',delim_whitespace=True,header=None)
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
plt.savefig(wd+'maf_diff_hist.png',dpi=600)

plt.plot(maf['MAF.1'],maf['MAF']-maf['MAF.1'],'.',ms=1,alpha=0.5)
plt.xlabel('1kg EUR maf')
plt.ylabel('SPARK EUR maf - 1kg EUR maf')
plt.savefig(wd+'1kg_eur_maf_vs_maf_diff.png',dpi=600)

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
plt.savefig(wd+'1kg_eur_vs_spark_EAF.png',dpi=600)




maf.loc[maf['A1']==maf['A1.1'],'1kgalignedAF'] = maf['MAF.1']
maf.loc[maf['A1']!=maf['A1.1'],'1kgalignedAF'] = 1-maf['MAF.1']
plt.plot(maf['1kgalignedAF'],maf['MAF'],'.',ms=5,alpha=1)
plt.plot([0,0.5],[0,0.5],'k--',alpha=0.5)
plt.plot([0.2,0.7],[0,0.5],'k:',alpha=0.2)
plt.plot([0,0.3],[0.2,0.5],'k:',alpha=0.2)
plt.xlabel(f'1kg EUR A1 freq')
plt.ylabel(f'SPARK EUR A1 freq')
plt.savefig(wd+'1kg_eur_vs_spark_A1F.png',dpi=600)






#compare EUR AF between SPARK and 1kg

af = pd.read_csv(wd+'SPARK_vs_1kg.eur_maf.frqx',sep='\t') #first set is from SPARK, second is 1kg eur
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


