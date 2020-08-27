#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:46:16 2020

Get SPARK sib pairs not in complete trios

@author: nbaya
"""

import pandas as pd


df0 = pd.read_csv('/Users/nbaya/Downloads/spark.eur.preimp_final.fam', delim_whitespace=True,
                 names=['FID','IID','PAT','MAT','SEX','PHEN'])
df0['FID'] = df0.FID.str.split('*',expand=True)[1]
parents = set(df0.PAT).union(df0.MAT).difference(['0'])

master = pd.read_csv('/Users/nbaya/Documents/lab/genotype-qc/spark/SPARK.30K.mastertable.20190208.csv')
master = master.rename(columns={'sfid':'FID','spid':'IID'})

df = df0.merge(master, on=['FID','IID'], how='left')


not_complete_trio = df[~(df.MAT.isin(parents)&df.PAT.isin(parents))]

not_parents_not_complete_trio = not_complete_trio[~not_complete_trio.IID.isin(parents)]

#not_parents_not_complete_trio = not_parents_not_complete_trio.merge(master, on=['FID','IID'])

fid_cts = not_parents_not_complete_trio.FID.value_counts()
two_or_more_fam_members = not_parents_not_complete_trio[not_parents_not_complete_trio.FID.isin(fid_cts[fid_cts>=2].index)]

# remove mothers and fathers
two_or_more_fam_members1 = two_or_more_fam_members[~two_or_more_fam_members.role.isin(['Father','Mother'])] 

# remove half siblings
two_or_more_fam_members2 = two_or_more_fam_members1[~two_or_more_fam_members1.role.isin(['Maternal_Half_Sibling','Paternal_Half_Sibling'])] 
#two_or_more_fam_members2 = two_or_more_fam_members1

# get fids of probands
proband_fids = two_or_more_fam_members2[two_or_more_fam_members2.PHEN==2].FID

# get df of unaffected siblings of probands
unaff_sib = two_or_more_fam_members2[(two_or_more_fam_members2.PHEN==1)&(two_or_more_fam_members2.FID.isin(proband_fids))]

# get df of sibling gropus with at least one proband and one unaffected sibling
sib_pairs = two_or_more_fam_members2[two_or_more_fam_members2.FID.isin(unaff_sib.FID)]

sib_pairs[['FID','IID','PAT','MAT','PHEN']]

# WARNING: Some samples from the same family have different PAT/MAT fields and 
# aren't in the master table, thus making it unclear how they are related
sib_pairs[(sib_pairs.PAT!='0')&(sib_pairs.MAT!='0')][['FID','IID','PAT','MAT','PHEN']]

sib_pairs.FID.value_counts().shape

