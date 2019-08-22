#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 07:46:28 2019

Compare hail liftover vs. Danfeng's liftover on the cluster

@author: nbaya
"""

import pandas as pd

wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/'

cluster_lift = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19.bim',
                           delim_whitespace=True,header=None,
                           names=['chr','varid','cm','bp','a1','a2'])
hail_lift = pd.read_csv(wd+'SPARK.27K.genotype.20190501.liftover.bim',
                        delim_whitespace=True,header=None,
                        names=['chr','varid','cm','bp','a1','a2'])

len(hail_lift[(hail_lift.chr != 'MT')&
              (hail_lift.chr != 'X')&
              (hail_lift.chr != 'Y')])

len(cluster_lift)

sorted([str(x) for x in list(set(hail_lift.chr))])
sorted([str(x) for x in list(set(cluster_lift.chr))])

hail_autosomes = set(hail_lift[(hail_lift.chr != 'MT')&(hail_lift.chr != 'X')&(hail_lift.chr != 'Y')].varid)
cluster_autosomes = set(cluster_lift.varid)
print(f'Number of SNPs in cluster autosomes: {len(cluster_autosomes)}') #613695
print(f'Number of SNPs in hail autosomes: {len(hail_autosomes)}') #613651

intersect = hail_autosomes.intersection(cluster_autosomes)
print(f'Number of SNPs in intersection of autosomes: {len(intersect)}') #613647

union = hail_autosomes.union(cluster_autosomes)
print(f'Number of SNPs in union of autosomes: {len(union)}') #613699

hail_diff_cluster = hail_autosomes.difference(cluster_autosomes)
print(f'Number of SNPs in hail but not cluster: {len(hail_diff_cluster)}') #4

cluster_diff_hail = cluster_autosomes.difference(hail_autosomes)
print(f'Number of SNPs in cluster but not hail: {len(cluster_diff_hail)}') #48


new_cluster_lift = cluster_lift.copy()
new_cluster_lift.loc[new_cluster_lift.varid.str.contains('GSA'),'varid'] = new_cluster_lift[new_cluster_lift.varid.str.contains('GSA')].varid.str.split('-',expand=True)[1]
new_cluster_lift.loc[new_cluster_lift.varid.str.contains('seq'),'varid'] = new_cluster_lift[new_cluster_lift.varid.str.contains('seq')].varid.str.split('-',expand=True)[1]




cluster_lift.loc[(cluster_lift.varid.str.contains('GSA'))&
                 (~cluster_lift.varid.str.contains('rs')),'GSA_match'] = (
                 cluster_lift[(cluster_lift.varid.str.contains('GSA'))&
                              (~cluster_lift.varid.str.contains('rs'))].varid.str.split('-',expand=True)[1]
                 ==(cluster_lift[(cluster_lift.varid.str.contains('GSA'))&(~cluster_lift.varid.str.contains('rs'))].chr.astype('str')+':'+
                                 cluster_lift[(cluster_lift.varid.str.contains('GSA'))&(~cluster_lift.varid.str.contains('rs'))].bp.astype('str')))

cluster_lift[cluster_lift.GSA_match==False]

hail_lift.loc[(hail_lift.varid.str.contains('GSA'))&
                 (~hail_lift.varid.str.contains('rs')),'GSA_match'] = (
                 hail_lift[(hail_lift.varid.str.contains('GSA'))&
                              (~hail_lift.varid.str.contains('rs'))].varid.str.split('-',expand=True)[1]
                 ==(hail_lift[(hail_lift.varid.str.contains('GSA'))&(~hail_lift.varid.str.contains('rs'))].chr.astype('str')+':'+
                                 hail_lift[(hail_lift.varid.str.contains('GSA'))&(~hail_lift.varid.str.contains('rs'))].bp.astype('str')))

hail_lift[hail_lift.GSA_match==False]



original = pd.read_csv(wd+'SPARK.27K.genotype.20190501.bim',
                           delim_whitespace=True,header=None,
                           names=['chr','varid','cm','bp','a1','a2']) #632015

all_original_diff_cluster = set(original.varid).difference(cluster_lift.varid)
print(f'Number of SNPs in original but not cluster: {len(all_original_diff_cluster )}') #56

all_original_diff_hail = set(original.varid).difference(hail_lift.varid)
print(f'Number of SNPs in original but not hail: {len(all_original_diff_hail )}') #56

original_autosomes = set(original[original.chr<=22].varid) #613701

original_diff_hail = original_autosomes.difference(hail_autosomes)
print(f'Number of autosomal SNPs in original but not hail: {len(original_diff_hail)}') #50

original_diff_cluster = original_autosomes.difference(cluster_autosomes)
print(f'Number of autosomal SNPs in original but not cluster: {len(original_diff_cluster)}') #6
original[original.varid.isin(original_diff_cluster)]


# check for variants on which hail and danfeng's versions _disagreed_ wrt the positions of the variants

merged = cluster_lift.merge(hail_lift, on=['varid'])

merged[(merged.bp_x!=merged.bp_y)][['varid','chr_x','bp_x','chr_y','bp_y']]


# check preimp3 version

hail_preimp3 = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hail_hg19_preimp3.bim',
                           delim_whitespace=True,header=None,
                           names=['chr','varid','cm','bp','a1','a2'])

cluster_preimp3 = pd.read_csv(wd+'SPARK.27K.genotype.20190501.cluster_hg19_preimp3.bim',
                           delim_whitespace=True,header=None,
                           names=['chr','varid','cm','bp','a1','a2'])

preimp3_cluster_diff_hail = set(cluster_preimp3.varid).difference(hail_preimp3.varid)
preimp3_hail_diff_cluster = set(hail_preimp3.varid).difference(cluster_preimp3.varid)

cluster_preimp3[cluster_preimp3.varid.isin(preimp3_cluster_diff_hail)]
hail_preimp3[hail_preimp3.varid.isin(preimp3_hail_diff_cluster)]