#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:02:36 2019

Get maximal independent set for parents in SPARK dataset

@author: nbaya
"""

import pandas as pd
import datetime
import numpy as np

wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/preimp3/'

ibd = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.parents.genome.gz',
                      delim_whitespace=True,compression='gzip')
ibd = ibd[['FID1','IID1','FID2','IID2','PI_HAT']]

fam = pd.read_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.parents.fam',
                  delimiter=' ',
                  names=['FID','IID','PID','MID','SEX','PHEN'])

def remove_related_alt(iid, remaining_ids, iid1iid2):
    #only remove individuals with one degree of relatedness to iid (neighboring nodes in network graph)
    remaining_ids.remove(iid)
    print(f'{iid}: {ct_dict[iid]}')
    rows_w_iid = [x for x in iid1iid2 if iid in x] #pairs of iids that have iid in them
    iids_in_rows = set.union(*rows_w_iid) #set of iids from rows that have iid in them
    iids_in_rows.remove(iid) #remove iid from set of iids
    related = list(iids_in_rows)
    print(related)
    remaining_ids = [x for x in remaining_ids if x not in related] #remove all related individuals to iid from remaining_ids
    return remaining_ids, iid1iid2

ibd_iid1 = ibd.IID1.tolist()
ibd_iid2 = ibd.IID2.tolist()
ibd_ids = list(set(ibd_iid1+ibd_iid2))
print('Starting to count connections...')
start_ct = datetime.datetime.now()
ct_ls = []
len_ibd_ids = len(ibd_ids)
for i, iid in enumerate(ibd_ids): #REALLY SLOW (~20-30 min?), ~7 min for preimp3
    if i%100==0:
        print(f'{round(100*i/len_ibd_ids,2)}% complete')
    ct = ibd_iid1.count(iid) + ibd_iid2.count(iid) #count occurrences of iid in both iid1 and iid2 lists
    ct_ls.append((iid,ct))
ct_ls0 = ct_ls.copy()
ct_ls_tmp = sorted(ct_ls, key=lambda x: x[1]) #sort to have people with the fewest connections at the start, and those with the most at the end. This should maximize the unrelated set, although it is slow.
ct_dict = dict(ct_ls_tmp)
print(f'Time to count connections: {round((datetime.datetime.now() - start_ct).seconds/60, 2)} minutes')


ibd_ids = [x[0] for x in ct_ls_tmp]

iid1iid2 = list(zip(ibd.IID1.tolist(), ibd.IID2.tolist()))
iid1iid2  = [set(x) for x in iid1iid2]

all_ids = set(fam.IID)
unrelated = list(set(all_ids).difference(ibd_ids)) #add people who are completely unrelated, i.e. not in .genome file, pihat<0.09375 with everyone else
remaining_ids = list(ibd_ids)
start = datetime.datetime.now()
print('Getting unrelated individuals...')
ii = 0
while len(remaining_ids)>0:
    print(f'\titeration {ii}')
    iid = remaining_ids[0] #get first element
    unrelated.append(iid)
    remaining_ids, iid1iid2 = remove_related_alt(iid=iid, 
                                                 remaining_ids=remaining_ids, 
                                                 iid1iid2=iid1iid2)
    ii += 1
print(f'Time for getting unrelated individuals: {round((datetime.datetime.now() - start).seconds/60, 2)} minutes')
#np.savetxt('unrelated',np.asarray(unrelated),fmt='%s')
n_completely_unrelated = len(set(all_ids).difference(ibd_ids))
print(f'Number of unrelated individuals: {len(unrelated)} = {n_completely_unrelated} completely unrelated + {len(unrelated)-n_completely_unrelated} from IMUS')

#check relatedness among unrelated individuals
start_check = datetime.datetime.now()
for i, iid in enumerate(unrelated):
    df_tmp = ibd[(ibd.IID1==iid)|(ibd.IID2==iid)]
    if len(df_tmp)>0:
        print(f'checking iid {iid} ({i+1} of {len(unrelated)})')
        related = set(df_tmp.IID1).union(df_tmp.IID2)
        related.remove(iid)
        assert all([related not in unrelated]), 'issue with iid {iid}'
print(f'Time for checking: {round((datetime.datetime.now() - start_check).seconds/60, 2)} minutes')



unrelated_parents = fam[fam.IID.isin(unrelated)]
unrelated_parents.to_csv(wd+'SPARK.27K.genotype.20190501.hg19_preimp3.parents.IMUS.fam',index=False,header=False,sep=' ')