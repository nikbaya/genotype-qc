#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 09:29:11 2019

@author: nbaya
"""

import pandas as pd

wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/'

df = pd.read_csv(wd+'test.genome', sep=' ')

pi_hat_threshold = 0.09375 #default threshold from picopili: 0.09375

ids = list(set(df.IID1).union(set(df.IID2)))

df1 = df[df.PI_HAT > pi_hat_threshold] #subset to individuals with at least one related individual

ids1 = list(set(df1.IID1).union(set(df1.IID2))) #ids of individuals with at least one related individual


def check_related(df, id_dict, iid):
    id_dict[iid] = True
    related = []
    if len(df1[df1.IID1 == iid]) > 0:
        related += df[df.IID1 == iid].IID2.values.tolist()
    if len(df1[df1.IID2 == iid]) > 0:
        related += df[df.IID2 == iid].IID1.values.tolist()
    if len(related) > 0:
        for r in related:
            if not id_dict[r]: #if the related individual hasn't been visited yet
                id_dict = check_related(df, id_dict, r)
    return id_dict
    
%%time
print(f'number of individuals: {len(ids)}')
print(f'number of relations: {len(df1)}')
id_dict = dict(zip(ids,[False]*len(ids)))
imus_related = [] #imus for individuals with pi-hat > threshold chosen in ricopili when creating .genome file
for iid in ids:
    if not id_dict[iid]:
        imus.append(iid)
        id_dict = check_related(df1 , id_dict, iid)
        
%%time
imus_set = set(imus)
n_related = []
#check result
for iid in imus:
    related1 = df[(df.IID1==iid) & (df.PI_HAT > pi_hat_threshold)].IID2.values.tolist()
    related2 = df[(df.IID2==iid) & (df.PI_HAT > pi_hat_threshold)].IID1.values.tolist()
    related=related1+related2
    n_related.append(len(related))
#    print(f'{iid}: {len(related)}')
    assert len(set(related).intersection(imus_set))==0, f'{iid} has related individuals in IMUS'

imus_array = np.asarray(imus)

np.savetxt(f'{wd}imus.tsv',imus_array,fmt='%s', delimiter='\t')

imus1 = imus.copy()

for iid in [i for i in imus1 if 'DUP' in i]:
    related1 = df[(df.IID1==iid) & (df.PI_HAT > pi_hat_threshold)].IID2.values.tolist()
    related2 = df[(df.IID2==iid) & (df.PI_HAT > pi_hat_threshold)].IID1.values.tolist()
    related=related1+related2
    print(f'{iid}: {related}')
    if iid.split('-')[0] in related:
        idx = imus1.index(iid)
        removed_iid = imus1.pop(idx)
        imus1.insert(idx,iid.split('-')[0])
        print(f'{removed_iid} -> {iid.split("-")[0]}')
        
imus2 = imus1



# Check relatedness estimates against .fam file
'''
FID = family ID, IID = individual ID, PID = paternal ID, MID = maternal ID, 
SEX = sex (1 = male, 2 = female, 0 = unknown), 
PHEN = phenotype value (1 = control, 2 = case, -9/0/non-numeric = missing data if case/control)
'''
fam = pd.read_csv(wd+'SPARK.30K.genotype.20190206.fam',delimiter=' ',names=['FID','IID','PID','MID','SEX','PHEN'])
fam = fam.sort_values(by='FID')
males = fam[fam.SEX == 1]
females = fam[fam.SEX == 2]

all_ids = fam.IID.values.tolist()

#check all trios:
missing_pid = []
missing_mid = []
pihat_pat_ls = []
pihat_mat_ls = []
trio = fam[(fam.PID != '0') & ((fam.MID != '0'))] #extract all individuals with both parents in dataset
for iid in trio.IID.values:
    pid, mid = str(trio[trio.IID == iid].PID.values[0]), str(trio[trio.IID == iid].MID.values[0])
    pat = ibd[((ibd.IID1 == pid) & (ibd.IID2 == iid)) | ((ibd.IID1 == iid) & (ibd.IID2 == pid))]
    assert len(pat)<=1, f'{pat}'
    if len(pat) == 0:
#        print(f'PID not in dataset ({pid})')
        pihat_pat = None
        missing_pid.append(pid)
    else:
        pihat_pat = pat.PI_HAT.values[0]
        pihat_pat_ls.append(pihat_pat)
    mat = ibd[((ibd.IID1 == mid) & (ibd.IID2 == iid)) | ((ibd.IID1 == iid) & (ibd.IID2 == mid))]
    assert len(mat)<=1, f'{mat}'
    if len(mat) == 0:
#        print(f'MID not in dataset ({mid})')
        pihat_mat = None
        missing_mid.append(mid)
    else:
        pihat_mat = mat.PI_HAT.values[0]
        pihat_mat_ls.append(pihat_mat)
    if ((pihat_pat!=None)&(pihat_mat!=None)):
        if ((abs(pihat_pat-0.5)>0.05) | (abs(pihat_mat-0.5)>0.05)):
            print(f'{iid}\tFather ({pid}): {pihat_pat}\t Mother ({mid}): {pihat_mat}')

#check number of unique fathers, mothers. If both numbers are the same, there are no half siblings
len(set(trio.PID))
len(set(trio.MID))

parents = set(trio.PID).union(set(trio.MID))

parents_w_gt = set(fam.IID).intersection(parents)
    
plt.hist(np.asarray(pihat_pat_ls),100, alpha=0.5)
plt.hist(np.asarray(pihat_mat_ls),100, alpha=0.5)
plt.xlabel('pi_hat')
plt.title('Relatedness estimates')
plt.legend(['Father-child pi_hat','Mother-child pi_hat'])
plt.xlim([-0.005, max(pihat_pat_ls+pihat_mat_ls)])
fig=plt.gcf()
fig.savefig(wd+'parentchild_pihat.png',dpi=600)
np.median(np.asarray(pihat_pat_ls))
np.median(np.asarray(pihat_mat_ls))
np.mean(np.asarray(pihat_pat_ls))
np.mean(np.asarray(pihat_mat_ls))

#check siblings
complete_iids = []
sib_dict = {}
for iid in trio.IID.values:
    pid, mid = str(trio[trio.IID == iid].PID.values[0]), str(trio[trio.IID == iid].MID.values[0])
    sibs = trio[(trio.PID==pid) & (trio.MID==mid) & (trio.IID!=iid)]
    if len(sibs) > 0:
        sib_ls = sibs.IID.values.tolist()
    else:
        sib_ls =  None
    sib_dict[iid] =  


genome_file = wd+'pca1.mepr.genome'
genome_tsv = wd+'pca1.mepr.genome.tsv'
subprocess.call(["column","-t",genome_tsv,"|","cat", ">",genome_tsv])
with open(genome_tsv, 'w') as output:
    server = subprocess.Popen(["column","-t",genome_file,"|","awk","'{"+'print $1"\t"$2"\t"$3"\t"$4"\t"$10'+" }'"], stdout=output)
    server.communicate()
df_tmp = pd.read_csv(wd+genome_tsv,sep=' ')
df_tmp = pd.read_csv('/Users/nbaya/Documents/lab/genotype-qc/spark/pca1.mepr.genome.tsv',sep='\t')

process = subprocess.Popen(("column "+genome_file+" > "+genome_tsv).split(),stdout=subprocess.PIPE)
output, error  = process.communicate()

# get IMUS from all individuals (not just individuals in .genome file)
%%time
print(f'number of individuals: {len(all_ids)}')
print(f'number of relations: {len(df1)}')
all_id_dict = dict(zip(all_ids,[False]*len(all_ids)))
trueimus = []
for iid in all_ids:
    if not all_id_dict[iid]:
        trueimus.append(iid)
        all_id_dict = check_related(df1 , all_id_dict, iid)
        
def get_imus():
    pruned_snps = pcadir+'ldprune'
    subprocess.check_call([plinkx,
                           "--bfile", args.bfile,
                           "--indep-pairwise", 200, 100, 0.2,
                           "--silent",
                           "--out", pruned_snps])
    pruned_bfile = pcadir+'pruned'
    subprocess.check_call([plinkx,
                           "--bfile", args.bfile,
                           "--extract", pruned_snps+".prune.in",
                           "--silent",
                           "--make-bed", 
                           "--out", pruned_bfile])
    min_pi_hat = 0.09375 #default in Ricopili: 0.09375
    subprocess.check_call([plinkx,
                           "--bfile", pruned_bfile,
                           "--genome",
                           "--silent",
                           "--min", min_pi_hat
                           "--make-bed", 
                           "--out", pruned_bfile])
    genomefile = pruned_bfile+".genome"
    ibd = pd.read_csv(genomefile,delim_whitespace=True)
    fam = pd.read_csv(args.bfile+".fam",delimiter=' ',names=['FID','IID','PID','MID','SEX','PHEN'])
    ibd_ids = set(ibd.IID1).union(set(ibd.IID2))
    all_ids = set(fam.IID)
    imus = list(set(all_ids).difference(ibd_ids)) #add people who are completely unrelated
    
    def check_related(df, id_dict, iid):
        id_dict[iid] = True
        related = {}
        df_tmp = df[(df.IID1 == iid) | (df.IID2 == iid)]
        if len(df_tmp) > 0:
            related = set(df_tmp.IID1.values.tolist()+df_tmp.IID2.values.tolist()).difference({iid})
        if len(related) > 0:
            for r in related:
                if not id_dict[r]: #if the related individual hasn't been visited yet
                    id_dict = check_related(df, id_dict, r)
        return id_dict    
    
#    all_id_dict = dict(zip(all_ids,[False]*len(all_ids)))
    iid1_set = set(ibd)
    for iid in ibd_ids:
        if not all_id_dict[iid]:
            imus.append(iid)
            all_id_dict = check_related(df1 , all_id_dict, iid)    

%%timeit
'SP0126799' in ibd_ids

%%timeit
ibd[(ibd.IID1 != 'SP0126799') & (ibd.IID2 != 'SP0126799')]

%%timeit
ibd_ids_copy = ibd_ids.copy()
if 'SP0126799' in 
ibd_ids_copy.remove('SP0126799')

%%time
import datetime
import pandas as pd
import numpy as np

print('Reading in files...')
ibd = pd.read_csv('pruned.genome',delim_whitespace=True)
fam = pd.read_csv('pruned.fam',delimiter=' ',names=['FID','IID','PID','MID','SEX','PHEN'])
print('Finished reading in files...')
def test_check_related(iid, remaining_ids, ibd_iid1, ibd_iid2):
#    print(f'removing {iid}')
    remaining_ids.remove(iid)
    related = []
    while iid in ibd_iid1:
        idx1 = ibd_iid1.index(iid)
        ibd_iid1.remove(iid)
        related.append(ibd_iid2.pop(idx1))
    while iid in ibd_iid2:
        idx2 = ibd_iid2.index(iid)
        ibd_iid2.remove(iid)
        related.append(ibd_iid1.pop(idx2))
#    print(f'related for {iid}: {related}')
    for r in related:
        if r in remaining_ids:
            remaining_ids, ibd_iid1, ibd_iid2 = test_check_related(r, 
                                                                   remaining_ids, 
                                                                   ibd_iid1, 
                                                                   ibd_iid2)
    return remaining_ids, ibd_iid1, ibd_iid2

start = datetime.datetime.now()
print('Getting IMUS...')
ibd_iid1 = list(ibd.IID1)
ibd_iid2 = list(ibd.IID2)
ibd_ids = set(ibd_iid1).union(set(ibd_iid2))
all_ids = set(fam.IID)
imus = list(set(all_ids).difference(ibd_ids)) #add people who are completely unrelated, i.e. not in .genome file
remaining_ids = list(ibd_ids)
while len(remaining_ids)>0:
    iid = remaining_ids[0] #get first element
    imus.append(iid)
    remaining_ids, ibd_iid1, ibd_iid2 = test_check_related(iid=iid, 
                                                           remaining_ids=remaining_ids, 
                                                           ibd_iid1=ibd_iid1, 
                                                           ibd_iid2=ibd_iid2)
elapsed = datetime.datetime.now() - start
print(f'Time for getting IMUS: {round(elapsed.seconds/60, 2)} minutes')
np.savetxt('imus',np.asarray(imus),fmt='%s')




plt.hist(imiss[imiss.F_MISS>0.02].F_MISS,50)
plt.xlabel("missing rate")
plt.ylabel("count")


#check maximum PC for Europeans
pc = pd.read_csv(wd+'spark.menv.mds_cov',delim_whitespace=True)
eur = pc[pc.FID.str.contains('eur')]
max(eur.C1)
asn = pc[pc.FID.str.contains('asn')]
max(asn.C2) #threshold at -0.0244
afr = pc[pc.FID.str.contains('afr')]
amr = pc[pc.FID.str.contains('amr')]
min(amr.C2)

aut_pc = pc[pc.FID.str.contains('aut')]

#get set of all complete trios with only European parents
eur_parents = pd.read_csv(wd+'SPARK.eur_parents',delim_whitespace=True,names=['FID','IID'])  #parents who are in a trio and are european (not necessarily a trio with only european parents)
eur_children = fam[(fam.PID.isin(eur_parents.IID) & fam.MID.isin(eur_parents.IID))] #children with only european parents
eur_both_parents = fam[(fam.IID.isin(eur_children.PID)) | (fam.IID.isin(eur_children.MID))] #subset of fam with parents of european trios
eur_trios = fam[(fam.PID.isin(eur_parents.IID) & fam.MID.isin(eur_parents.IID)) | ((fam.IID.isin(eur_children.PID)) | (fam.IID.isin(eur_children.MID)))]

np.savetxt(wd+'eur_trios.keep', eur_trios[['FID','IID']].values, fmt='%s')

mix_eur_children = fam[((fam.PID.isin(eur_parents.IID) & ~fam.MID.isin(eur_parents.IID)) | #children in trios who have exactly one parent who is European
        (~fam.PID.isin(eur_parents.IID) & fam.MID.isin(eur_parents.IID)))]

#get set of Asian parents
asn_parents = aut_pc[aut_pc.C2<-0.0244][['FID','IID']]
asn_children = fam[(fam.PID.isin(asn_parents.IID) & fam.MID.isin(asn_parents.IID))]
mix_asn_children =  fam[((fam.PID.isin(asn_parents.IID) & ~fam.MID.isin(asn_parents.IID)) | 
        (~fam.PID.isin(asn_parents.IID) & fam.MID.isin(asn_parents.IID)))]
asn_both_parents = fam[(fam.IID.isin(asn_children.PID)) | (fam.IID.isin(asn_children.MID))]
asn_mother_eur_father_children =  fam[(fam.PID.isin(eur_parents.IID) & fam.MID.isin(asn_parents.IID))]
asn_father_eur_mother_children =  fam[(fam.PID.isin(asn_parents.IID) & fam.MID.isin(eur_parents.IID))]





## SPARK May 2019 version

wd = '/Users/nbaya/Documents/lab/genotype-qc/spark/'

# get parents

preimp1 = pd.read_csv(wd+'SPARK.27K.genotype.20190501_preimp1.fam',
                      delim_whitespace=True,names=['fid','iid','pat','mat','sex','phe'])

parents_ls = list(set(preimp1[preimp1.pat!='0'].pat).union(set(preimp1[preimp1.mat!='0'].mat)))

parents = preimp1[preimp1.iid.isin(parents_ls)]

parents.to_csv(wd+'SPARK.27K.genotype.20190501_preimp1parents.fam',index=False,header=False,sep=' ')

founders = preimp1[(preimp1.pat=='0')&(preimp1.mat=='0')]


# get unrelated parents

ibd = pd.read_csv(wd+'SPARK.27K.genotype.20190501_preimp1parents.pihat',
                      delim_whitespace=True)

#ibd = ibd.assign(freq_IID1=ibd.groupby('IID1')['IID1'].transform('count'))
#ibd = ibd.sort_values(by='freq_IID1',ascending=True) #sort so that individuals with many connections are left to the end of the algorithm (this should maximize the unrelated set)

fam = pd.read_csv(wd+'SPARK.27K.genotype.20190501_preimp1parents.fam',delimiter=' ',names=['FID','IID','PID','MID','SEX','PHEN'])

def remove_related(iid, remaining_ids, ibd_iid1, ibd_iid2):
    #only remove individuals with one degree of relatedness to iid (neighboring nodes in network graph)
#    print(f'removing {iid}')
    remaining_ids.remove(iid)
    related = []
#    print(f'ct of {iid} in iid1: {ibd_iid1.count(iid)}')
#    print(f'ct of {iid} in iid2: {ibd_iid2.count(iid)}')
    while iid in ibd_iid1:                  # for all iid entries in iid1
        idx1 = ibd_iid1.index(iid)          # get index of first instance of iid in ibd_iid1
        ibd_iid1.remove(iid)                # remove instance of iid from ibd_iid1
        related.append(ibd_iid2.pop(idx1))  # pop iid2 corresponding to iid instance in iid1 and append to list of relateds
    while iid in ibd_iid2:                  # for all iid entries in iid2
        idx2 = ibd_iid2.index(iid)          # get index of first instance of iid in ibd_iid2
        ibd_iid2.remove(iid)                # remove instance of iid for ibd_iid2
        related.append(ibd_iid1.pop(idx2))  # pop iid1 corresponding to iid instance in iid2 and append to list of relateds
#    print(f'related for {iid}: {related}')
    remaining_ids = [x for x in remaining_ids if x not in related] #remove all related individuals to iid from remaining_ids
    return remaining_ids, ibd_iid1, ibd_iid2

ibd_iid1 = ibd.IID1.tolist()
ibd_iid2 = ibd.IID2.tolist()
ibd_ids = list(set(ibd_iid1+ibd_iid2))
#ct_ls = []
#len_ibd_ids = len(ibd_ids)
#for i, iid in enumerate(ibd_ids):
#    if i%100==0:
#        print(f'{round(100*i/len_ibd_ids,2)}% complete')
#    ct = ibd_iid1.count(iid) + ibd_iid2.count(iid) #count occurrences of iid in both iid1 and iid2 lists
#    ct_ls.append((iid,ct))
#ct_ls0 = ct_ls.copy()
#ct_ls_tmp = sorted(ct_ls, key=lambda x: x[1]) #sort to have people with the fewest connections at the start, and those with the most at the end. This should maximize the unrelated set, although it is slow.

ibd_ids = [x[0] for x in ct_ls_tmp]

all_ids = set(fam.IID)
unrelated = list(set(all_ids).difference(ibd_ids)) #add people who are completely unrelated, i.e. not in .genome file
remaining_ids = list(ibd_ids)
start = datetime.datetime.now()
print('Getting unrelated individuals...')
ii = 0
while len(remaining_ids)>0:
    print(f'iteration {ii}')
    iid = remaining_ids[0] #get first element
    unrelated.append(iid)
    remaining_ids, ibd_iid1, ibd_iid2 = remove_related(iid=iid, 
                                                       remaining_ids=remaining_ids, 
                                                       ibd_iid1=ibd_iid1, 
                                                       ibd_iid2=ibd_iid2)
    ii += 1
elapsed = datetime.datetime.now() - start
print(f'Time for getting unrelated individuals: {round(elapsed.seconds/60, 2)} minutes')
np.savetxt('unrelated',np.asarray(unrelated),fmt='%s')

