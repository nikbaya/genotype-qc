#!/usr/bin/env python3

import pandas as pd

id_map = pd.read_csv('/well/lindgren/UKBIOBANK/DATA/WES/bridge_11867_12788.csv', dtype={'eid_11867':str,'eid_12788':str})

df = pd.read_csv('/well/lindgren/UKBIOBANK/stefania/RelGroups/2020_10_05/QC.txt',
                 delim_whitespace=True, dtype={'fid':str, 'eid':str, 'pid':str, 'mid':str})

# create ukb11867 (Lindgren lab) ID to ukb12788 (Gil's) ID dictionary
id_dict = id_map.set_index('eid_11867').to_dict()['eid_12788']

for id_field in ['eid','pid','mid']:
  df[id_field+'_new'] = df[id_field].map(id_dict)

df = df.fillna('0')

df['SEX'] = '0'
df['PHEN'] = '0'

out='/gpfs3/well/lindgren/UKBIOBANK/nbaya/resources/ukb12788_simplified_gatk_pedigree.ped'
df = df[['fid','eid_new','pid_new','mid_new','SEX','PHEN']].rename(columns={'fid':'FID','eid_new':'IID', 'pid_new':'PAT', 'mid_new':'MAT'})
df.to_csv(out, sep='\t', index=False)
