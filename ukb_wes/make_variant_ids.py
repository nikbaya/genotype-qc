#!/usr/bin/env python3
#
# Make variant IDs in the format of the UKB WES team
#

import pandas as pd
import numpy as np
import subprocess

wd = '/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink'

for chrom in range(1,23):
  df = pd.read_csv(f'{wd}/test_ukb_wes_chr{chrom}.bim', sep='\t', names=['chr','varid','cm','pos','a1','a2'])

  print(df.head())

  df['a1_len'] = df.a1.str.len()
  df['a2_len'] = df.a2.str.len()

  df['varid'] = df.chr.astype(str) + ':' + df.pos.astype(str) + ':'
  df.loc[df.a1_len==df.a2_len, 'varid'] = df.varid + df.a2 + ':' + df.a1
  df.loc[df.a1_len>df.a2_len, 'varid'] = df.varid + 'I:'+(df.a1_len-df.a2_len).abs().astype(str)
  df.loc[df.a1_len<df.a2_len, 'varid'] = df.varid + 'D:'+(df.a1_len-df.a2_len).abs().astype(str)

  print(df.head())

  varids = df.groupby('varid')['varid']
  suffix = varids.cumcount()+1
  repeats = varids.transform('size')

  df['new_varid'] = np.where(repeats > 1, df['varid']+'.'+suffix.map(str).str.zfill(2), df['varid'])

  df[['chr','new_varid','cm','pos','a1','a2']].to_csv(f'{wd}/test_ukb_wes_chr{chrom}.bim_new', index=False, sep='\t', header=False)
