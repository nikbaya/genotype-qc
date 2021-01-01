#!/usr/bin/env python3
#
# Make variant IDs in the format of the UKB WES team
#
# Author: Nikolas Baya (2020-01-01)

import argparse
import pandas as pd
import numpy as np
import subprocess

wd = '/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink'

def make_variant_ids(chr):
  df = pd.read_csv(f'{wd}/test_ukb_wes_chr{chr}.bim', sep='\t', names=['chr','varid','cm','pos','a1','a2'])

  df['a1_len'] = df.a1.str.len()
  df['a2_len'] = df.a2.str.len()

  df['varid'] = df.chr.astype(str) + ':' + df.pos.astype(str) + ':'
  df.loc[df.a1_len==df.a2_len, 'varid'] = df.varid + df.a2 + ':' + df.a1
  df.loc[df.a1_len>df.a2_len, 'varid'] = df.varid + 'I:'+(df.a1_len-df.a2_len).abs().astype(str)
  df.loc[df.a1_len<df.a2_len, 'varid'] = df.varid + 'D:'+(df.a1_len-df.a2_len).abs().astype(str)

  varids = df.groupby('varid')['varid']
  suffix = varids.cumcount()+1
  repeats = varids.transform('size')

  df['new_varid'] = np.where(repeats > 1, df['varid']+'.'+suffix.map(str).str.zfill(2), df['varid'])

  df[['chr','new_varid','cm','pos','a1','a2']].to_csv(f'{wd}/test_ukb_wes_chr{chr}.bim_new', index=False, sep='\t', header=False)


def main(args):
  if args.chr=='all':
    for chr in range(1,23):
      make_variant_ids(chr)
  else:
    make_variant_ids(args.chr)

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--chr', help='Which chromosome to run (default: "all")', default='all')
  args = parser.parse_args()

  main(args)
