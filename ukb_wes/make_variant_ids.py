#!/usr/bin/env python3
#
# Make variant IDs in the format of the UKB WES team
#
# Author: Nikolas Baya (2020-01-01)

import argparse
import pandas as pd
import numpy as np
import subprocess

def make_variant_ids(bfile, chr):
  df = pd.read_csv(f'{wd}/${bfile}{chr}.bim', sep='\t', names=['chr','varid','cm','pos','a1','a2'])

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

  df[['chr','new_varid','cm','pos','a1','a2']].to_csv(f'{wd}/${bfile}{chr}.bim_new', index=False, sep='\t', header=False)

def main(args):
  wd = args.wd
  if args.chr=='all':
    for chr in range(1,23):
      make_variant_ids(bfile=args.bfile, chr=chr)
  else:
    make_variant_ids(bfile=args.bfile, chr=args.chr)

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--chr', help='Which chromosome to run (default: "all")', default='all')
  parser.add_argument('--bfile', help='bfile prefix, assuming it ends with "chr"', default='ukb_wes_chr')
  parser.add_argument('--wd', help='directory containing PLINK files', default='/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink')
  args = parser.parse_args()

  main(args)
