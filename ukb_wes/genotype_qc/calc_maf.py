#!/usr/bin/env python3

import pandas as pd
import argparse

def get_a1_freq(df):
    total_ct = 2*df['C(HOM A1)']+2*df['C(HET)']+2*df['C(HOM A2)']+df['C(HAP A1)']+df['C(HAP A2)']
    A1_ct = 2*df['C(HOM A1)']+df['C(HET)']+df['C(HAP A1)']
    A1_freq = A1_ct/total_ct
    return A1_freq

def get_maf(df):
    if 'maf' in df.columns:
        return df.maf
    if not 'A1_freq' in df.columns:
        df['A1_freq'] = get_a1_freq(df=df)
    df['not_A1_freq'] = 1-df['A1_freq']
    maf = df.apply(func=lambda x: x[['A1_freq','not_A1_freq']].min(), axis=1)
    return maf

def main(args):
    fname = args.frqx
    df = pd.read_csv(f'{wd}/{fname}', compression='gzip', sep='\t')
    df['MAF'] = get_maf(df)
    df.to_csv(f'{wd}/{fname.replace(".frqx",".frqx_v2")}', compression='gzip', index=False, sep='\t')

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--frqx', help='path to frqx file')
    args = parser.parse_args()

    main(args)
