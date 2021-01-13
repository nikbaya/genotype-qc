#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 08:14:12 2021

Plot QC statistics for UKB WES

@author: nbaya
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import re
import matplotlib.pyplot as plt

wd = '/Users/nbaya/Downloads'

def get_title_str(title: str):
    title = title.replace(" ","_").replace(",","")
    title = title.replace("(","_").replace(")","_")
    title = title.replace(">","gt").replace("<","lt")
    title = re.sub(r'(_)(?=\1)', '', title)
    title = title.rstrip('_').lower()
    return title

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

def plot_a1_freq_distr(df, title, logy=False):
    r'''
    Plot A1 allele frequencies.
    '''
    A1_freq = get_a1_freq(df=df)
    plt.figure(figsize=(6,4))
    plt.hist(A1_freq,np.linspace(0,max(0.5, A1_freq.max()),50))
    plt.xlabel('A1 frequency')
    plt.ylabel('density')
    plt.title(title)
    plt.xlim([0, max(0.5, A1_freq.max())])
    if logy:
        plt.yscale('log')
    plt.savefig(f'{wd}/{get_title_str(title)}.a1_frq.png',dpi=300)

def plot_exp_obs_het(df, title, test, maf=None):
    r'''
    Plot expected vs. observed site-level heterozygosity
    '''
    df = df[df.TEST==test]    
    if maf!=None:
        df = df.loc[df.maf>=maf, :]            
    plt.figure(figsize=(6,4))
    plt.plot(df['E(HET)'], df['O(HET)'],'.')
    plt.plot(*[[0,df['E(HET)'].max()]]*2,'k--')
    plt.xlabel('E(HET), site level')
    plt.ylabel('O(HET), site level')
    plt.tight_layout()
    plt.title(title+f'\n{f"maf≥{maf}" if maf!=None else "" }({df.shape[0]} variants)')
    plt.savefig(f'{wd}/exp_obs_site_het.{get_title_str(title)}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}.png',dpi=300)

def plot_het_distr(df, title, test, maf=None, logy=False):
    r'''
    Plot distribution of observed heterozygosity
    '''
    df = df[df.TEST==test]
    if maf!=None:
        df = df.loc[df.maf>=maf, :]            
    plt.figure(figsize=(6,4))
    plt.hist(df['O(HET)'],np.linspace(0,max(0.5, df['O(HET)'].max()),50))
    if logy:
        plt.yscale('log')
    plt.xlabel('O(HET), site level')
    plt.ylabel('density')
    plt.xlim(0)
    plt.title(title+f'\n{f"maf≥{maf}" if maf!=None else "" }({df.shape[0]} variants)')
    plt.savefig(f'{wd}/site_het_distr.{get_title_str(title)}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}.png',dpi=300)

def plot_hwe_pval_qq(df, title, test, maf=None, ymax=None):
    r'''
    Plot QQ plot of HWE p-values
    '''
    df = df[df.TEST==test]
    if maf!=None:
        if 'maf' not in df.columns:
            df['maf'] = get_maf(df)
        df = df.loc[df.maf>=maf, :]            
    obs = -np.log10(df.sort_values(by='P',ascending=False).P)
    exp = -np.log10(np.linspace(start=1,stop=1/df.shape[0],num=df.shape[0]))
    plt.figure(figsize=(6,4))
    plt.plot(exp, obs, '.')
    plt.plot(*[[0,exp.max()]]*2,'k--')
    plt.xlabel('Expected(-log10(p))')
    plt.ylabel('Observed(-log10(p))')
    if ymax!=None:
        plt.ylim([-ymax/20, ymax])
    plt.title(title+f'\n{f"maf≥{maf}" if maf!=None else "" }({df.shape[0]} variants)')
    plt.savefig(f'{wd}/hwe_qq.{title.replace(" ","_")}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}{f".ymax_{ymax}" if ymax!=None else ""}.png',dpi=300)
    
def plot_sample_het_hist(df, title, logy=False):
    df['O(HET)'] = 1-  df['O(HOM)']/df['N(NM)']
    plt.figure(figsize=(6,4))
    plt.hist(df['O(HET)'],np.linspace(0,max(0.5, df['O(HET)'].max()),50))
    plt.xlabel('O(HET), sample level')
    plt.ylabel('density')
    plt.xlim(0)
    plt.title(title)
    if logy:
        plt.yscale('log')
    plt.savefig(f'{wd}/sample_het_distr.{get_title_str(title)}.png',dpi=300)
    
def plot_exp_obs_sample_het(df, title):
    for x in ['O','E']:
        df[f'{x}(HET)'] = 1 - df[f'{x}(HOM)']/df['N(NM)']
    plt.figure(figsize=(6,4))
    plt.plot(df['E(HET)'], df['O(HET)'], '.', alpha=0.05)
    plt.plot(*[[df['E(HET)'].min(),df['E(HET)'].max()]]*2,'k--')
    # result = stats.linregress(x=df['E(HET)'],y=df['O(HET)'])
    plt.xlabel('E(HET), sample level')
    plt.ylabel('O(HET), sample level')
    # plt.ylim([df['O(HET)'].min(), df['O(HET)'].max()])
    plt.tight_layout()
    plt.title(title)
    plt.savefig(f'{wd}/exp_obs_sample_het.{get_title_str(title)}.png',dpi=300,
                bbox_inches='tight')
    
def plot_fhet_ind(df, title, logy=False):
    plt.figure()
    plt.hist(df.F, 50)
    if logy:
        plt.yscale('log')
    plt.xlabel('F')
    plt.ylabel('density')
    plt.title(title+f'\n({df.shape[0]} samples)')
    title_str = get_title_str(title)
    plt.savefig(f'{wd}/fhet_ind.{title_str}{".logscale" if logy else ""}.png',dpi=300)
    
def plot_imiss_hist(df, title, logy=True):
    plt.figure()
    plt.hist(df.F_MISS, np.linspace(0,1,100))
    if logy:
        plt.yscale('log')
    plt.xlim([-0.01,1.01])
    plt.xlabel('Genotype missingness rate per individual')
    plt.ylabel('density')
    plt.title(title)
    plt.savefig(f'{wd}/imiss_hist.{get_title_str(title)}.png',dpi=300)

def plot_lmiss_hist(df, title, logy=True):
    plt.figure()
    plt.hist(df.F_MISS, np.linspace(0,1,100))
    if logy:
        plt.yscale('log')
    plt.xlim([-0.01,1.01])
    plt.xlabel('Genotype missingness rate per variant')
    plt.ylabel('density')
    plt.title(title)
    plt.savefig(f'{wd}/lmiss_hist.{get_title_str(title)}{".logscale" if logy else ""}.png',dpi=300)
    
def plot_maf_vs_lmiss(frqx, lmiss, title):
    plt.figure()
    merge = frqx.merge(lmiss, on=['CHR', 'SNP'])
    plt.plot(get_maf(merge), merge.F_MISS, '.', alpha=0.05)
    plt.xlabel('MAF')
    plt.ylabel('Genotype missingness rate per variant')
    plt.title(title)
    plt.savefig(f'{wd}/maf_vs_lmiss.{get_title_str(title)}.png',dpi=300)
    

def make_all_plots(bfile, title):
    frqx = pd.read_csv(f'{wd}/{bfile}.frqx.gz', sep='\t', compression='gzip')
    frqx['maf'] = get_maf(frqx)
    # plot_a1_freq_distr(df=frqx, title=title, logy=True)
    
    het = pd.read_csv(f'{wd}/{bfile}.het', delim_whitespace=True)
    plot_sample_het_hist(df=het, title=title, logy=True)
    plot_exp_obs_sample_het(df=het, title=title)
    plot_fhet_ind(df=het, title=title, logy=True)
        
    imiss = pd.read_csv(f'{wd}/{bfile}.imiss.gz', delim_whitespace=True, compression='gzip')
    plot_imiss_hist(df=imiss, title=title, logy=True)
    
    lmiss = pd.read_csv(f'{wd}/{bfile}.lmiss.gz', delim_whitespace=True, compression='gzip')
    plot_lmiss_hist(df=lmiss, title=title, logy=True)
    
    hwe = pd.read_csv(f'{wd}/{bfile}.hwe.gz', delim_whitespace=True, compression='gzip')
    plot_het_distr(df=hwe, title=title, test='ALL(NP)', logy=True)
    if hwe.shape[0]==frqx.shape[0]:
        hwe = pd.concat((hwe, frqx[['maf']]), axis=1)
    # hwe = hwe.merge(frqx, on=['CHR', 'SNP', 'A1', 'A2'])
    # plot_hwe_pval_qq(df=hwe, title=title, test='ALL(NP)', maf=0.01, ymax=None)

    

make_all_plots(bfile='ukb_wes_200k_pre_qc', title='UKB WES 200k (pre QC)')
make_all_plots(bfile='ukb_wes_200k', title='UKB WES 200k (post QC)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_autosomes', title='UKB WES 200k autosomes (pre QC)')

make_all_plots(bfile='ukb_wes_200k_autosomes', title='UKB WES 200k autosomes (post QC)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb', title='UKB WES 200k white British (pre QC)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb_maf0.01', title='UKB WES 200k white British (pre QC, MAF>0.01)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb_maf0.05', title='UKB WES 200k white British (pre QC, MAF>0.05)')
