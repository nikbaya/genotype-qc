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
        if 'MAF' not in df.columns:
            df['MAF'] = get_maf(df)
        df = df.loc[df.MAF>=maf, :]            
    obs = -np.log10(df.sort_values(by='P',ascending=False).P)
    exp = -np.log10(np.linspace(start=1,stop=1/df.shape[0],num=df.shape[0]))
    plt.figure(figsize=(6,4))
    plt.plot(exp, obs, '.')
    plt.plot(*[[0,exp.max()]]*2,'k--')
    plt.xlabel('Expected(-log10(p))')
    plt.ylabel('Observed(-log10(p))')
    if ymax!=None:
        plt.ylim([-ymax/20, ymax])
    plt.title(title+f'\n{f"MAF≥{maf}" if maf!=None else "" }({df.shape[0]} variants)')
    plt.savefig(f'{wd}/hwe_qq.{get_title_str(title)}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}{f".ymax_{ymax}" if ymax!=None else ""}.png',dpi=300,
                bbox_inches='tight')
    
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
    # plt.tight_layout()
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
    
def plot_maf_vs_hwe(hwe, title):
    zero = hwe[hwe.P==0]
    nonzero = hwe[hwe.P!=0]

    fig, axs = plt.subplots(2,1, figsize=(6,8))
    plt.xlim([0,0.5])
    axs[0].plot(nonzero.MAF, -np.log10(nonzero.P), '.', alpha=0.05, ms=0.5)
    axs[0].set_title(r'MAF vs. -log$_{10}$(HWE p-value)')
    axs[0].set_ylabel(r'-log$_{10}$(p)')
    axs[0].set_xlim([0,0.5])
    ct, bins = np.histogram(zero.MAF, np.linspace(0,0.5,21))
    axs[1].bar((bins[:-1]+bins[1:])/2, ct*(max(-np.log10(nonzero.P))/max(ct)/2), width=((bins[:-1]+bins[1:]))[0])
    axs[1].set_title('MAF distribution of HWE p-value=0 variants')
    axs[1].set_ylabel('Density')
    fig.suptitle(title, fontsize=15)
    # plt.xscale('symlog', linthreshx=1e-3)
    plt.xlabel('MAF')
    
    plt.savefig(f'{wd}/maf_vs_hwe_pval.{get_title_str(title)}.png',dpi=300,
                bbox_inches='tight')
    
def plot_hwe_indels(hwe, title):
    
    
def plot_ohet_ehet_vs_hwe(hwe, title, logy=False, logx=False):
    # print(hwe[(hwe.ratio_minus_1>-1e-1)&(hwe.ratio_minus_1<-2e-2)&(hwe.P>5e-3)&(hwe.P<2e-2)].sort_values(by='ratio_minus_1').to_string(index=False))
    # hwe = hwe[(hwe.ratio_minus_1>-1e-1)&(hwe.ratio_minus_1<-2e-2)&(hwe.P>5e-3)&(hwe.P<2e-2)].sort_values(by='ratio_minus_1')
    
    # hwe0['HOM(A1)'] = hwe0.GENO.str.split('/', expand=True)[0].astype(int)
    # hwe0['HOM(A2)'] = hwe0.GENO.str.split('/', expand=True)[2].astype(int)
    
    # hwe = hwe0[hwe0['HOM(A1)']<=10]
    # hwe = hwe0[hwe0['HOM(A2)']<=10]
    if logx:
        plt.plot(hwe['O(HET)']/hwe['E(HET)']-1, hwe.P, '.', alpha=0.04, ms=0.5)
    else:
        plt.plot(hwe['O(HET)']/hwe['E(HET)'], hwe.P, '.', alpha=0.5, ms=5)
    # plt.xlim([-10**-1,-10**-2])
    # plt.ylim([5e-3,2e-2])
    plt.ylabel('HWE p-value')
    if logy:
        plt.yscale('symlog', linthreshy=1e-5)
    if logx:
        plt.xscale('symlog', linthreshx=1e-2)
        plt.xlabel('O(HET)/E(HET)-1, variant level')
    else:
        plt.xlabel('O(HET)/E(HET), variant level')
    plt.title(title)
    # plt.savefig(f'{wd}/ohet_ehet_vs_hwe_pval.a1hom_le_10.{get_title_str(title)}{".logscaley" if logy else ""}{".logscaley" if logy else ""}.png',dpi=300,
    #             bbox_inches='tight')

    # hwe['ratio_minus_1'] = hwe['O(HET)']/hwe['E(HET)']-1
    
    plt.plot(hwe.GENO.str.split('/', expand=True)[1].astype(int), hwe['O(HET)']/hwe['E(HET)'], '.')
    plt.xscale('log')
    plt.xlabel('N(HET)')
    plt.ylabel('O(HET)/E(HET)')
    
    plt.plot(hwe.GENO.str.split('/', expand=True)[1].astype(int), hwe.ratio_minus_1, '.')
    plt.xscale('log')
    plt.xlabel('N(HET)')
    plt.ylabel('O(HET)/E(HET)-1')
    
    
    
    
    plt.savefig(f'{wd}/ohet_ehet_vs_hwe_pval.{get_title_str(title)}{".logscaley" if logy else ""}{".logscaley" if logy else ""}.png',dpi=300,
                bbox_inches='tight')


def make_all_plots(bfile, title):
    try:
        frqx = pd.read_csv(f'{wd}/{bfile}.frqx_v2.gz', sep='\t', compression='gzip')
    except:
        frqx = pd.read_csv(f'{wd}/{bfile}.frqx.gz', sep='\t', compression='gzip')
    frqx['MAF'] = get_maf(frqx)
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
    # plot_het_distr(df=hwe, title=title, test='ALL(NP)', logy=True)
    if hwe.shape[0]==frqx.shape[0]:
        hwe = pd.concat((hwe, frqx[['MAF']]), axis=1)
    # hwe = hwe.merge(frqx, on=['CHR', 'SNP', 'A1', 'A2'])
    plot_hwe_pval_qq(df=hwe, title=title, test='ALL(NP)', ymax=None)
    

    

make_all_plots(bfile='ukb_wes_200k_pre_qc', title='UKB WES 200k (pre QC)')
make_all_plots(bfile='ukb_wes_200k', title='UKB WES 200k (post QC)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_autosomes', title='UKB WES 200k autosomes (pre QC)')

make_all_plots(bfile='ukb_wes_200k_autosomes', title='UKB WES 200k autosomes (post QC)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb', title='UKB WES 200k white British (pre QC)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb_maf0.01', title='UKB WES 200k white British (pre QC, MAF>0.01)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb_maf0.05', title='UKB WES 200k white British (pre QC, MAF>0.05)')

make_all_plots(bfile='ukb_wes_200k_pre_qc_wb_maxmaf0.01', title='UKB WES 200k white British (pre QC, MAF<0.01)')

make_all_plots(bfile='test_ukb_wes_200k_clean_wb', title='UKB WES 200k white British (post QC with F-het filter)')

make_all_plots(bfile='ukb_wes_200k_wb', title='UKB WES 200k white British (post QC)')

