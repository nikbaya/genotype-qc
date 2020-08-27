#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 08:11:16 2020

Code for evaluating chr X QC and imputation

@author: nbaya
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats

wd = '/Users/nbaya/Documents/lab/genotype-qc/spark'

#chrX_male_hwe = pd.read_csv(f'{wd}/chrX.males.hwe', delim_whitespace=True) # invalid because heterozygosity on chrX (non-PAR) is impossible

def get_title_str(title: str):
    title = title.replace(" ","_").replace(",","")
    title = title.replace("(","_").replace(")","_")
    title = title.replace(">","gt").replace("<","lt")
    return title

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
    plt.title(title+f'\n{f"maf≥{maf}" if maf!=None else "" }({df.shape[0]} SNPs)')
    plt.savefig(f'/Users/nbaya/Downloads/exp_obs_site_het.{title.replace(" ","_")}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}.png',dpi=300)

def plot_het_distr(df, title, test, maf=None):
    r'''
    Plot distribution of observed heterozygosity
    '''
    df = df[df.TEST==test]
    if maf!=None:
        df = df.loc[df.maf>=maf, :]            
    plt.figure(figsize=(6,4))
    plt.hist(df['O(HET)'],np.linspace(0,max(0.5, df['O(HET)'].max()),50))
    plt.xlabel('O(HET), site level')
    plt.ylabel('density')
    plt.xlim(0)
    plt.title(title+f'\n{f"maf≥{maf}" if maf!=None else "" }({df.shape[0]} SNPs)')
    plt.savefig(f'/Users/nbaya/Downloads/site_het_distr.{get_title_str(title)}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}.png',dpi=300)

def plot_hwe_pval_qq(df, title, test, maf=None, ymax=None):
    r'''
    Plot QQ plot of HWE p-values
    '''
    df = df[df.TEST==test]
    if maf!=None:
        df = df.loc[df.maf>=maf, :]            
    obs = -np.log10(df.sort_values(by='P',ascending=False).P)
    exp = -np.log10(np.linspace(start=1,stop=1/df.shape[0],num=df.shape[0]))
    plt.figure(figsize=(6,4))
    plt.plot(exp, obs, '.')
    plt.plot(*[[0,exp.max()]]*2,'k--')
    plt.xlabel('Exp(-log10(p))')
    plt.ylabel('Obs(-log10(p))')
    if ymax!=None:
        plt.ylim([-ymax/20, ymax])
    plt.title(title+f'\n{f"maf≥{maf}" if maf!=None else "" }({df.shape[0]} SNPs)')
    plt.savefig(f'/Users/nbaya/Downloads/hwe_qq.{title.replace(" ","_")}{"."+test if test != None else ""}{f".maf_{maf}" if maf!=None else ""}{f".ymax_{ymax}" if ymax!=None else ""}.png',dpi=300)


chrX_females_hwe = pd.read_csv(f'{wd}/chrX.females.hwe', delim_whitespace=True)
title = 'SPARK chr X females'
plot_exp_obs_het(df=chrX_females_hwe, title=title)
plot_het_distr(df=chrX_females_hwe, title=title)
plot_hwe_pval_qq(df=chrX_females_hwe, title=title)

ukb_X_females_hwe = pd.read_csv(f'{wd}/ukb_X.females.hwe', delim_whitespace=True)
title = 'UKB chr X females'
plot_exp_obs_het(df=ukb_X_females_hwe, title=title)
plot_het_distr(df=ukb_X_females_hwe, title=title)
plot_hwe_pval_qq(df=ukb_X_females_hwe, title=title)

ukb_XY_females_hwe = pd.read_csv(f'{wd}/ukb_XY.females.hwe', delim_whitespace=True)
title = 'UKB PAR females'
plot_exp_obs_het(df=ukb_XY_females_hwe, title=title)
plot_het_distr(df=ukb_XY_females_hwe, title=title)
plot_hwe_pval_qq(df=ukb_XY_females_hwe, title=title)

ukb_XY_males_hwe = pd.read_csv(f'{wd}/ukb_XY.males.hwe', delim_whitespace=True)
title = 'UKB PAR males'
plot_exp_obs_het(df=ukb_XY_males_hwe, title=title)
plot_het_distr(df=ukb_XY_males_hwe, title=title)
plot_hwe_pval_qq(df=ukb_XY_males_hwe, title=title)

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

def plot_a1_freq_distr(df, title):
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
    plt.savefig(f'/Users/nbaya/Downloads/{title.replace(" ","_")}.a1_frq.png',dpi=300)
    
def plt_a1_freq_comparison(females, males, title, alpha=0.5):
    r'''
    Plot A1 allele frequency of both sexes
    '''
    init_ct_f, init_ct_m = females.shape[0], males.shape[0]
    females = females[females.SNP.isin(males.SNP)]
    males = males[males.SNP.isin(females.SNP)]
    post_ct_f, post_ct_m = females.shape[0], males.shape[0]
    if post_ct_f < init_ct_f or post_ct_m < init_ct_m:
        assert post_ct_f==post_ct_m
        assert post_ct_f>0
        print(f'SNPs in common: {post_ct_f}')
    A1_freq_list = []
    for df in [females, males]:
        A1_freq = get_a1_freq(df=df)
        A1_freq_list.append(A1_freq)
    plt.figure(figsize=(6,4))
    plt.plot(*A1_freq_list, '.', alpha=alpha)
    plt.plot(*[[0, 0.5]]*2,'k--')
    plt.xlabel('A1 frequency, females')
    plt.ylabel('A1 frquency, males')
    plt.title(title)
    plt.savefig(f'/Users/nbaya/Downloads/{title.replace(" ","_").replace(",","_")}.a1_frq.png', dpi=300)
    
chrX_females_con_frqx = pd.read_csv(f'{wd}/chrX.females.controls.frqx', sep='\t')
title='SPARK chr X female controls'
plot_a1_freq_distr(df=chrX_females_con_frqx, title=title)

chrX_males_con_frqx = pd.read_csv(f'{wd}/chrX.males.controls.frqx', sep='\t')
title='SPARK chr X male controls'
plot_a1_freq_distr(df=chrX_males_con_frqx , title=title)

chrY_males_con_frqx = pd.read_csv(f'{wd}/chrY.males.controls.frqx', sep='\t')
title='SPARK chr Y male controls'
plot_a1_freq_distr(df=chrY_males_con_frqx , title=title)

title = 'SPARK chr X controls, female vs male'
plt_a1_freq_comparison(females=chrX_females_con_frqx, 
                       males=chrX_males_con_frqx,
                       title=title)

ukb_X_females_frqx = pd.read_csv(f'{wd}/ukb_X.females.frqx', sep='\t')
title='UKB chr X female'
plot_a1_freq_distr(df=ukb_X_females_frqx, title=title)

ukb_X_males_frqx = pd.read_csv(f'{wd}/ukb_X.males.frqx', sep='\t')
title='UKB chr X male'
plot_a1_freq_distr(df=ukb_X_males_frqx, title=title)

title = 'UKB chr X, female vs male'
plt_a1_freq_comparison(females=ukb_X_females_frqx, 
                       males=ukb_X_males_frqx,
                       title=title,
                       alpha=0.5)

ukb_XY_females_frqx = pd.read_csv(f'{wd}/ukb_XY.females.frqx', sep='\t')
title='UKB PAR female'
plot_a1_freq_distr(df=ukb_XY_females_frqx, title=title)

ukb_XY_males_frqx = pd.read_csv(f'{wd}/ukb_XY.males.frqx', sep='\t')
title='UKB PAR male'
plot_a1_freq_distr(df=ukb_XY_males_frqx, title=title)

title = 'UKB PAR, female vs male'
plt_a1_freq_comparison(females=ukb_XY_females_frqx, 
                       males=ukb_XY_males_frqx,
                       title=title,
                       alpha=0.5)

ukb_chrY_males_frqx = pd.read_csv(f'{wd}/ukb_Y.males.frqx', sep='\t')
title='UKB chr Y male'
plot_a1_freq_distr(df=ukb_chrY_males_frqx, title=title)



def plot_sample_het_hist(df, title):
    df['O(HET)'] = 1-  df['O(HOM)']/df['N(NM)']
    plt.figure(figsize=(6,4))
    plt.hist(df['O(HET)'],np.linspace(0,max(0.5, df['O(HET)'].max()),50))
    plt.xlabel('O(HET), sample level')
    plt.ylabel('density')
    plt.xlim(0)
    plt.title(title)
    plt.savefig(f'/Users/nbaya/Downloads/sample_het_distr.{title.replace(" ","_")}.png',dpi=300)
    
def plot_exp_obs_sample_het(df, title):
    for x in ['O','E']:
        df[f'{x}(HET)'] = 1 - df[f'{x}(HOM)']/df['N(NM)']
    plt.figure(figsize=(6,4))
    plt.plot(df['E(HET)'], df['O(HET)'], '.')
    plt.plot(*[[df['E(HET)'].min(),df['E(HET)'].max()]]*2,'k--')
    plt.xlabel('E(HET), sample level')
    plt.ylabel('O(HET), sample level')
    plt.title(title)
    plt.savefig(f'/Users/nbaya/Downloads/exp_obs_sample_het.{title.replace(" ","_")}.png',dpi=300)
    

ukb_XY_females_het = pd.read_csv(f'{wd}/ukb_XY.females.het.gz', delim_whitespace=True, compression='gzip')
title='UKB PAR female'
plot_sample_het_hist(df=ukb_XY_females_het, title=title)
plot_exp_obs_sample_het(df=ukb_XY_females_het, title=title)

ukb_XY_males_het = pd.read_csv(f'{wd}/ukb_XY.males.het.gz', delim_whitespace=True, compression='gzip')
title='UKB PAR male'
plot_sample_het_hist(df=ukb_XY_males_het, title=title)
plot_exp_obs_sample_het(df=ukb_XY_males_het, title=title)





# 7/13/20
# 1. Remove SNPs with MAF<0.01 and re-plot the SNP level heterozygosity in non-par portion of X chr in females
ukb_X_f_hwe = pd.read_csv(f'{wd}/ukb_X.females.hwe', delim_whitespace=True)
ukb_X_f_hwe_frqx = pd.read_csv(f'{wd}/ukb_X.females.frqx', sep='\t')
ukb_X_f_hwe_frqx['A1_freq'] = get_a1_freq(df=ukb_X_f_hwe_frqx )
ukb_X_f_hwe_frqx['maf'] = get_maf(df=ukb_X_f_hwe_frqx )
ukb_X_f_hwe_frqx = ukb_X_f_hwe.merge(ukb_X_f_hwe_frqx, on=['CHR','SNP'])

title = 'UKB chr X females'
plot_exp_obs_het(df=ukb_X_f_hwe_frqx, title=title, test='ALL(NP)', maf=0.01)
plot_het_distr(df=ukb_X_f_hwe_frqx, title=title, test='ALL(NP)', maf=0.01)


ukb_sex_chrom_imp_f_hwe = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_female.hwe', delim_whitespace=True)
ukb_sex_chrom_imp_f_frqx = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_female.frqx', sep='\t')
ukb_sex_chrom_imp_f_frqx['A1_freq'] = get_a1_freq(df=ukb_sex_chrom_imp_f_frqx )
ukb_sex_chrom_imp_f_frqx['maf'] = get_maf(df=ukb_sex_chrom_imp_f_frqx )
ukb_sex_chrom_imp_f_hwe_frqx = ukb_sex_chrom_imp_f_hwe.merge(ukb_sex_chrom_imp_f_frqx, on=['CHR','SNP'])
ukb_sex_chrom_imp_f_hwe_frqx = ukb_sex_chrom_imp_f_hwe_frqx[ukb_sex_chrom_imp_f_hwe_frqx.CHR==23] # filter to X chrom

title = 'UKB chr X inferred females'
plot_exp_obs_het(df=ukb_sex_chrom_imp_f_hwe_frqx, title=title, test='ALL(NP)', maf=0.01)
plot_het_distr(df=ukb_sex_chrom_imp_f_hwe_frqx, title=title, test='ALL(NP)', maf=0.01)


chrX_females_hwe = pd.read_csv(f'{wd}/chrX.females.hwe', delim_whitespace=True)
chrX_females_frqx = pd.read_csv(f'{wd}/chrX.females.controls.frqx', sep='\t')
chrX_females_frqx['A1_freq'] = get_a1_freq(df=chrX_females_frqx )
chrX_females_frqx['maf'] = get_maf(df=chrX_females_frqx )
chrX_females_hwe_frqx = chrX_females_hwe.merge(chrX_females_frqx, on=['CHR','SNP'])

title = 'SPARK chr X females'
plot_exp_obs_het(df=chrX_females_hwe_frqx, title=title, test='ALL', maf=0.01)
plot_het_distr(df=chrX_females_hwe_frqx, title=title, test='ALL', maf=0.01)

# 2. Same as 1 for males


ukb_sex_chrom_imp_m_hwe = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_male.hwe', delim_whitespace=True)
ukb_sex_chrom_imp_m_frqx = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_male.frqx', sep='\t')
ukb_sex_chrom_imp_m_frqx['A1_freq'] = get_a1_freq(df=ukb_sex_chrom_imp_m_frqx )
ukb_sex_chrom_imp_m_frqx['maf'] = get_maf(df=ukb_sex_chrom_imp_m_frqx )
ukb_sex_chrom_imp_m_hwe_frqx = ukb_sex_chrom_imp_m_hwe.merge(ukb_sex_chrom_imp_m_frqx, on=['CHR','SNP'])
ukb_sex_chrom_imp_m_hwe_frqx = ukb_sex_chrom_imp_m_hwe_frqx[ukb_sex_chrom_imp_m_hwe_frqx.CHR==23] # filter to X chrom
title = 'UKB chr X inferred males'

plot_exp_obs_het(df=ukb_sex_chrom_imp_m_hwe_frqx, title=title, test='ALL(NP)', maf=0.01)
plot_het_distr(df=ukb_sex_chrom_imp_m_hwe_frqx, title=title, test='ALL(NP)', maf=0.01)

# 3a. Make QQ plot of HWE p-val after removing SNPs with MAF<0.01 
# 3b. Add a QQ plot for UKB with max on y-axis of 20

ukb_X_f_hwe = pd.read_csv(f'{wd}/ukb_X.females.hwe', delim_whitespace=True)
ukb_X_f_frqx = pd.read_csv(f'{wd}/ukb_X.females.frqx', sep='\t')
ukb_X_f_frqx ['A1_freq'] = get_a1_freq(df=ukb_X_f_frqx )
ukb_X_f_frqx ['maf'] = get_maf(df=ukb_X_f_frqx )
ukb_X_f_hwe_frqx = ukb_X_f_hwe.merge(ukb_X_f_frqx, on='SNP')
title = 'UKB chr X females'
plot_hwe_pval_qq(ukb_X_f_hwe_frqx, title, test='ALL(NP)', maf=0.1, ymax=20)
plot_hwe_pval_qq(ukb_X_f_hwe_frqx, title, test='ALL(NP)', ymax=5)


# NOTE: HWE p-val is 1 for males
#ukb_X_m_hwe = pd.read_csv(f'{wd}/ukb_X.males.hwe', delim_whitespace=True)
#ukb_X_m_frqx = pd.read_csv(f'{wd}/ukb_X.males.frqx', sep='\t')
#ukb_X_m_frqx['A1_freq'] = get_a1_freq(df=ukb_X_m_frqx )
#ukb_X_m_frqx['maf'] = get_maf(df=ukb_X_m_frqx )
#ukb_X_m_hwe_frqx = ukb_X_m_hwe.merge(ukb_X_m_frqx, on='SNP')
#title = 'UKB chr X males'
#plot_hwe_pval_qq(ukb_X_m_hwe_frqx, title, test='ALL(NP)', maf=0.01)
#plot_hwe_pval_qq(ukb_X_m_hwe_frqx, title, test='ALL(NP)', ymax=20)


# 4. Make Fhet distributions for individuals heterozygosity

def plot_fhet_ind(df, title='', logscale=False):
    plt.figure()
    plt.hist(df.F, 50)
    if logscale:
        plt.yscale('symlog',linthrehsholdy=1)
    plt.xlabel('F')
    plt.ylabel('density')
    plt.title(title+f'\n({df.shape[0]} samples)')
    title_str = get_title_str(title)
    plt.savefig(f'/Users/nbaya/Downloads/fhet_ind.{title_str}{".logscale" if logscale else ""}.png',dpi=300)

ukb_XY_females_het = pd.read_csv(f'{wd}/ukb_XY.females.het.gz', delim_whitespace=True, compression='gzip')
plot_fhet_ind(df=ukb_XY_females_het, title='UKB PAR females')
plot_fhet_ind(df=ukb_XY_females_het, title='UKB PAR females', logscale=True)

ukb_XY_males_het = pd.read_csv(f'{wd}/ukb_XY.males.het.gz', delim_whitespace=True, compression='gzip')
plot_fhet_ind(df=ukb_XY_males_het, title='UKB PAR males')
plot_fhet_ind(df=ukb_XY_males_het, title='UKB PAR males', logscale=True)

ukb_sex_chrom_imp_f_het = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_female.het.gz', delim_whitespace=True, compression='gzip')
plot_fhet_ind(df=ukb_sex_chrom_imp_f_het, title='UKB PAR inferred females')
plot_fhet_ind(df=ukb_sex_chrom_imp_f_het, title='UKB PAR inferred females', logscale=True)

ukb_sex_chrom_imp_m_het = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_male.het.gz', delim_whitespace=True, compression='gzip')
plot_fhet_ind(df=ukb_sex_chrom_imp_m_het, title='UKB PAR inferred males')
plot_fhet_ind(df=ukb_sex_chrom_imp_m_het, title='UKB PAR inferred males', logscale=True)



# 5. Plot sample-level Fhet distributions for individuals with 0 genotype missingness (after removing SNPs >0.01 missingness), SNPs with MAF<0.01 removed
ukb_XY_zeromiss_het = pd.read_csv(f'{wd}/ukb_XY.geno_0.01.mind_0.maf_0.01.het.gz', delim_whitespace=True, compression='gzip')
title = 'UKB PAR zero missingness'
plot_fhet_ind(df=ukb_XY_zeromiss_het, title=title)
plot_fhet_ind(df=ukb_XY_zeromiss_het, title=title, logscale=True)
plot_sample_het_hist(df=ukb_XY_zeromiss_het, title=title)
plot_exp_obs_sample_het(df=ukb_XY_zeromiss_het, title=title)


ukb_XY_zeromiss_females_het = pd.read_csv(f'{wd}/ukb_XY.geno_0.01.mind_0.maf_0.01.females.het.gz', delim_whitespace=True, compression='gzip')
title = 'UKB PAR zero missingness females'
plot_fhet_ind(df=ukb_XY_zeromiss_females_het, title=title)
plot_fhet_ind(df=ukb_XY_zeromiss_females_het, title=title, logscale=True)

ukb_XY_zeromiss_males_het = pd.read_csv(f'{wd}/ukb_XY.geno_0.01.mind_0.maf_0.01.males.het.gz', delim_whitespace=True, compression='gzip')
title = 'UKB PAR zero missingness males'
plot_fhet_ind(df=ukb_XY_zeromiss_males_het, title=title)
plot_fhet_ind(df=ukb_XY_zeromiss_males_het, title=title, logscale=True)



# 6. Get Z value for males vs. females MAF (after removal of SNPs with MAF<0.01)

def get_maf_diff_t(male, female):
    for df in [male, female]:
        df['maf'] = get_maf(df)
        df['s2'] = df.maf*(1-df.maf)
        df['n'] = df[['C(HOM A1)','C(HET)','C(HOM A2)','C(HAP A1)', 'C(HAP A2)']].sum(axis=1)
    assert (male.drop_duplicates('SNP').shape==male.shape) & (female.drop_duplicates('SNP').shape==female.shape)
    merge = male.merge(female, on='SNP', suffixes=('_m','_f'))
    s2_div_n_m = merge.s2_m/merge.n_m 
    s2_div_n_f = merge.s2_f/merge.n_f 
    merge['sd_diff'] = np.sqrt(s2_div_n_m + s2_div_n_f)
    merge['df'] = ((s2_div_n_m + s2_div_n_f)**2)/((s2_div_n_m)**2/(merge.n_m) + (s2_div_n_f)**2/(merge.n_f))
    merge['t'] = (merge.maf_m-merge.maf_f)/merge.sd_diff
    merge['pval'] = 2*stats.t.cdf(x=-merge.t.abs(), df=merge.df)
    return merge

def plot_maf_diff_test_qq(merge, title='', maf=None, ymax=None):
    if maf!=None:
        for sex in ['m','f']:
            merge = merge.loc[merge[f'maf_{sex}']>=maf, :]
    zero_pval_ct = merge[merge.pval==0].shape[0]
    if zero_pval_ct>0:
        zero_pval_proxy = 1e-300
        assert zero_pval_proxy < merge[merge.pval!=0].pval.min()
        print(f'{zero_pval_ct} SNPs with p-val=0, setting pval={zero_pval_proxy}')
#        merge.loc[merge.pval==0, 'pval'] = zero_pval_proxy
    obs = -np.log10(merge.sort_values(by='pval',ascending=False).pval)
    exp = -np.log10(np.linspace(start=1,stop=1/merge.shape[0],num=merge.shape[0]))    
    plt.figure(figsize=(6,4))
    plt.plot(exp, obs, '.')
    plt.plot(*[[0,exp.max()]]*2,'k--')
    plt.xlabel('Exp(-log10(p))')
    plt.ylabel('Obs(-log10(p))')
    if ymax!=None:
        plt.ylim([-ymax/20, ymax])
    plt.title(title+f'\n{f"maf in both sexes≥{maf} " if maf!=None else "" }({merge.shape[0]} SNPs{f", {zero_pval_ct} SNPs with p-val=0" if zero_pval_ct!=0 else ""})')
    plt.savefig(f'/Users/nbaya/Downloads/maf_diff_qq.{get_title_str(title)}{f".maf_{maf}" if maf!=None else ""}{f".ymax_{ymax}" if ymax!=None else ""}.png',dpi=300)

def plot_t_manhattan(merge, yaxis='abs_z', title='', maf=None):
    if maf!=None:
        for sex in ['m','f']:
            merge = merge.loc[merge[f'maf_{sex}']>=maf, :]
    if yaxis == 'abs_z':
        y = merge.t.abs()
        xlabel = 'abs(z)'
        zero_pval_ct=0
    elif yaxis=='nlog10p':
        y = -np.log10(merge.pval)
        xlabel = '-log10(p)'
        zero_pval_ct = merge[merge.pval==0].shape[0]
        if zero_pval_ct>0:
            zero_pval_proxy = 1e-300
            assert zero_pval_proxy < merge[merge.pval!=0].pval.min()
            print(f'{zero_pval_ct} SNPs with p-val=0, setting pval={zero_pval_proxy}')
    plt.figure()
    plt.plot(merge.POS, y,'.')
    plt.xlabel('base pair position')
    plt.ylabel(xlabel)
#    plt.yscale('symlog',linthresholdy=1e-9)
    plt.title(title+f'\n{f"maf in both sexes≥{maf} " if maf!=None else "" }({merge.shape[0]} SNPs{f", {zero_pval_ct} SNPs with p-val=0" if zero_pval_ct!=0 else ""})')
    plt.savefig(f'/Users/nbaya/Downloads/maf_ttest_manhattan.{yaxis}.{get_title_str(title)}{f".maf_{maf}" if maf!=None else ""}.png',dpi=300)
    
ukb_X_f_frqx = pd.read_csv(f'{wd}/ukb_X.females.frqx', sep='\t')
ukb_X_m_frqx = pd.read_csv(f'{wd}/ukb_X.males.frqx', sep='\t')
merge = get_maf_diff_t(male=ukb_X_f_frqx, female=ukb_X_m_frqx)

plot_maf_diff_test_qq(merge=merge, title='UKB chr X')
plot_maf_diff_test_qq(merge=merge, title='UKB chr X', maf=0.01)
plot_maf_diff_test_qq(merge=merge, title='UKB chr X', ymax=20)
plot_maf_diff_test_qq(merge=merge, title='UKB chr X', maf=0.01, ymax=20)

bim = pd.read_csv(f'{wd}/ukb_snp_chrX_v2.bim', delim_whitespace=True, names=['CHR','SNP','CM','POS','A1','A2'])
merge = merge.merge(bim, on='SNP')
plot_t_manhattan(merge=merge, yaxis='abs_t', title='UKB chr X', maf=0.01)
plot_t_manhattan(merge=merge, yaxis='nlog10p', title='UKB chr X', maf=0.01)


ukb_XY_f_frqx = pd.read_csv(f'{wd}/ukb_XY.females.frqx', sep='\t')
ukb_XY_m_frqx = pd.read_csv(f'{wd}/ukb_XY.males.frqx', sep='\t')
merge = get_maf_diff_t(male=ukb_XY_f_frqx, female=ukb_XY_m_frqx)

plot_maf_diff_test_qq(merge=merge, title='UKB PAR')
plot_maf_diff_test_qq(merge=merge, title='UKB PAR', maf=0.01)


ukb_X_imp_f_frqx = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_female.frqx', sep='\t')
ukb_X_imp_m_frqx = pd.read_csv(f'{wd}/ukb_sex_chrom.imputed_male.frqx', sep='\t')
merge = get_maf_diff_t(male=ukb_X_imp_f_frqx, female=ukb_X_imp_m_frqx)
merge = merge[merge.CHR_m==23]

plot_maf_diff_test_qq(merge=merge, title='UKB X inferred sex')
plot_maf_diff_test_qq(merge=merge, title='UKB X inferred sex', maf=0.01)





# Genotype meeting 11/08/2020

# 1. Histogram of chr X site-level het
# 1.1 SPARK
for data, label in [('spark.geno',' (genotyped)'),('spark',' (imputed)')]:
    spark_chrX_hwe = pd.read_csv(f'{wd}/{data}.chrX.females.hwe', delim_whitespace=True)
    plot_het_distr(df=spark_chrX_hwe, title=f'SPARK chr X females{label}', test='ALL')
    #plot_het_distr(df=spark_chrX_hwe, title=f'SPARK chr X females', test='AFF')



# 1.2 UKB non-PAR
ukb_chrX_hwe = pd.read_csv(f'{wd}/ukb_X.hwe', delim_whitespace=True)
plot_het_distr(df=ukb_chrX_hwe, title='UKB chr X', test='ALL(NP)')


# 1.3 UKB PAR
ukb_PAR_hwe = pd.read_csv(f'{wd}/ukb_XY.hwe', delim_whitespace=True)
plot_het_distr(df=ukb_PAR_hwe, title='UKB PAR', test='ALL(NP)')




# 2. Sample-level Fhet as calculated by PLINK's --check-sex

# 2.1 SPARK chr X
# spark_chrX_sexcheck = pd.read_csv(f'{wd}/spark.chrX.sexcheck', delim_whitespace=True)

# 2.1a Females 
# only reported females
#spark_f = spark_chrX_sexcheck[spark_chrX_sexcheck.PEDSEX==2]
for data, label in [('spark.geno',' (genotyped)'),('spark',' (imputed)')]:
    spark_f = pd.read_csv(f'{wd}/{data}.chrX.females.sexcheck', delim_whitespace=True)
    plot_fhet_ind(df=spark_f, title=f'SPARK chr X reported female{label}', logscale=False)
    
    # only females confirmed w/ Fhet<0.5
    #spark_f_fhet_lt_05 = spark_f[spark_f.F<0.5]
    spark_f_fhet_lt_05 = pd.read_csv(f'{wd}/{data}.chrX.threshold1.females.sexcheck', delim_whitespace=True)
    plot_fhet_ind(df=spark_f_fhet_lt_05, title=f'SPARK chr X reported female, Fhet < 0.5{label}', logscale=False)
    
    # only females confirmed w/ Fhet<0.2
    #spark_f_fhet_lt_02 = spark_f[spark_f.F<0.2]
    spark_f_fhet_lt_02 = pd.read_csv(f'{wd}/{data}.chrX.threshold2.females.sexcheck', delim_whitespace=True)
    plot_fhet_ind(df=spark_f_fhet_lt_02, title=f'SPARK chr X reported female, Fhet < 0.2{label}', logscale=False)
    
    # 2.1b Males
    # only reported males
    #spark_m = spark_chrX_sexcheck[spark_chrX_sexcheck.PEDSEX==1]
    spark_m = pd.read_csv(f'{wd}/{data}.chrX.males.sexcheck', delim_whitespace=True)
    plot_fhet_ind(df=spark_m, title=f'SPARK chr X reported male{label}', logscale=False)
    #plot_fhet_ind(df=spark_m, title='SPARK chr X reported male', logscale=True)
    
    # only males confirmed w/ Fhet>0.5
    #spark_m_fhet_gt_05 = spark_m[spark_m.F>0.5]
    spark_m_fhet_gt_05 = pd.read_csv(f'{wd}/{data}.chrX.threshold1.males.sexcheck', delim_whitespace=True)
    plot_fhet_ind(df=spark_m_fhet_gt_05, title=f'SPARK chr X reported male, Fhet > 0.5{label}', logscale=False)
    
    # only males confirmed w/ Fhet>0.8
    #spark_m_fhet_gt_08 = spark_m[spark_m.F>0.8]
    spark_m_fhet_gt_08 = pd.read_csv(f'{wd}/{data}.chrX.threshold2.males.sexcheck', delim_whitespace=True)
    plot_fhet_ind(df=spark_m_fhet_gt_08 , title=f'SPARK chr X reported male, Fhet > 0.8{label}', logscale=False)



# 2.2 UKB chr X
ukb_chrX = pd.read_csv(f'{wd}/ukb_X.sexcheck', delim_whitespace=True)

# 2.2a Females
# only reported females
ukb_chrX_f = ukb_chrX[ukb_chrX.PEDSEX==2]
plot_fhet_ind(df=ukb_chrX_f, title='UKB chr X reported female', logscale=False)

# only females confirmed w/ Fhet<0.5
ukb_chrX_f_fhet_lt_05 = ukb_chrX_f[ukb_chrX_f.F<0.5]
plot_fhet_ind(df=ukb_chrX_f_fhet_lt_05, title='UKB chr X reported female, Fhet < 0.5', logscale=False)

# only females confirmed w/ Fhet<0.2
ukb_chrX_f_fhet_lt_02 = ukb_chrX_f[ukb_chrX_f.F<0.2]
plot_fhet_ind(df=ukb_chrX_f_fhet_lt_02, title='UKB chr X reported female, Fhet < 0.2', logscale=False)

# 2.2b Males
# only reported males
ukb_chrX_m = ukb_chrX[ukb_chrX.PEDSEX==1]
plot_fhet_ind(df=ukb_chrX_m, title='UKB chr X reported male', logscale=False)
plot_fhet_ind(df=ukb_chrX_m, title='UKB chr X reported male', logscale=True)

# NOTE: All Fhet values are 1 for UKB males


# 2.3 UKB PAR
ukb_PAR = pd.read_csv(f'{wd}/ukb_XY.sexcheck', delim_whitespace=True)

# 2.3a Females
# only reported females
ukb_PAR_f = ukb_PAR[ukb_PAR.PEDSEX==2]
plot_fhet_ind(df=ukb_PAR_f, title='UKB PAR reported female', logscale=False)

# only females confirmed w/ Fhet<0.5
ukb_PAR_f_fhet_lt_05 = ukb_PAR_f[ukb_PAR_f.F<0.5]
plot_fhet_ind(df=ukb_PAR_f_fhet_lt_05, title='UKB PAR reported female, Fhet < 0.5', logscale=False)

# only females confirmed w/ Fhet<0.2
ukb_PAR_f_fhet_lt_02 = ukb_PAR_f[ukb_PAR_f.F<0.2]
plot_fhet_ind(df=ukb_PAR_f_fhet_lt_02, title='UKB PAR reported female, Fhet < 0.2', logscale=False)

# 2.3b Males
# only reported males
ukb_PAR_m = ukb_PAR[ukb_PAR.PEDSEX==1]
plot_fhet_ind(df=ukb_PAR_m, title='UKB PAR reported male', logscale=False)

# only males confirmed with Fhet>0.5
ukb_PAR_m_fhet_gt_05 = ukb_PAR_m[ukb_PAR_m.F>0.5]
plot_fhet_ind(df=ukb_PAR_m_fhet_gt_05, title='UKB PAR reported male, Fhet > 0.5', logscale=False)

# only males confirmed w/ Fhet>0.8
ukb_PAR_m_fhet_gt_08 = ukb_PAR_m[ukb_PAR_m.F>0.8]
plot_fhet_ind(df=ukb_PAR_m_fhet_gt_08, title='UKB PAR reported male, Fhet > 0.8', logscale=False)

## Another way to get Fhet instead of running --check-sex on PAR
#ukb_PAR_het = pd.read_csv(f'{wd}/ukb_XY.het', delim_whitespace=True)
#ukb_PAR_m = ukb_PAR_het.merge(ukb_chrX_m.drop('F',axis=1), on=['FID','IID'])
#plot_fhet_ind(df=ukb_PAR_m, title='UKB PAR reported male', logscale=False)


# 3. QQ plot and Manhattan plot of Z values for allele freq differences

# 3.1 SPARK chr X
for data, label in [('spark.geno',' (genotyped)'),('spark',' (imputed)')]:
    spark_chrX_f_frqx = pd.read_csv(f'{wd}/{data}.chrX.founders.females.frqx', sep='\t')
    spark_chrX_m_frqx = pd.read_csv(f'{wd}/{data}.chrX.founders.males.frqx', sep='\t')
    merge = get_maf_diff_t(male=spark_chrX_m_frqx , female=spark_chrX_f_frqx )
    
    plot_maf_diff_test_qq(merge=merge, title=f'SPARK chr X founders{label}')
    
    bim = pd.read_csv(f'{wd}/spark.eur.chr23.{"cleaned" if data=="spark" else "geno"}.bim', delim_whitespace=True, names=['CHR','SNP','CM','POS','A1','A2'])
    merge = merge.merge(bim, on='SNP')
    
    plot_t_manhattan(merge=merge, yaxis='abs_z', title=f'SPARK chr X founders{label}')



# 3.2 UKB chr X
ukb_chrX_f_frqx = pd.read_csv(f'{wd}/ukb_X.females.frqx', sep='\t')
ukb_chrX_m_frqx = pd.read_csv(f'{wd}/ukb_X.males.frqx', sep='\t')
merge = get_maf_diff_t(male=ukb_chrX_m_frqx , female=ukb_chrX_f_frqx )

plot_maf_diff_test_qq(merge=merge, title='UKB chr X')

bim = pd.read_csv(f'{wd}/ukb_snp_chrX_v2.bim', delim_whitespace=True, names=['CHR','SNP','CM','POS','A1','A2'])
merge = merge.merge(bim, on='SNP')

plot_t_manhattan(merge=merge, yaxis='abs_z', title='UKB chr X')

# 3.3 UKB PAR
ukb_PAR_f_frqx = pd.read_csv(f'{wd}/ukb_XY.females.frqx', sep='\t')
ukb_PAR_m_frqx = pd.read_csv(f'{wd}/ukb_XY.males.frqx', sep='\t')
merge = get_maf_diff_t(male=ukb_PAR_m_frqx , female=ukb_PAR_f_frqx )

plot_maf_diff_test_qq(merge=merge, title='UKB PAR')

bim = pd.read_csv(f'{wd}/ukb_XY.cleaned.bim', delim_whitespace=True, names=['CHR','SNP','CM','POS','A1','A2'])
merge = merge.merge(bim, on='SNP')

plot_t_manhattan(merge=merge, yaxis='abs_z', title='UKB PAR')


