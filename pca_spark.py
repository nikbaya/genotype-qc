#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 08:25:02 2019

Calculate PCs for SPARK parents IMUS + HGDP (merged, post-Ricopili PCA dataset)

@author: nbaya
"""

import hail as hl
import ast
import subprocess
import requests
url = 'https://raw.githubusercontent.com/nikbaya/ldscsim/master/ldscsim.py'
r = requests.get(url).text
exec(r)

hl.init(log='/tmp/foo.log')

#wd = 'gs://qc-nbaya/spark/array_May2019/preimputation/preimp7_imus/'
wd = 'gs://qc-nbaya/spark/array_May2019/preimputation/'

def plink_to_mt(bfile):
    mt_path = bfile+'.mt'
    try: #check if the file has been created
        subprocess.check_output(
            ['gsutil','ls',mt_path]) is not None
        mt = hl.read_matrix_table(mt_path)
        print(f'\n... PLINK to mt conversion already completed ...\n')
    except BaseException:
        print(f'\n... Importing PLINK bfile ...\npath: {bfile}\n')
        mt = hl.import_plink(bed=bfile+'.bed',
                             bim=bfile+'.bim',
                             fam=bfile+'.fam')
        
        mt.write(mt_path)
        
        print(f'\n... Importing matrix table ...\npath: {bfile}.mt\n')
        mt = hl.read_matrix_table(mt_path)
    return mt

def run_pca(bfile):
    mt = plink_to_mt(bfile)
    
#    mt = mt.annotate_rows(all_gt_defined = hl.agg.stats(hl.is_defined(mt.GT.n_alt_alleles())).mean)
#    mt = mt.filter_rows(mt.all_gt_defined ==1)
    
#    mt = mt.filter_cols(mt.fam_id.matches('HGDP'),keep=False) # remove HGDP individuals
    
    print(f'\ncount: {mt.count()}\n')
    
    eigenvalues, scores, row_loadings = hl.hwe_normalized_pca(mt.GT,
                                                              k=20,
                                                              compute_loadings=True)
    
    print(f'\neigenvalues:\n{eigenvalues}\n')
    
    scores_path = bfile+'.scores.founders.tsv'
    print('\n... Exporting table with column scores ...\n')
    print(f"Exporting to: {scores_path}\n")
    scores.export(scores_path)
    
    loadings_path = bfile+'.loadings.founders.tsv'
    print('\n... Exporting table with row loadings ...\n')
    print(f"Exporting to: {loadings_path}\n")
    row_loadings.export(loadings_path)
    
def calc_scores(loadings, mt, path):
    '''
    loadings: a Hail table, keyed by locus and alleles, with a field called "loadings"
    containing arrays of PC loadings for each variant.
    mt: dataset for which to calculate PC scores, with rows keyed by locus and alleles
    '''
    mt1 = mt.annotate_rows(loadings = loadings[mt.locus,mt.alleles]['loadings'])
    mt1 = normalize_genotypes(mt1.GT)
    pc_ls = [{f'pc{pc}': hl.agg.sum(mt1.norm_gt*mt1.loadings[pc-1])} for pc in range(1,21)]
    mt2 = mt1.annotate_cols(**{k: v for d in pc_ls for k, v in d.items()})
    mt2.cols().export(path)


if __name__ == '__main__':
    
    ## Run PCA on SPARK founders IMUS + HGDP post-PCA merged dataset
#    bfile = wd+'rp_pca_hgdp_v3/pcaer_preimp7.founders.imus.hgdp_v3/preimp7.founders.imus.hgdp_v3.menv' # mixed population SPARK + HGDP
    bfile1 =  wd+'eur_pca2/pcaer_spark.eur.founders/spark.eur.founders.menv' # EUR SPARK founders
#    run_pca(bfile=bfile1)
    
    ## Convert dataset with all PLINK individuals + post-Ricopili PCA merged dataset
    ## to a matrix table
#    bfile2 = wd+'admixture_hgdp/spark.all.admixture_tmp1'
    bfile2 = wd+'eur_preimp1/qc_spark.eur.preimp3/spark.eur.preimp3.qc_snp'
    spark = plink_to_mt(bfile=bfile2)
    
    ## Calculate scores for all SPARK individuals
#    loadings_path = wd+'rp_pca_hgdp_v3/pcaer_preimp7.founders.imus.hgdp_v3/preimp7.founders.imus.hgdp_v3.menv.loadings.tsv'
    loadings_path = bfile1+'.loadings.founders.tsv'
    loadings0 = hl.import_table(loadings_path,
                                impute=True)
    loadings1 = loadings0.to_pandas()
    loadings2 = loadings1.copy()
    loadings_ls = [ast.literal_eval(x) for x in loadings1['loadings']]
    alleles_ls = [ast.literal_eval(x) for x in loadings1['alleles']]
    loadings2['loadings'] = loadings_ls
    loadings2['alleles'] = alleles_ls
    loadings3 = hl.Table.from_pandas(loadings2)
    loadings4 = loadings3.annotate(locus = hl.parse_locus(loadings3.locus))
    loadings5 = loadings4.key_by('locus','alleles')
#    pca_scores_path = wd+'admixture_hgdp/spark.all.admixture_tmp1.scores.v2.tsv.bgz'
    pca_scores_path = bfile2+'.scores.tsv.bgz'
#    
    calc_scores(loadings=loadings5,
                mt=spark,
                path=pca_scores_path)