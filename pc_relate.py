#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 18:24:51 2019

Run PC-relate on SPARK data

@author: nbaya
"""

import hail as hl
hl.init(log='/tmp/foo.log')


wd = 'gs://qc-nbaya/spark/array_May2019/preimputation/spark_preimp7/'

bfile = wd+'SPARK.27K.genotype.20190501.hg19_preimp7.founders'
#print(f'Using bfile: {bfile}')
#mt = hl.import_plink(bed=bfile+'.bed',
#                     bim=bfile+'.bim',
#                     fam=bfile+'.fam')
#
#mt = mt.checkpoint(bfile+'.mt')

mt = hl.read_matrix_table(bfile+'.mt')

min_kinship=0.09375/2

pcrelate = hl.pc_relate(call_expr=mt.GT, min_individual_maf=0.01, k=20, 
                        min_kinship=min_kinship, statistics='kin')

ct = pcrelate.count()

print('\n############\ncount:{ct}\n############')

pcrelate.export(bfile+'.pc_relate.v2.tsv.bgz')