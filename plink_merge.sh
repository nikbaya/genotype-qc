#!/usr/bin/env bash
#
# This script merges datasets together using PLINK.
# The intended application is to merge a cohort with a reference panel
# to run PCA w/ reference samples to determine continental ancestry.
#
# Example: Merge datasets named "data" and "ref"
#
# ./plink_merge.sh "data" "ref"
#
# NOTE: Check SNP IDs of both datasets to be sure that they have
# similar formats, otherwise the merge may miss some SNPs in common.
#
# NOTE: Depending on which SNPs where genotyped or imputed in the non-reference
# dataset, there may not be many SNPs shared between the two datasets.
#
# Author: Nikolas Baya (2020-09-19)

bfile1=$1 	# bfile prefix of 1st dataset
bfile2=$2 	# bfile prefix of 2nd dataset
out="merged" 	# bfile prefix of merged output

plink \
	--bfile $bfile1 \
	--bmerge $bfile2 \
	--make-bed \
	--out $out
