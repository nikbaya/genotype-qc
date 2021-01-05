#!/usr/bin/env bash
#
# This script is for merging the per-chromosome PLINK files
#
# Author: Nik Baya (2021-01-05)
#
#$ -N plink_merge
#$ -o ./plink_merge.log
#$ -e ./plink_merge.errors.log
#$ -q short.qf
#$ -l h_rt=60:00:00
#$ -pe shmem 20
#$ -V
#$ -P lindgren.prjc

wd=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink

mergelist=${wd}/mergelist.txt
out=${wd}/ukb_wes_200k

# create merge list file
if [ ! -f mergelist.txt ]; then
  for chr in {1..24}; do
    if [ $chr -eq 23 ]; then
      chr="X"
    elif [ $chr -eq 24 ]; then
      chr="Y"
    fi
    echo "ukb_wes_chr${chr}" >> ${mergelist}
  done
else
  echo "Warning: ${mergelist} already exists, skipping file creation"
fi

if [ ! -f ${out}.bed ]; then
  plink --merge-list ${mergelist} \
    --make-bed \
    --memory 12000 \
    --out ${out}
  if [ -f ${out}.bed ]; then
    echo "
else
  echo "Warning: ${out}.bed already exists, skipping PLINK merge"
fi
