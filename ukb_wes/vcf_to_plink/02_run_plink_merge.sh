#!/usr/bin/env bash
#
# This script is for merging the per-chromosome PLINK files
#
# Author: Nik Baya (2021-01-05)
#
#$ -N plink_merge
#$ -o ./plink_merge.log
#$ -e ./plink_merge.errors.log
#$ -q long.qf
#$ -l h_rt=60:00:00
#$ -pe shmem 20
#$ -V
#$ -P lindgren.prjc

readonly WD=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink

readonly MERGELIST=${WD}/mergelist.txt
readonly OUT=${WD}/ukb_wes_200k_pre_qc

# create merge list file
if [ ! -f ${MERGELIST} ]; then
  for chr in {1..24}; do
    if [ $chr -eq 23 ]; then
      chr="X"
    elif [ $chr -eq 24 ]; then
      chr="Y"
    fi
    echo "ukb_wes_chr${chr}" >> ${MERGELIST}
  done
else
  echo "Warning: ${MERGELIST} already exists, skipping file creation"
fi

echo -e "\nstarting plink merge (job id: ${JOB_ID}, $( date ))\n"
cat ${MERGELIST}

if [ ! -f ${OUT}.bed ]; then
  plink --merge-list ${MERGELIST} \
    --make-bed \
    --memory 12000 \
    --out ${OUT}
  if [ -f ${OUT}.bed ]; then
    echo "Error: ${OUT}.bed was not successfully written."
    exit 1
else
  echo "Warning: ${OUT}.bed already exists, skipping PLINK merge"
fi
