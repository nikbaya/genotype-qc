#!/usr/bin/env bash
#
# This script scatters the jobs of subsetting the original per-chrom pVCFs into the 50k and 150k cohorts
#
# Author: Nik Baya (2021-01-19)
#
#$ -N subset_wes
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/subset_wes.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/subset_wes.errors.log
#$ -q short.qf
#$ -pe shmem 20
#$ -V
#$ -P lindgren.prjc
#$ -t 1-24

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly WD="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly IN="/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr${CHR}.vcf.gz"
readonly OUT="${WD}/tmp" # output directory
readonly SAMPLES="${OUT}/samples_50k_150k_split_chr${CHR}.txt"

readonly VCF_PREFIX="ukb_wes_oqfe_pvcf_"
readonly VCF_SUFFIX="_chr${CHR}"

mkdir -p ${OUT}

if [ ! -f ${SAMPLES} ]; then
  for subset in {50k,150k}; do
    sample_list="${WD}/ukb12788_wes_${subset}.sample_ids.txt"
    n_samples=$( cat ${sample_list} | wc -l )
    paste ${sample_list} \
      <( yes ${VCF_PREFIX}${subset}${VCF_SUFFIX} | head -n ${n_samples} ) >> ${SAMPLES}
  done
fi

if [ ! -f ${OUT}/${VCF_PREFIX}${subset}${VCF_SUFFIX} ]; then
  bcftools +split \
    -S ${SAMPLES} \
    ${IN} \
    -o ${OUT} \
    --threads 20
fi
