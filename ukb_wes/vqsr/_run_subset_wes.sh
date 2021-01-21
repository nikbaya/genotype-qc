#!/usr/bin/env bash
#
# This is an "internal" script for subsetting to the 50k and 150k cohorts
#
# Author: Nik Baya (2021-01-20)
#
#$ -N _subset_wes
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/subset_wes.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/subset_wes.errors.log
#$ -q short.qf
#$ -pe shmem 40
#$ -V
#$ -P lindgren.prjc

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly SUBSET=${1?Error: _run_subset_wes.sh requires the subset cohort as 1st arg} # subset cohort options: "50k" or "150k"
readonly OUT_DIR=${2?Error: _run_subset_wes.sh requires the output directory as 2nd arg} # output directory

readonly SAMPLES="/well/lindgren/UKBIOBANK/nbaya/resources/ukb12788_wes_${SUBSET}.sample_ids.txt"
readonly IN="/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr${CHR}.vcf.gz" # input VCF
readonly OUT="${OUT_DIR}/ukb_wes_oqfe_pvcf_${SUBSET}_chr${CHR}.vcf.gz" # output VCF

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

vcf_check() {
  if [ ! -f $1 ]; then
    raise_error "$1 does not exist."
  elif [ $( bcftools view -h $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    raise_error "$1 may be truncated"
  fi
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

vcf_check ${IN}

SECONDS=0

if [ ! -f ${OUT} ]; then
  set -x
  bcftools view \
    ${IN} \
    --samples-file ${SAMPLES} \
    --force-samples \
    -Oz \
    --threads 40 \
    -o ${OUT}
  set +x
fi

vcf_check ${OUT}

duration=${SECONDS}
echo "${SUBSET} chr${CHR} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
