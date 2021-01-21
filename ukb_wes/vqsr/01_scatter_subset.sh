#!/usr/bin/env bash
#
# This script scatters the jobs of subsetting the original per-chrom pVCFs into the 50k and 150k cohorts
#
# Author: Nik Baya (2021-01-19)
#
#$ -N scatter_subset
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/scatter_subset.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/scatter_subset.errors.log
#$ -q test.qc
#$ -V
#$ -P lindgren.prjc
#$ -t 1-2

if [ ${SGE_TASK_ID} -eq 1 ]; then
  SUBSET="50k"
elif [ ${SGE_TASK_ID} -eq 2 ]; then
  SUBSET="150k"
fi
readonly SUBSET

readonly OUT_DIR="/well/lindgren/UKBIOBANK/nbaya/resources/ukb_wes_${SUBSET}-backup" # output directory
readonly SUBSET_WES_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_subset_wes.sh"

readonly QUEUE="long.qf" # queue to use for scattered subsetting (default: long.qf)
readonly N_CORES=10 # number of cores (or "slots") to use (default: ?)

qsub \
  -q ${QUEUE} \
  -pe shmem ${N_CORES} \
  -t 1 \
  -N "_${SUBSET}_subset" \
  ${SUBSET_WES_SCRIPT} \
  ${SUBSET} \
  ${OUT_DIR} \
