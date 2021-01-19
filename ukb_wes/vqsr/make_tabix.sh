#!/usr/bin/env bash
#
# This script merges the per-chromosome pVCFs
#
# Author: Nik Baya (2021-01-18)
#
#$ -N make_tabix
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/make_tabix.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/make_tabix.errors.log
#$ -q long.qe
#$ -V
#$ -P lindgren.prjc

# input and output
readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/vcf"
readonly OUT="${WD}/ukb_wes_oqfe_pvcf.vcf.gz"


raise_error() {
  >&2 echo -e "Error: $1. Exiting.\n"
  exit 1
}

vcf_check() {
  if [ ! -f $1 ]; then
    raise_error "$1 does not exist"
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    raise_error "$1 may be truncated"
  fi
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}


SECONDS=0

ln -f -s ${OUT} ${OUT}-tmp.vcf.gz
if [ ! -f ${OUT}.tbi ]; then
  tabix -p vcf ${OUT}-tmp.vcf.gz
fi

duration=${SECONDS}
echo "make_tabix finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}, $( date )"
