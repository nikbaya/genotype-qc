#!/usr/bin/env bash
#
# This script merges the per-chromosome pVCFs
#
# Author: Nik Baya (2021-01-18)
#
#$ -N merge_vcfs
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/merge_vcfs.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/merge_vcfs.errors.log
#$ -q long.qf
#$ -pe shmem 40
#$ -V
#$ -P lindgren.prjc

# input and output
readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/vcf"
readonly OUT="${WD}/ukb_wes_oqfe_pvcf.vcf.gz"

readonly FILE_LIST="${WD}/merge_file_list.txt"

raise_error() {
  >&2 echo -e "Error: $1. Exiting.\n"
  exit 1
}

vcf_check() {
  if [ ! -f $1 ]; then
    raise_error "$1 does not exist"
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    raise_error"$1 may be truncated"
    exit 1
  fi
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}


SECONDS=0

if [ ! -f ${FILE_LIST} ]; then
  echo "/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr1.vcf.gz" > ${FILE_LIST}
  for chr in {2..24}; do
    if [ ${chr} -eq 23 ]; then
      chr="X"
    elif [ ${chr} -eq 24 ]; then
      chr="Y"
    fi
    echo "/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr${chr}.vcf.gz" >> ${FILE_LIST}
  done
fi

if [ ! -f ${OUT} ]; then
  # options invoked
  # -f, --file-list : file containing names of VCFs to merge
  # -n, --naive : Concatenates VCF files without recompression. Requires all VCFs to have same headers
  # -Oz : write output as compressed VCF
  # --threads : number of threads for multithreading
  # -o : output file path
  bcftools concat \
    -f ${FILE_LIST} \
    --naive \
    -Oz \
    --threads 40 \
    -o ${OUT}
fi

if [ ! -f ${OUT}.tbi ]; then
  tabix -p vcf ${OUT}
fi

duration=${SECONDS}
echo "merge_vcfs finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}, $( date )"
