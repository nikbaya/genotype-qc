#!/usr/bin/env bash
#
# This is an "internal" script for concatenating files
#
# Author: Nik Baya (2021-01-27)
#
#$ -N _concat
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/concat.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/concat.errors.log
#$ -q short.qe
#$ -V
#$ -P lindgren.prjc

readonly FILE_LIST=${1?Error: _run_concat.sh requires the files to be concatenated as 1st arg} # input VCFs to concatenate
readonly OUT=${2?Error: _run_concat.sh requires the output path as 2nd arg} # output path
#readonly TMP_DIR=${3?Error: _run_concat.sh requires the temporary directory as 3rd arg} # temporary directory

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

make_tabix() {
  tabix -f -p vcf $1 \
  || ( EXIT_CODE=$? \
  && raise_error "tabix did not exit properly (exit code ${EXIT_CODE}, job id: ${JOB_ID} $( date ))" )
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

while read f; do
  vcf_check $f
done < ${FILE_LIST}

SECONDS=0

if [ ! -f ${OUT} ]; then
  set -x
  bcftools concat \
    --file-list ${FILE_LIST} \
    --allow-overlaps \
    --rm-dups all \
    -Oz \
    -o ${OUT} \
  || ( EXIT_CODE=$? \
  && raise_error "concat did not exit properly (exit code: ${EXIT_CODE}, job id: ${JOB_ID} $( date ))" )

  echo "starting tabix for ${OUT} (job id: ${JOB_ID} $( date ))"
  make_tabix ${OUT}
fi

vcf_check ${OUT}

if [ ! -f ${OUT}.tbi ]; then
  make_tabix ${OUT}
fi

duration=${SECONDS}
echo "concat ${FILE_LIST} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
