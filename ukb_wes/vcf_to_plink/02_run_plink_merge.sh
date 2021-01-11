#!/usr/bin/env bash
#
# This script is for merging the per-chromosome PLINK files
#
# Author: Nik Baya (2021-01-05)
#
#$ -N plink_merge
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/scripts/plink_merge.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/scripts/plink_merge.errors.log
#$ -q short.qf
##$ -l h_rt=10:00:00
#$ -pe shmem 40
#$ -V
#$ -P lindgren.prjc

readonly WD=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink

readonly BFILE=${WD}/ukb_wes_200k_chr
readonly MERGELIST=${WD}/mergelist.txt
readonly OUT=${WD}/ukb_wes_200k_pre_qc

plink_files_exist() {
  [ ` ls -1 $1.{bed,bim,fam} 2> /dev/null | wc -l ` -eq 3 ]
}

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

# create merge list file
if [ ! -f ${MERGELIST} ]; then
  for chr in {1..24}; do
    if [ ${chr} -eq 23 ]; then
      chr="X"
    elif [ ${chr} -eq 24 ]; then
      chr="Y"
    fi
    if plink_files_exist ${BFILE}${chr}; then
      echo "${BFILE}${chr}" >> ${MERGELIST}
    else
      raise_error "${BFILE}${chr}.{bed,bim,fam} files are not successfully written."
    fi
  done
else
  echo "Warning: ${MERGELIST} already exists, skipping file creation"
  # if merge list exists, check that each entry is a valid bfile before proceeding
  while read line; do
    if ! plink_files_exist ${line}; then
      raise_error "${line}.{bed,bim,fam} files do not exist.\nPlease double-check the contents of ${MERGELIST}"
    fi
  done < ${MERGELIST}
fi

echo -e "\nstarting plink merge (job id: ${JOB_ID}, $( date ))\n"
echo -e "Files to merge:\n$( cat ${MERGELIST} )"

SECONDS=0

# create list of withdrawn samples
readonly WITHDRAWN="${WD}/withdrawn_samples.txt"
if [ ! -f ${WITHDRAWN} ]; then
  awk '{ print $1,$2 }' ${BFILE}1.fam | grep "-" > ${WITHDRAWN} # assume fam files are identical across chroms, so we can use the chr1.fam file
fi
echo "Removing withdrawn samples (n=$( cat ${WITHDRAWN} | wc -l))"

if ! plink_files_exist ${OUT}; then
  plink --merge-list ${MERGELIST} \
    --remove ${WITHDRAWN} \
    --make-bed \
    --memory 160000 \
    --out ${OUT}
  if ! plink_files_exist ${OUT}; then
    raise_error "${OUT}.{bed,bim,fam} was not successfully written."
  fi
else
  echo "Warning: ${OUT}.{bed,bim,fam} already exist, skipping PLINK merge"
fi

duration=${SECONDS}
echo "plink merge finished, elapsed time: $( echo "scale=2; ${duration}/3600" | bc -l ) hrs (job id: ${JOB_ID}, $( date )"
