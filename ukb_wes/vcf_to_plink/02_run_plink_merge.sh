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
#$ -l h_rt=100:00:00
#$ -pe shmem 20
#$ -V
#$ -P lindgren.prjc

readonly WD=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink

readonly BFILE=${WD}/ukb_wes_200k_chr
readonly MERGELIST=${WD}/mergelist.txt
readonly OUT=${WD}/ukb_wes_200k_pre_qc

plink_files_exist() {
  [ ` ls -1 $1.{bed,bim,fam} 2> /dev/null | wc -l ` -eq 3 ]
}

# create merge list file
if [ ! -f ${MERGELIST} ]; then
  for chr in {1..24}; do
    if [ ${chr} -eq 23 ]; then
      chr="X"
    elif [ ${chr} -eq 24 ]; then
      chr="Y"
    fi
    if plink_files_exist ${BFILE}${CHR}; then
      echo "${BFILE}${CHR}" >> ${MERGELIST}
    fi
  done
else
  echo "Warning: ${MERGELIST} already exists, skipping file creation"
  while read line; do
    if ! plink_files_exist ${line}; then
      echo -e "Error: ${line}.{bed,bim,fam} files do not exist.\nPlease double-check the contents of ${MERGELIST}"
      exit 1
    fi
  done < ${MERGELIST}
fi

echo -e "\nstarting plink merge (job id: ${JOB_ID}, $( date ))\n"
cat ${MERGELIST}

exit 0

if ! plink_files_exist ${OUT}; then
  plink --merge-list ${MERGELIST} \
    --make-bed \
    --memory 160000 \
    --out ${OUT}
  if ! plink_files_exist ${OUT}; then
    echo "Error: ${OUT}.{bed,bim,fam} was not successfully written."
    exit 1
  fi
else
  echo "Warning: ${OUT}.{bed,bim,fam} already exist, skipping PLINK merge"
fi
