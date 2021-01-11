#!/usr/bin/env bash
#
# This script is for running simple genotype QC
#
# Author: Nik Baya (2021-01-05)
#
#$ -N genotype_qc
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/genotype_qc/scripts/genotype_qc.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/genotype_qc/scripts/genotype_qc.errors.log
#$ -q short.qf
##$ -l h_rt=:00:00
#$ -pe shmem 40
#$ -V
#$ -P lindgren.prjc

readonly OVERWRITE=0 # 0 = false, 1 = true

# input
readonly BFILE="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink/ukb_wes_200k_pre_qc"

# output
readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/genotype_qc/plink"
readonly TMP="${WD}/tmp-ukb_wes_200k"
readonly OUT="${WD}/ukb_wes_200k"

# filter parameters
# NOTE: We choose less stringent filters than the default to only remove the truly bad outliers
readonly GENO1=0.10 # include only variants with missing-rate < $GENO1 (before sample filter), important for post merge of multiple platforms (default: 0.05)
readonly MIND=0.10     # include only samples with missing-rate < $MIND (default: 0.02)
#readonly FHET=0.2   # include only samples with -$FHET < FHET < $FHET (default: 0.2)
readonly GENO2=0.10     # include only variants with missing-rate < $GENO2 (default: 0.02)
# HWE_TH=1e-10  # include only variants with HWE p-val > $HWE_TH
# MAF=0.005     # include only variants with minor allele frequency < $MAF

time_check() {
  echo -e "\n########\n$1 (job id: ${JOB_ID}.${SGE_TASK_ID}, $(date))\n########"
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

plink_files_exist() {
  [ ` ls -1 $1.{bed,bim,fam} 2> /dev/null | wc -l ` -eq 3 ]
}

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

check_plink_files() {
  if ! plink_files_exist $1; then
    raise_error "$1.{bed,bim,fam} files do not all exist."
  fi
}

if [ ${OVERWRITE} -ne 0 ] && [ ${OVERWRITE} -ne 1 ]; then
  echo -e "Error: OVERWRITE must be either 0 (do not overwrite) or 1 (overwrite). Exiting."
  exit 1
fi

check_plink_files ${BFILE}

time_check "Starting genotype qc"
echo -e "\n###########
BFILE: ${BFILE}
TMP:   ${TMP}
OUT:   ${OUT}

GENO1: ${GENO1}
MIND:  ${MIND}
GENO2: ${GENO2}
###########"

SECONDS=0

## SNP call rate 1st pass filtering
echo "Starting SNP call rate 1st pass filtering ($GENO1)"
plink --bfile ${BFILE} \
  --geno ${GENO1} \
  --make-bed \
  --out ${TMP}_geno1
check_plink_files ${TMP}_geno1

## Sample call rate filtering
echo "Starting sample call rate filtering ($MIND)"
plink --bfile ${TMP}_geno1 \
  --mind ${MIND} \
  --make-bed \
  --out ${TMP}_mind
check_plink_files ${TMP}_mind

## Inbreeding coefficient (Fhet)
#echo "...Sample inbreeding coefficient calculation ($FHET)..."
#plink --bfile ${TMP}_mind \
#  --het \
#  --out ${TMP}
#awk -v FHET=${FHET} '{ if ($6 < -FHET || $6> FHET) print $1, $2, $6 }' ${TMP}.het > ${TMP}.remove.fhet.txt
#plink --bfile ${TMP}_mind \
#  --remove ${TMP}.remove.fhet.txt \
#  --make-bed \
#  --out ${TMP}_fhet

## SNP call rate 2nd pass filtering
echo "Starting SNP call rate 2nd pass filtering ($GENO2)"
plink --bfile ${TMP}_mind \
  --geno ${GENO2} \
  --make-bed \
  --out ${OUT}
check_plink_files ${OUT}

duration=${SECONDS}
echo "genotype qc finished, elapsed time: $( echo "scale=2; ${duration}/3600" | bc -l ) hrs (job id: ${JOB_ID}, $( date )"
