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
#$ -l h_rt=00:05:00
#$ -pe shmem 1
#$ -V
#$ -P lindgren.prjc

# input
readonly BFILE=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/plink/ukb_wes_200k_pre_qc

# output
readonly WD=/well/lindgren/UKBIOBANK/nbaya/wes_200k/genotype_qc/plink
readonly TMP=${WD}/tmp-ukb_wes_200k
readonly OUT=${WD}/ukb_wes_200k

# filter parameters
readonly GENO1=0.05 # include only variants with missing-rate < $GENO1 (before sample filter), important for post merge of multiple platforms
readonly MIND=0.02     # include only samples with missing-rate < $MIND
readonly FHET=0.2   # include only samples with -$FHET < FHET < $FHET
readonly GENO2=0.02     # include only variants with missing-rate < $GENO2
# HWE_TH=1e-10  # include only variants with HWE p-val > $HWE_TH
# MAF=0.005     # include only variants with minor allele frequency < $MAF

for f in ${BFILE}.{bed,bim,fam}; do
  test ! -f $f && echo "Error: $f does not exist. Exiting" >&2 && exit 1
done

echo -e "\n###########
starting genotype qc (job id: ${JOB_ID}.${SGE_TASK_ID}, $(date))

BFILE: ${BFILE}
TMP:   ${TMP}
OUT:   ${OUT}

GENO1: ${GENO1}
MIND:  ${MIND}
FHET:  ${FHET}
GENO2: ${GENO2}
###########"

## SNP call rate 1st pass filtering
echo "...SNP call rate 1st pass filtering ($GENO1)..."
plink --bfile ${BFILE} \
  --geno ${GENO1} \
  --make-bed \
  --out ${TMP}_geno1

## Sample call rate filtering
echo "...Sample call rate filtering ($MIND)..."
plink --bfile ${TMP}_geno1 \
  --mind ${MIND} \
  --make-bed \
  --out ${TMP}_mind

## Inbreeding coefficient (Fhet)
echo "...Sample inbreeding coefficient filtering ($FHET)..."
plink --bfile ${TMP}_mind \
  --het \
  --out ${TMP}
awk -v FHET=${FHET} '{ if ($6 < -FHET || $6> FHET) print $1, $2, $6 }' ${TMP}.het > ${TMP}.remove.fhet.txt
plink --bfile ${TMP}_mind \
  --remove ${TMP}.remove.fhet.txt \
  --make-bed \
  --out ${TMP}_fhet

## SNP call rate 2nd pass filtering
echo "...SNP call rate 2nd pass filtering ($GENO2)..."
plink --bfile ${TMP}_fhet \
  --geno ${GENO2} \
  --make-bed \
  --out ${OUT}
