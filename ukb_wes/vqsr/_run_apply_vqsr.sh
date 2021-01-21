#!/usr/bin/env bash
#
# This is an "internal" script for running GATK's ApplyVQSR tool
# This script is called by 03_run_vqsr.sh
#
# Author: Nik Baya (2021-01-21)
#
#$ -N _apply_vqsr
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/apply_vqsr.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/apply_vqsr.errors.log
#$ -q short.qf
#$ -pe shmem 1
#$ -V
#$ -P lindgren.prjc

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

readonly IN=${1?Error: _run_apply_vqsr.sh requires the VCF to be filtered by ApplyVQSR as 1st arg}
readonly OUT=${2?Error: _run_apply_vqsr.sh requires the path of the output VCF as 2nd arg}
readonly RECAL_SNP=${3?Error: _run_apply_vqsr.sh requires the snp recal path prefix as 3rd arg}
readonly RECAL_INDEL=${4?Error: _run_apply_vqsr.sh requires the indel recal path prefix as 4th arg}
readonly MEM=${5?Error: _run_apply_vqsr.sh requires the memory allocated to the GATK JVM as the 5th arg}

readonly REF="/well/lindgren/UKBIOBANK/nbaya/resources/ref" # directory containing reference files

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

readonly TMP_APPLY_VCF="${RECAL_SNP}-tmp.vcf.gz" # intermediate VCF

if [ ! -f ${OUT} ]; then
  if [ ! -f ${TMP_APPLY_VCF} ]; then
    # filter input file using SNP recalibration results
    gatk --java-options "-Xmx${MEM}g -Xms${MEM}g -XX:-UseParallelGC" ApplyVQSR \
      -V ${IN} \
      --recal-file ${RECAL_SNP}.recal \
      --tranches-file ${RECAL_SNP}.tranches \
      --truth-sensitivity-filter-level 99.7 \
      --create-output-variant-index true \
      -mode SNP \
      -O ${TMP_APPLY_VCF}
  else
    echo "Warning: ${TMP_APPLY_VCF} already exists, skipping first stage of ApplyVQSR"
  fi

  vcf_check ${TMP_APPLY_VCF}

  # filter intermediate file using INDEL recalibration results
  gatk --java-options "-Xmx${MEM}g -Xms${MEM}g -XX:-UseParallelGC" ApplyVQSR \
    -V ${TMP_APPLY_VCF} \
    --recal-file ${RECAL_INDEL}.recal \
    --tranches-file ${RECAL_INDEL}.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -O ${OUT}
else
  echo "Warning: ${OUT} already exists, skipping ApplyVQSR"
fi

vcf_check ${OUT}

rm ${TMP_APPLY_VCF} && echo "Removed

duration=${SECONDS}
echo "finished ${VARIANT_TYPE} variant recal, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
