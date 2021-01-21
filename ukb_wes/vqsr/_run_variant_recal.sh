#!/usr/bin/env bash
#
# This is an "internal" script for running GATK's VariantRecalibrator tool
# This script is called by 03_run_vqsr.sh
#
# Author: Nik Baya (2021-01-21)
#
#$ -N _variant_recal
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/variant_recal.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/variant_recal.errors.log
#$ -q short.qf
#$ -pe shmem 1
#$ -V
#$ -P lindgren.prjc

readonly IN=${1?Error: _run_variant_recal.sh requires the VCF to run VQSR on as 1st arg}
readonly VARIANT_TYPE=${2?Error: _run_variant_recal.sh requires the variant type (snp or indel) as 2nd arg}
readonly RECAL_PATH=${3?Error: _run_variant_recal.sh requires the output path prefix as 3rd arg} # output path
readonly MAX_GAUSS=${4?Error: _run_variant_recal.sh requires the max gaussian number as 4th arg} # expected number of clusters
readonly MEM_RECAL=${5?Error: _run_variant_recal.sh requires the memory allocated to the GATK JVM as the 5th arg}

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

if [ "${VARIANT_TYPE}" == "snp" ]; then
  gatk --java-options "-Xmx${MEM_RECAL}g -Xms${MEM_RECAL}g -XX:-UseParallelGC" VariantRecalibrator \
     -V ${IN} \
     --trust-all-polymorphic \
     -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
     -an QD -an FS -an SOR \
     -mode SNP \
     --max-gaussians ${MAX_GAUSS} \
     --rscript-file ${RECAL_PATH}_plots.R \
     -resource:hapmap,known=false,training=true,truth=true,prior=15 ${REF}/hapmap_3.3.hg38.vcf.gz \
     -resource:omni,known=false,training=true,truth=true,prior=12 ${REF}/1000G_omni2.5.hg38.vcf.gz \
     -resource:1000G,known=false,training=true,truth=false,prior=10 ${REF}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${REF}/Homo_sapiens_assembly38.dbsnp138.vcf \
     -O ${RECAL_PATH}.recal \
     --tranches-file ${RECAL_PATH}.tranches
  readonly EXIT_CODE=$?

elif [ $"${VARIANT_TYPE}" == "indel" ]; then
  gatk --java-options "-Xmx${MEM}g -Xms${MEM}g  -XX:-UseParallelGC" VariantRecalibrator \
      -V ${IN} \
      --trust-all-polymorphic \
      -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
      -an FS -an QD -an SOR \
      -mode INDEL \
      --max-gaussians ${MAX_GAUSS} \
      --rscript-file ${RECAL_PATH}_plots.R \
      -resource:mills,known=false,training=true,truth=true,prior=12 ${REF}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${REF}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${REF}/Homo_sapiens_assembly38.dbsnp138.vcf \
      -O ${RECAL_PATH}.recal \
      --tranches-file ${RECAL_PATH}.tranches
  readonly EXIT_CODE=$?

else
  raise_error "Invalid VARIANT_TYPE: ${VARIANT_TYPE}\nVARIANT_TYPE must be either \"snp\" or \"indel\""
fi

if [ $( ls -1 ${RECAL_PATH}.{recal,tranches,recal.idx} | wc -l ) -ne 3 ]; then
  raise_error "Not all ${RECAL_PATH}.{recal,tranches,recal.idx} files were completely written (GATK exit code: ${EXIT_CODE})"
fi

duration=${SECONDS}
echo "finished ${VARIANT_TYPE} variant recal, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
