#!/usr/bin/env bash
#
# This script runs all steps related to VQSR:
# 1) Merge scattered annotated genomic location sites-only VCFs into single sites-only VCF
# 2) Filter by Excesshet
# 3) Submit VariantRecalibrator jobs
# 4) Submit ApplyVQSR jobs
#
# The pipeline is based on suggestions from this page:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
#
# Author: Nik Baya (2021-01-20)
#
#$ -N vqsr
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/vqsr.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/vqsr.errors.log
#$ -q long.qf
#$ -pe shmem 10
#$ -V
#$ -P lindgren.prjc

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11
module load R/3.6.2-foss-2019b

readonly SUBSET="200k"

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/vcf/ukb_wes_${SUBSET}"

# scripts being called at the end
readonly VARIANT_RECAL_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_variant_recal.sh"
readonly APPLY_VQSR_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_apply_vqsr.sh"

# paths for merge step
readonly SCATTER_DIR_PREFIX="${WD}/scatter_annot_chr" # prefix to directories containing scatter annotation VCFs
readonly MERGED="${WD}/ukb_wes_${SUBSET}_sitesonly.vcf.gz" # name of merged sites-only VCF to create

# maximum gaussian args (expected number of clusters in data)
readonly MAX_GAUSS_SNP=6 # default: 6
readonly MAX_GAUSS_INDEL=4 # default: 4

# file prefixes for output of VariantRecalibrator
readonly RECAL_SNP="${WD}/ukb_wes_${SUBSET}_snp_maxgauss${MAX_GAUSS_SNP}"
readonly RECAL_INDEL="${WD}/ukb_wes_${SUBSET}_indel_maxgauss${MAX_GAUSS_INDEL}"

readonly MEM=30 # memory in GB used for GATK Java options

time_check() {
  echo -e "\n########\n$1 (job id: ${JOB_ID}.${SGE_TASK_ID}, $(date))\n########"
}

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

vcf_check() {
  if [ ! -f $1 ]; then
    raise_error "$1 does not exist"
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    raise_error "$1 may be truncated"
  fi
}

SECONDS=0



# merge genomic interval VCFs
if [ ! -f ${MERGED} ]; then
  time_check "Starting merge for ${SUBSET} cohort"
  bcftools concat \
    --file-list <( ls -1 ${SCATTER_DIR_PREFIX}*/*vcf.gz | sort -V ) \
    --naive \
    --Oz \
    --threads 20 \
    -o ${MERGED}
else
  time_check "Warning: ${MERGED} already exists, skipping merge step"
fi

vcf_check ${MERGED}

readonly TMP_EXCESSHET="${WD}/tmp-ukb_wes_${SUBSET}_sitesonly_excesshet.vcf.gz" # intermediate VCF filtered by excess heterozygosity



# filter by excess heterozygosity
readonly EXCESSHET_MAX=54.69 # maximum ExcessHet value allowed, any variants with ExcessHet greater than this threshold will be filtered, i.e. removed (default: 54.69,  NOTE: this default comes from the GATK default pipeline and correspnods to a z-score of -4.5)

if [ ! -f ${TMP_EXCESSHET} ]; then
  time_check "Starting ExcessHet filter for ${SUBSET} cohort"
  gatk --java-options "-Xmx${MEM}g -Xms${MEM}g -XX:-UseParallelGC" VariantFiltration \
    -V ${MERGED} \
    --filter-expression "ExcessHet > ${EXCESS_HET}" \
    --filter-name ExcessHet \
    -O ${TMP_EXCESSHET}
else
  time_check "Warning: ${TMP_EXCESSHET} already exists, skipping removal of variants with ExcessHet>${EXCESSHET_MAX}"
fi

vcf_check ${TMP_EXCESSHET}



# submit VariantRecalibrator jobs

if [ $( ls -1 ${RECAL_SNP}.{recal,tranches,recal.idx ) -ne 3 ]; then
  time_check "Submitting SNP VariantRecalibrator job for ${SUBSET} cohort"
  readonly SNP_JOB_ID=$( qsub -terse
