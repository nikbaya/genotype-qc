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
#$ -q short.qf
#$ -pe shmem 1
#$ -V
#$ -P lindgren.prjc

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

readonly SUBSET="200k"

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/vcf/ukb_wes_${SUBSET}"

# scripts being called at the end
readonly VARIANT_RECAL_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_variant_recal.sh"
readonly APPLY_VQSR_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_apply_vqsr.sh"

# paths for merge step
readonly SCATTER_DIR_PREFIX="${WD}/scatter_annot_chr" # prefix to directories containing scatter annotation VCFs
readonly MERGED="${WD}/ukb_wes_${SUBSET}_sitesonly.vcf.gz" # name of merged sites-only VCF to create

# input (pre-filter) and output (post-filter) VCF paths
# NOTE: The input VCF to be filtered does not need to be the same VCF upon which the VQSR runs
readonly IN=$1 # VCF to be filtered using VQSR output
readonly OUT=$2 # output path for VCF after filter

# maximum gaussian args (expected number of clusters in data)
readonly MAX_GAUSS_SNP=6 # default: 6
readonly MAX_GAUSS_INDEL=4 # default: 4

# file prefixes for output of VariantRecalibrator
readonly RECAL_SNP="${WD}/ukb_wes_${SUBSET}_snp_maxgauss${MAX_GAUSS_SNP}"
readonly RECAL_INDEL="${WD}/ukb_wes_${SUBSET}_indel_maxgauss${MAX_GAUSS_INDEL}"

# queues to be used by jobs submitted by this script
readonly QUEUE_RECAL="short.qe"
readonly QUEUE_APPLY="short.qf"

# number of cores to be used by jobs submitted by this script
readonly N_CORES_RECAL=2
readonly N_CORES_APPLY=2

# memory in GB to be used by GATK as Java limits during various parts of the pipeline
# NOTE: MEM_RECAL and MEM_APPLY should depend on what queues and numbers of cores are specified above
readonly MEM_EXCESSHET=3 # this memory should be determined by the memory allocated to this script (3g per qf slot, 10g per qe slot)
readonly MEM_RECAL=20
readonly MEM_APPLY=6

time_check() {
  echo -e "\n########\n$1 (job id: ${JOB_ID}, $(date))\n########"
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

get_chr_str() {
  local CHR=$1
  if [ ${CHR} -eq 23 ]; then
    CHR="X"
  elif [ ${CHR} -eq 24 ]; then
    CHR="Y"
  fi
  echo ${CHR}
}

tabix_check() {
  local VCF=$1
  if [ ! -f ${VCF}.tbi ]; then
    tabix -p vcf ${VCF}
    if [[ $? -ne 0 ]]; then
      raise_error "tabix of ${VCF} did not successfully complete"
    fi
  else
    time_check "Warning: ${VCF}.tbi already exists, skipping tabix step"
  fi
}

SECONDS=0



# two-step merge genomic interval VCFs
if [ ! -f ${MERGED} ]; then
  time_check "Starting merge for ${SUBSET} cohort"
  # first step of merge: genomic intervals -> per-chrom VCFs
  readonly TMP_MERGE_PREFIX="${WD}/tmp-merge_chr" # prefix of path for intermediate merged chromosome VCF
  for chr_idx in {1..24}; do
    chr=$( get_chr_str ${chr_idx} )
    if [ ! -f ${TMP_MERGE_PREFIX}${chr}.vcf.gz ]; then
      bcftools concat \
        --file-list <( ls -1 ${SCATTER_DIR_PREFIX}${chr}/*vcf.gz | sort -V ) \
        --allow-overlaps \
        --rm-dups all \
        -Oz \
        -o ${TMP_MERGE_PREFIX}${chr}.vcf.gz
    fi
  done

  # second stage
  bcftools concat \
    --file-list <( ls -1 ${TMP_MERGE_PREFIX}*.vcf.gz ) \
    --naive \
    -Oz \
    -o ${MERGED} \
  && echo "Removing intermediate VCFs with paths matching ${TMP_MERGE_PREFIX}*" \
  && rm ${TMP_MERGE_PREFIX}*.vcf.gz
else
  time_check "Warning: ${MERGED} already exists, skipping merge step"
fi

vcf_check ${MERGED}
tabix_check ${MERGED}

# filter by excess heterozygosity
readonly TMP_EXCESSHET="${WD}/tmp-ukb_wes_${SUBSET}_sitesonly_excesshet.vcf.gz" # intermediate VCF filtered by excess heterozygosity
readonly EXCESSHET_MAX=54.69 # maximum ExcessHet value allowed, any variants with ExcessHet greater than this threshold will be filtered, i.e. removed (default: 54.69,  NOTE: this default comes from the GATK default pipeline and correspnods to a z-score of -4.5)

if [ ! -f ${TMP_EXCESSHET} ]; then
  time_check "Starting ExcessHet filter for ${SUBSET} cohort (remove variants with ExcessHet>${EXCESSHET_MAX})"
  gatk --java-options "-Xmx${MEM_EXCESSHET}g -Xms${MEM_EXCESSHET}g -XX:-UseParallelGC" VariantFiltration \
    -V ${MERGED} \
    --filter-expression "ExcessHet > ${EXCESSHET_MAX}" \
    --filter-name ExcessHet \
    -O ${TMP_EXCESSHET}
else
  time_check "Warning: ${TMP_EXCESSHET} already exists, skipping removal of variants with ExcessHet>${EXCESSHET_MAX}"
fi

vcf_check ${TMP_EXCESSHET}
tabix_check ${TMP_EXCESSHET}

# submit VariantRecalibrator jobs

submit_recal() {
  local VARIANT_TYPE=$1 # "snp" or "indel"
  local RECAL_PATH=$2 # prefix for recal output
  local MAX_GAUSS=$3 # max gaussian argument
  if [ $( ls -1 ${RECAL_PATH}.{recal,tranches,recal.idx} 2> /dev/null | wc -l ) -ne 3 ]; then
    local JOB_NAME="_${SUBSET}_${VARIANT_TYPE}_recal"
    qsub -N ${JOB_NAME} \
      -q ${QUEUE_RECAL} \
      -pe shmem ${N_CORES_RECAL} \
      ${VARIANT_RECAL_SCRIPT} \
      ${TMP_EXCESSHET} \
      ${VARIANT_TYPE} \
      ${RECAL_PATH} \
      ${MAX_GAUSS} \
      ${MEM_RECAL}
  fi
  echo ${JOB_NAME}
}


JOB_NAME_SNP=$( submit_recal "snp" ${RECAL_SNP} ${MAX_GAUSS_SNP} )
JOB_NAME_INDEL=$( submit_recal "indel" ${RECAL_INDEL} ${MAX_GAUSS_INDEL} )


# submit ApplyVQSR job

if [[ ! -z "${IN}" && ! -z "${OUT}" ]]; then
  if [[ ${IN} != *".vcf.gz" && ${OUT} != *".vcf.gz"  ]]; then # if both paths don't end with ".vcf.gz"
    if [[ ${IN} == *"chr" && ${OUT} == *"chr" ]]; then # if both paths end with "chr"
      qsub -N "_${SUBSET}_apply_vqsr" \
        -t 1:24 \
        -hold_jid ${JOB_NAME_SNP},${JOB_NAME_INDEL} \
        -q ${QUEUE_APPLY} \
        -pe shmem ${N_CORES_APPLY} \
        ${APPLY_VQSR_SCRIPT} \
        ${IN} \
        ${OUT} \
        ${RECAL_SNP} \
        ${RECAL_INDEL} \
        ${MEM_APPLY}
    elif [[ ${IN} != *"chr" && ${OUT} != *"chr" ]]; then
      qsub -N "_${SUBSET}_apply_vqsr" \
        -hold_jid ${JOB_NAME_SNP},${JOB_NAME_INDEL} \
        -q ${QUEUE_APPLY} \
        -pe shmem ${N_CORES_APPLY} \
        ${APPLY_VQSR_SCRIPT} \
        ${IN} \
        ${OUT} \
        ${RECAL_SNP} \
        ${RECAL_INDEL} \
        ${MEM_APPLY}
    else
      raise_error "IN and OUT args must either 1) both end in \"chr\"f or 2) both not end in \"chr\""
    fi
  else
    raise_error "IN and OUT args cannot end with \".vcf.gz\", this suffix is added automatically"
  fi
else
  raise_error "Include IN and OUT file paths as 1st and 2nd args, respectively, in order to run ApplyVQSR step"
fi

duration=${SECONDS}
echo "finished submitinng all VQSR jobs for ${SUBSET} cohort, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
