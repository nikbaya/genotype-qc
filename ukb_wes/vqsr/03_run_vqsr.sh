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

# options: ukb_wes_200k, haplo_50, ukb_wes_150k, ukb_wes_50k
readonly SUBSET="ukb_wes_200k"
#readonly SUBSET="haplo_50"
#readonly SUBSET="haplo_50"

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/vcf/${SUBSET}"

# scripts being called at the end
readonly VARIANT_RECAL_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_variant_recal.sh"
readonly APPLY_VQSR_SCRIPT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/_run_apply_vqsr.sh"

# paths for merge step
readonly SCATTER_DIR_PREFIX="${WD}/scatter_annot_chr" # prefix to directories containing scatter annotation sites-only VCFs
readonly PER_CHROM_MERGE_PREFIX="${WD}/tmp-merge_chr"
#readonly PER_CHROM_MERGE_PREFIX="${WD}/annot/haplo_50_gvcf_chr" # prefix of path for intermediate merged chromosome VCF, default: ${WD}/tmp-merge_chr (can also be input paths for sites-only VCFs, if scatter annotation wasn't used; in this case, SCATTER_DIR_PREFIX is ignored by the script)
readonly MERGED="${WD}/${SUBSET}_sitesonly.vcf.gz" # name of merged sites-only VCF to create

# set of annotations to use during VariantRecalibrator
readonly ANNOT_SET="limited" # options: full or limited (3 annotations)

# input (pre-filter) and output (post-filter) VCF paths
# NOTE: The input VCF to be filtered does not need to be the same VCF upon which the VQSR runs
readonly FILTER_LEVEL="99.7" # filter level for ApplyVQSR step (default: 99.7)
readonly IN="${WD}/${SUBSET}_sitesonly.vcf.gz"  # VCF to be filtered using VQSR output
readonly OUT="${WD}/${SUBSET}_sitesonly_recal_annot${ANNOT_SET}_filter${FILTER_LEVEL}.vcf.gz" # output path for VCF after filter (ApplyVQSR step)
#readonly IN="${WD}/${SUBSET}_chr2_sitesonly.vcf.gz"
#readonly OUT="${WD}/${SUBSET}_chr2_sitesonly_recal_fl${FILTER_LEVEL}.vcf.gz"
#readonly IN=${MERGED}
#readonly OUT="${WD}/${SUBSET}_sitesonly_recal_annot${ANNOT_SET}_filter${FILTER_LEVEL}.vcf.gz"
#readonly IN="${WD}/${SUBSET}_sitesonly_chr10.vcf.gz"
#readonly OUT="${WD}/${SUBSET}_sitesonly_chr10_recal_annot${ANNOT_SET}_filter${FILTER_LEVEL}.vcf.gz"
#readonly IN="/well/lindgren/UKBIOBANK/saskia/haplo_50/ggvcf_chr21.vcf.gz"
#readonly OUT="${WD}/${SUBSET}_sitesonly_chr21_recal_annot${ANNOT_SET}_filter${FILTER_LEVEL}.vcf.gz"

# maximum gaussian args (expected number of clusters in data)
readonly MAX_GAUSS_SNP=6 # default: 6
readonly MAX_GAUSS_INDEL=4 # default: 4

# file prefixes for output of VariantRecalibrator
readonly RECAL_SNP="${WD}/${SUBSET}_snp_annot${ANNOT_SET}_maxgauss${MAX_GAUSS_SNP}"
readonly RECAL_INDEL="${WD}/${SUBSET}_indel_annot${ANNOT_SET}_maxgauss${MAX_GAUSS_INDEL}"

# queues to be used by jobs submitted by this script
readonly QUEUE_RECAL="short.qe"
readonly QUEUE_APPLY="short.qf"

# number of cores to be used by jobs submitted by this script
readonly N_CORES_RECAL=2
readonly N_CORES_APPLY=10

# memory in GB to be used by GATK as Java limits during various parts of the pipeline
# NOTE: MEM_RECAL and MEM_APPLY should depend on what queues and numbers of cores are specified above
readonly MEM_EXCESSHET=3 # this memory should be determined by the memory allocated to this script (3g per qf slot, 10g per qe slot)
readonly MEM_RECAL=20
readonly MEM_APPLY=35

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
  for chr_idx in {1..24}; do
    chr=$( get_chr_str ${chr_idx} )
    if [ ! -f ${PER_CHROM_MERGE_PREFIX}${chr}.vcf.gz ]; then
      bcftools concat \
        --file-list <( ls -1 ${SCATTER_DIR_PREFIX}${chr}/*vcf.gz | sort -V ) \
        --allow-overlaps \
        --rm-dups all \
        -Oz \
        -o ${PER_CHROM_MERGE_PREFIX}${chr}.vcf.gz
    else
      echo "Warning: ${PER_CHROM_MERGE_PREFIX}${chr}.vcf.gz already exists, skipping merge of genomic interval VCFs"
    fi
    vcf_check ${PER_CHROM_MERGE_PREFIX}${chr}.vcf.gz
  done

  # second stage
  bcftools concat \
    --file-list <( ls -1 ${PER_CHROM_MERGE_PREFIX}*.vcf.gz ) \
    --naive \
    -Oz \
    -o ${MERGED} \

  if [[ "${PER_CHROM_MERGE_PREFIX}" == *"tmp-"* ]]; then
    echo "Removing intermediate VCFs with paths matching ${PER_CHROM_MERGE_PREFIX}*" \
    && rm ${PER_CHROM_MERGE_PREFIX}*.vcf.gz
  fi
else
  time_check "Warning: ${MERGED} already exists, skipping merge step"
fi

vcf_check ${MERGED}
tabix_check ${MERGED}

# filter by excess heterozygosity
readonly TMP_EXCESSHET="${WD}/tmp-${SUBSET}_sitesonly_excesshet.vcf.gz" # intermediate VCF filtered by excess heterozygosity
readonly EXCESSHET_MAX=54.69 # maximum ExcessHet value allowed, any variants with ExcessHet greater than this threshold will be filtered, i.e. removed (default: 54.69,  NOTE: this default comes from the GATK default pipeline and correspnods to a z-score of -4.5)

if [ ! -f ${TMP_EXCESSHET} ]; then
  time_check "Starting ExcessHet filter for ${SUBSET} cohort (remove variants with ExcessHet>${EXCESSHET_MAX})"
  gatk --java-options "-Xmx${MEM_EXCESSHET}g -Xms${MEM_EXCESSHET}g -XX:-UseParallelGC" VariantFiltration \
    -V ${MERGED} \
    --filter-expression "ExcessHet > ${EXCESSHET_MAX}" \
    --filter-name ExcessHet \
    -O ${TMP_EXCESSHET} \
  || raise_error "GATK ExcessHet filter on ${MERGED} failed"
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
  local JOB_NAME="_${SUBSET}_${VARIANT_TYPE}_recal_${ANNOT_SET}"
  if [ $( ls -1 ${RECAL_PATH}.{recal,tranches,recal.idx} 2> /dev/null | wc -l ) -ne 3 ]; then
    qsub -N ${JOB_NAME} \
      -terse \
      -q ${QUEUE_RECAL} \
      -pe shmem ${N_CORES_RECAL} \
      ${VARIANT_RECAL_SCRIPT} \
      ${TMP_EXCESSHET} \
      ${VARIANT_TYPE} \
      ${ANNOT_SET} \
      ${RECAL_PATH} \
      ${MAX_GAUSS} \
      ${MEM_RECAL} > /dev/null
  fi
  echo ${JOB_NAME}
}


JOB_NAME_SNP=$( submit_recal "snp" ${RECAL_SNP} ${MAX_GAUSS_SNP} )
JOB_NAME_INDEL=$( submit_recal "indel" ${RECAL_INDEL} ${MAX_GAUSS_INDEL} )


# submit ApplyVQSR job

if [[ ${IN} == *"chr" && ${OUT} == *"chr" ]]; then # if both paths end with "chr"
  set -x
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
    ${FILTER_LEVEL} \
    ${MEM_APPLY}
  set +x
else
  vcf_check ${IN}
  set -x
  qsub -N "_${SUBSET}_apply_vqsr" \
    -hold_jid ${JOB_NAME_SNP},${JOB_NAME_INDEL} \
    -q ${QUEUE_APPLY} \
    -pe shmem ${N_CORES_APPLY} \
    ${APPLY_VQSR_SCRIPT} \
    ${IN} \
    ${OUT} \
    ${RECAL_SNP} \
    ${RECAL_INDEL} \
    ${FILTER_LEVEL} \
    ${MEM_APPLY}
  set +x
fi

duration=${SECONDS}
echo "finished submitinng all VQSR jobs for ${SUBSET} cohort, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
