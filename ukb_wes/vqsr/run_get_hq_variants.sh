#!/usr/bin/env bash
#
# This script gets sets of high-quality (or high-confidence) variants
#
# Author: Nik Baya (2021-02-01)
#
#$ -N get_hq_variants
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/get_hq_variants.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/get_hq_variants.errors.log
#$ -q short.qf
#$ -V
#$ -P lindgren.prjc
#$ -t 1-24

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

if [[ "${QUEUE}" != "long."* ]] && [[ ${CHR} -eq 1 || ${CHR} -eq 2 || ${CHR} -eq 19 ]]; then
  readonly QUEUE="long.qf"
  set -x
  qsub -q ${QUEUE} \
    -N "${JOB_NAME}_resubmit" \
    -t ${SGE_TASK_ID} \
    -pe shmem ${NSLOTS} \
    /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/run_get_hq_variants.sh $1 \
  && echo "resubmitting chr${CHR} hq variant job using ${QUEUE}" \
  && exit 1
fi

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr"

readonly SUBSET="$1" # options: 200k, 150k, 50k

if [[ "${SUBSET}" == "200k" ]]; then
  # full 200k WES cohort
  readonly IN="/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr${CHR}.vcf.gz"
elif [[ "${SUBSET}" == "150k" || "${SUBSET}" == "50k" ]]; then
  # 150k or 50k subset cohorts
  readonly IN="/well/lindgren/UKBIOBANK/nbaya/resources/ukb_wes_${SUBSET}/ukb_wes_oqfe_pvcf_${SUBSET}_chr${CHR}.vcf.gz"
else
  >&2 echo "Error: SUBSET must be either 200k, 150k or 50k. Exiting" && exit 1
fi

readonly VERSION="1.2" # options: 1.1, 1.2
readonly OUT_DIR="${WD}/stats/ukb_wes_${SUBSET}" # output directory
readonly OUT="${OUT_DIR}/ukb_wes_oqfe_pvcf_${SUBSET}_sitesonly_hq_v${VERSION}_chr${CHR}.vcf.gz" # output VCF name

if [[ -z $IN || -z ${OUT_DIR} || -z ${OUT} ]]; then
  >&2 "Error: Variables IN, OUT_DIR, OUT_PREFIX must all be defined" && exit 1
fi

raise_error() {
  >&2 echo -e "Error: $1. Exiting." && exit 1
}

vcf_check() {
  if [ ! -f $1 ]; then
    raise_error "$1 does not exist"
  elif [ ! -s $1 ]; then
    raise_error "$1 exists but is empty"
  elif [ $( bcftools view -h $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    raise_error "$1 may be truncated"
  fi
}

make_tabix() {
  module load BCFtools/1.10.2-GCC-8.3.0
  tabix -f -p vcf $1 \
  || ( exit_code=$? \
  && raise_error "tabix of $1 did not exit successfully (exit code: ${exit_code}, job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))" )
  module unload BCFtools/1.10.2-GCC-8.3.0
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

SECONDS=0

if [ ! -f ${OUT} ]; then
  mkdir -p ${OUT_DIR}
  vcf_check ${IN}
  if [[ "${VERSION}" == "1.1" ]]; then
    readonly F_MISS=0.01
    readonly MIN_AF=0.001
    bcftools view ${IN} \
      -G \
      -m2 \
      -M2 \
      --types snps \
      -i "F_MISS<${F_MISS} & INFO/AF>${MIN_AF}" \
      --threads 10 \
      -o ${OUT}
  elif [[ "${VERSION}" == "1.2" ]]; then
    readonly F_MISS=0.01
    readonly MIN_AF=0.001
    readonly MAX_AF=0.999
    bcftools view ${IN} \
      -G \
      -m2 \
      -M2 \
      --types snps \
      -i "F_MISS<${F_MISS} & INFO/AF>${MIN_AF} & INFO/AF<${MAX_AF}" \
      --threads 10 \
      -o ${OUT}
  fi
fi

vcf_check ${OUT}

if [ ! -f ${OUT}.tbi ]; then
  make_tabix ${OUT}
fi

duration=${SECONDS}
echo "finished ${SUBSET} chr${CHR} v${VERSION} high-quality variant set creation, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
