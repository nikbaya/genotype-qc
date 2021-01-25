#!/usr/bin/env bash
#
# This is an "internal" script for running GATK annotation on genomic interval "chunks"
#
# Author: Nik Baya (2021-01-16)
#
#$ -N _annot_chunk
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/annot_chunk.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/annot_chunk.errors.log
#$ -V
#$ -P lindgren.prjc

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

readonly CHUNK_IDX=${SGE_TASK_ID}
readonly CHR=${1?Error: _run_annot_chunk.sh requires the chromosome as 1st arg} # take CHR as argumen
readonly IN=${2?Error: _run_annot_chunk.sh requires the input VCF as 2nd arg} # input VCF to annotate
readonly OUT_PREFIX=${3?Error: _run_annot_chunk.sh requires the output path prefix as 3rd arg} # output path prefix, e.g. /path/to/directory/ukb_wes_chr12
readonly INTERVALS=${4?Error: _run_annot_chunk.sh requires the intervals file or range as 4th arg} # intervals file or range, e.g. chr1:43209-5899930 (if a range, the arg must being with the string "chr")
readonly MEM=${5?Error: _run_annot_chunk.sh requires the GATK memory limit as 5th arg} # memory in gb used for gat

# NOTE: In the following .ped file parents have different FIDs than offspring, but since GATK only looks at whether
# the PAT/MAT fields are both nonzero when determing founder status, the founders will still be identified correctly
readonly PED="/gpfs3/well/lindgren/UKBIOBANK/nbaya/resources/ukb12788_simplified_gatk_pedigree.ped" # pedigree file


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


if [[ "$INTERVALS" == "chr"* && ! -f ${INTERVALS} ]]; then # if INTERVALS starts with "chr" and is not a file (redundant check)
  readonly OUT="${OUT_PREFIX}.vcf.gz" # output VCF
elif [ -f ${INTERVALS} ]; then
  readonly OUT="${OUT_PREFIX}.${CHUNK_IDX}of${SGE_TASK_LAST}.vcf.gz" # output
else
  raise_error "INTERVALS variable (4th input arg) must either be 1) an intervals file with rows of the form chr1:10000393-2000003382 or 2) a string of the same form"
fi

vcf_check ${IN}

SECONDS=0

if [ -f ${INTERVALS} ]; then
  readonly INTERVAL=$( sed "${CHUNK_IDX}q;d"  ${INTERVALS} )
else
  readonly INTERVAL=${INTERVALS}
fi

echo "starting ${INTERVAL} (${CHUNK_IDX}/${SGE_TASK_LAST}) (GATK mem: ${MEM}g, job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"

readonly OUT_FAILED="${OUT}-failed" # path for failed output VCF

if [ -f ${OUT_FAILED} ]; then
  echo "Removing existing failed file: ${OUT_FAILED}"
  rm ${OUT_FAILED}
fi

if [ ! -f ${OUT} ]; then

  set -x
  # NOTE: 2nd row of annotations are annotations that will not be calculated in the 200k cohort because of missing fields
  gatk --java-options "-Xmx${MEM}g -Xms${MEM}g -XX:-UseParallelGC" \
    VariantAnnotator \
    -V ${IN} \
    -L ${INTERVAL} \
    -A ExcessHet -A InbreedingCoeff -A StrandOddsRatio -A QualByDepth -A FisherStrand \
    -A ReadPosRankSumTest -A MappingQualityRankSumTest -A Coverage -A RMSMappingQuality \
    -ped ${PED} \
    --sites-only-vcf-output \
    -O ${OUT}
  set +x

  if [ $( bcftools view -h ${OUT} 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    mv ${OUT} ${OUT_FAILED}
    raise_error "GATK VariantAnnotator for chr${CHR} (${CHUNK_IDX}/${SGE_TASK_LAST}) did not successfully write output VCF, failed output VCF has been tagged with suffix \"-failed\" (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
  fi
fi

vcf_check ${OUT}

duration=${SECONDS}
echo "${INTERVAL} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
