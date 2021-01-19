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
readonly OUT=${2?Error: _run_annot_chunk.sh requires the output directory as 2nd arg} # output directory
readonly MEM=${3?Error: _run_annot_chunk.sh requires the GATK memory limit as 3rd arg} # memory in gb used for gat

readonly IN="/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr${CHR}.vcf.gz"
readonly OUT_CHUNK="${OUT}/ukb_wes_oqfe_pvcf_chr${CHR}.${CHUNK_IDX}of${SGE_TASK_LAST}.vcf.gz" # VCF of chunk output

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

vcf_check() {
  if [ ! -f $1 ]; then
    raise_error "$1 does not exist."
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    raise_error "$1 may be truncated"
  fi
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

vcf_check ${IN}

SECONDS=0

readonly INTERVALS="${OUT}/intervals_chr${CHR}.txt" # intervals file

if [ -f ${INTERVALS} ]; then
  readonly INTERVAL=$( sed "${CHUNK_IDX}q;d"  ${INTERVALS} )
else
  raise_error "${INTERVALS} file does not exist."
fi

echo "starting ${INTERVAL} (${CHUNK_IDX}/${SGE_TASK_LAST}) (GATK mem: ${MEM}g, job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"

if [ ! -f ${OUT_CHUNK} ]; then

  set -x
  gatk --java-options "-Xmx${MEM}g -Xms${MEM}g" VariantAnnotator \
    -V ${IN} \
    -L ${INTERVAL} \
    -A ExcessHet -A InbreedingCoeff -A StrandOddsRatio -A QualByDepth -A FisherStrand \
    -O ${OUT_CHUNK}

  readonly EXIT_CODE=$?
  set +x

  if [ ${EXIT_CODE} -ne 0 ]; then
    mv ${OUT_CHUNK} ${OUT_CHUNK}-failed
    raise_error "GATK VariantAnnotator had exit code ${EXIT_CODE}, failed output file has been tagged (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
  fi
fi

vcf_check ${OUT_CHUNK}

duration=${SECONDS}
echo "${INTERVAL} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
