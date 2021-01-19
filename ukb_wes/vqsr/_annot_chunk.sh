#!/usr/bin/env bash
#
# This is an "internal" script for running GATK annotation on genomic interval "chunks"
#
# Author: Nik Baya (2021-01-16)
#
#$ -N _test_annot_chunk
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/test_annot_chunk.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/test_annot_chunk.errors.log
#$ -q short.qf
#$ -pe shmem 2
#$ -V
#$ -P lindgren.prjc

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

readonly CHUNK_IDX=${SGE_TASK_ID}
readonly CHR=$1 # take CHR as argumen

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test"
readonly IN="${WD}/ukbb-wes-oqfe-pvcf-chr${CHR}.vcf.gz"
readonly OUT="${WD}/test_scatter_annot_chr${CHR}" # output directory
readonly OUT_CHUNK="${OUT}/ukb_wes_oqfe_pvcf_chr${CHR}.${CHUNK_IDX}of${SGE_TASK_LAST}.vcf.gz" # VCF of chunk output

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

vcf_check() {
  if [ ! -f $1 ]; then
    >&2 "Error: $1 does not exist. Exiting.\n"
    exit 1
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    >&2 "Error: $1 may be truncated. Exiting.\n"
    exit 1
  fi
}

vcf_check ${IN}

SECONDS=0

readonly INTERVALS="${OUT}/intervals_chr${CHR}.txt" # intervals file
readonly INTERVAL=$( sed "${CHUNK_IDX}q;d"  ${INTERVALS} )

echo "starting ${INTERVAL} (${CHUNK_IDX}/${SGE_TASK_LAST})"

readonly MEM=8 # memory in gb used for gatk

if [ ! -f ${OUT_CHUNK} ]; then
  gatk --java-options "-Xmx${MEM}g -Xms${MEM}g" VariantAnnotator \
    -V ${IN} \
    -L ${INTERVAL} \
    -A ExcessHet -A InbreedingCoeff -A StrandOddsRatio -A QualByDepth -A FisherStrand \
    -O ${OUT_CHUNK}
fi

vcf_check ${OUT_CHUNK}

duration=${SECONDS}
echo "${INTERVAL} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
