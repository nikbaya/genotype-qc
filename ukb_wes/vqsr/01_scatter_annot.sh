#!/usr/bin/env bash
#
# Author: Nik Baya (2021-01-19)
#
#$ -N test_scatter_annot
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/test_scatter_annot.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/test_scatter_annot.errors.log
#$ -q test.qc
#$ -V
#$ -P lindgren.prjc
#$ -t 13

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test"
readonly OUT="${WD}/test_scatter_annot_chr${CHR}" # output directory
readonly ANNOT_CHUNK_SCRIPT="${WD}/_test_run_annot_chunk.sh"

SECONDS=0

# split unique base pair positons of variants into equal chunks
readonly BIM="/well/ukbb-wes/calls/oqfe/ukbb-wes-oqfe-calls-chr${CHR}.bim"
readonly MAX_CHUNK_SIZE=10000 # maximum number of base pair positions per chunk
readonly SPLIT_VARIANTS_DIR="${OUT}/split_variants" # directory to place split files containing variants
readonly SPLIT_PREFIX="${SPLIT_VARIANTS_DIR}/split_variants_chr${CHR}_" # prefix for split command output files

mkdir -p ${SPLIT_VARIANTS_DIR}  # -p flag will not throw an error if directory already exists
split -l ${MAX_CHUNK_SIZE} <( cut -f1,4 ${BIM} | uniq ) ${SPLIT_PREFIX}

# create intervals
readonly INTERVALS="${OUT}/intervals_chr${CHR}.txt"

rm ${INTERVALS} 2> /dev/null # remove file if it exists
# iterate through the list of split files
while read split_file; do
  interval=$( paste <( head -1 $split_file ) <( tail -n1 $split_file | cut -f2 ) | awk '{ print "chr"$1":"$2"-"$3 }' )
  echo ${interval} >> ${INTERVALS}
done < <( ls -1 ${SPLIT_PREFIX}* )

N_CHUNKS=$( cat ${INTERVALS} | wc -l )

qsub -t 1:${N_CHUNKS} ${ANNOT_CHUNK_SCRIPT} ${CHR}

duration=${SECONDS}
echo "chr${CHR} scatter_annot finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
