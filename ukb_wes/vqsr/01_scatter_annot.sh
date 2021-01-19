#!/usr/bin/env bash
#
# This script scatters the jobs of annotating chunks for each chromosome
#
# Author: Nik Baya (2021-01-19)
#
#$ -N scatter_annot
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/scatter_annot.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/scripts/scatter_annot.errors.log
#$ -q test.qc
#$ -V
#$ -P lindgren.prjc
#$ -t 9-12

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr"
readonly OUT="${WD}/vcf/scatter_annot_chr${CHR}" # output directory
readonly ANNOT_CHUNK_SCRIPT="${WD}/scripts/_run_annot_chunk.sh"

# split unique base pair positons of variants into equal chunks
readonly BIM="/well/ukbb-wes/calls/oqfe/ukbb-wes-oqfe-calls-chr${CHR}.bim"
readonly MAX_CHUNK_SIZE=10000 # maximum number of base pair positions per chunk
readonly SPLIT_VARIANTS_DIR="${OUT}/split_variants" # directory to place split files containing variants
readonly SPLIT_PREFIX="${SPLIT_VARIANTS_DIR}/split_variants_chr${CHR}_" # prefix for split command output files

mkdir -p ${SPLIT_VARIANTS_DIR}  # -p flag will not throw an error if directory already exists
split -l ${MAX_CHUNK_SIZE} <( cut -f1,4 ${BIM} | uniq ) ${SPLIT_PREFIX}

# create intervals
readonly INTERVALS="${OUT}/intervals_chr${CHR}.txt"

if [ ! -f ${INTERVALS} ]; then
  # iterate through the list of split files
  while read split_file; do
    interval=$( paste <( head -1 $split_file ) <( tail -n1 $split_file | cut -f2 ) | awk '{ print "chr"$1":"$2"-"$3 }' )
    echo ${interval} >> ${INTERVALS}
  done < <( ls -1 ${SPLIT_PREFIX}* )

  if [[  "${CHR}" = "X" || "${CHR}" = "Y" ]]; then
    sed -i "s/chr${SGE_TASK_ID}/chr${CHR}/g" ${INTERVALS} # note that SGE_TASK_ID is 23 or 24 for chr X and Y, respectively
  fi
fi

readonly QUEUE="short.qe" # queue to use for scattered annotation (default: short.qf)
readonly N_CORES=1 # number of cores (or "slots") to use (default: 1)
readonly N_CHUNKS=$( cat ${INTERVALS} | wc -l ) # number of chunks (i.e. intervals) for the given chromosome
readonly MEM=10 # memory in gb used for gat (default: 4 for 1 qc core, 8 for 2 qf cores)

qsub -q ${QUEUE} -pe shmem ${N_CORES} -t 1:${N_CHUNKS} -N "_c${CHR}_annot_chunk" ${ANNOT_CHUNK_SCRIPT} ${CHR} ${OUT} ${INTERVALS} ${MEM}
