#!/usr/bin/env bash
#
# This script is for converting VCFs to PLINK.
#
# 1) Create temporary intermediate VCF
#    - Multiallelic sites are split
#    - Indels are left-aligned and normalised
# 2) VCF is converted to PLINK
# 3) Variant IDs are updated
#
# Author: Nik Baya (2021-01-01)
#
#$ -N vcf_to_plink
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/scripts/vcf_to_plink.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/scripts/vcf_to_plink.errors.log
#$ -q long.qf
#$ -l h_rt=100:00:00
#$ -V
#$ -pe shmem 20
#$ -P lindgren.prjc
#$ -t 24

readonly CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi

readonly WD=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink

# input paths
readonly VCF=/well/lindgren/UKBIOBANK/saskia/vcf_Filter/ukbb-wes-oqfe-pvcf-chr${CHR}_filtered.vcf.gz
readonly FASTA=/well/lindgren/saskia/WES/trial/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# output paths
readonly PREFIX=ukb_wes_200k_chr${CHR}
readonly TMP_VCF=${WD}/vcf/tmp-${PREFIX}.vcf.gz # intermediate VCF before going to PLINK
readonly BFILE=${WD}/plink/${PREFIX} # PLINK bfile prefix

time_check() {
  echo -e "\n########\n$1 (job id: ${JOB_ID}.${SGE_TASK_ID}, $(date))\n########"
}

# check that VCFs are not truncated
vcf_check() {
  if [ ! -f $1 ]; then
    >&2 time_check "Error: $1 does not exist. Exiting.\n"
    exit 1
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    >&2 time_check "Error: $1 may be truncated. Exiting.\n"
    exit 1
  fi
}

time_check "chr${CHR} start"

if [ ! -f ${TMP_VCF} ]; then
  vcf_check ${VCF}
  time_check "chr${CHR} start temporary VCF writing"
  bcftools norm -Oz -f ${FASTA} -m -any --threads 20 -o ${TMP_VCF} ${VCF}
  vcf_check ${TMP_VCF}
  time_check "chr${CHR} temporary VCF finished writing"
else
  vcf_check ${TMP_VCF}
  time_check "Warning: ${TMP_VCF} already exists, skipping multiallelic split and indel normalisation"
fi

# NOTE: We only check for the bed file for simplicity
if [ ! -f ${BFILE}.bed ]; then
  time_check "chr${CHR} start PLINK files writing"
  plink2 --vcf ${TMP_VCF} \
    --keep-allele-order \
    --double-id \
    --allow-extra-chr \
    --memory 40000 \
    --make-bed \
   --out ${BFILE}

  if [ -f ${BFILE}.bed ]; then
    time_check "chr${CHR} PLINK files finished writing"
  else
    >&2 time_check "Error: chr${CHR} PLINK files were not written successfully. Exiting.\n"
    exit 1
  fi
else
  time_check "Warning: ${BFILE}.bed already exists, skipping VCF to PLINK conversion"
fi

# create variant IDs
if [ -f ${BFILE}.bim ]; then
  if [ -f ${BFILE}.bim_new ]; then
    >&2 time_check "Error: ${BFILE}.bim_new already exists. Exiting.\n"
    exit 0
  else
    python3 ${WD}/scripts/make_variant_ids.py --bfile ${BFILE}
  fi
else
  >&2 time_check "Error: ${BFILE}.bim does not exist, new variant IDs cannot be created. Exiting.\n"
  exit 1
fi

if [ $? -eq 0 ]; then
  if [[ ! -f ${BFILE}.bim_old  && -f ${BFILE}.bim_new ]]; then # add a somewhat redundant check to be sure that .bim_new exists
    mv ${BFILE}.bim ${BFILE}.bim_old && \
      ln -s -f ${BFILE}.bim_new ${BFILE}.bim # create symlink instead of changing filename to indicate that we're using the .bim_new file as a replacement
  else
    echo "Error: Check that ${BFILE}.bim_old does not exist and that ${BFILE}.bim_new exists"
  fi
fi

time_check "chr${CHR} VCF to PLINK script complete"
