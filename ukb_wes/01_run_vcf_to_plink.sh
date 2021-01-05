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
#$ -pe shmem 10
#$ -P lindgren.prjc
#$ -t 4,8

if [ ${SGE_TASK_ID} -eq 23 ]; then
  chrom="X"
else
 chrom=${SGE_TASK_ID}
fi

wd=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink

# input paths
vcf=/well/lindgren/UKBIOBANK/saskia/vcf_Filter/ukbb-wes-oqfe-pvcf-chr${chrom}_filtered.vcf.gz
fasta=/well/lindgren/saskia/WES/trial/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# output paths
prefix=ukb_wes_200k_chr${chrom}
tmp_vcf=${wd}/vcf/tmp-${prefix}.vcf.gz # intermediate VCF before going to PLINK
bfile=${wd}/plink/${prefix} # PLINK bfile prefix

time_check() {
  echo -e "\n########\n$1 (job id: ${JOB_ID}.${SGE_TASK_ID}, $(date))\n########"
}

# check that VCFs are not truncated
vcf_check() {
  if [ ! -f $1 ]; then
    >&2 time_check "Error: $1 does not exist. Exiting.\n"
    exit 1
  else
    if [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
      >&2 time_check "Error: $1 may be truncated. Exiting.\n"
      exit 1
    fi
  fi
}

time_check "chr${chrom} start"

if [ ! -f ${tmp_vcf} ]; then
  set -x
  vcf_check ${vcf}
  bcftools norm -Ou -m -any ${vcf} | \
    bcftools norm -Oz -f ${fasta} -o ${tmp_vcf} --threads 10 # left-align and normalise indels using fasta file
  vcf_check ${tmp_vcf}
  time_check "chr${chrom} temporary VCF finished writing"
  set +x
else
  vcf_check ${tmp_vcf}
  time_check "Warning: ${tmp_vcf} already exists, skipping multiallelic split and indel normalisation"
fi

# NOTE: We only check for the bed file for simplicity
if [ ! -f ${bfile}.bed ]; then
  plink2 --vcf ${tmp_vcf} \
    --keep-allele-order \
    --double-id \
    --allow-extra-chr \
    --memory 40000 \
    --make-bed \
   --out ${bfile}

  if [ -f ${bfile}.bed ]; then
    time_check "chr${chrom} PLINK files finished writing"
  else
    >&2 time_check "Error: chr${chrom} PLINK files were not written successfully. Exiting.\n"
    exit 1
  fi
else
  time_check "Warning: ${bfile}.bed already exists, skipping VCF to PLINK conversion"
fi

# create variant IDs
if [ -f ${bfile}.bim ]; then
  if [ -f ${bfile}.bim_new ]; then
    >&2 time_check "Error: ${bfile}.bim_new already exists. Exiting.\n"
    exit 0
  else
    python3 ${wd}/scripts/make_variant_ids.py --bfile ${bfile}
  fi
else
  >&2 time_check "Error: ${bfile}.bim does not exist, new variant IDs cannot be created. Exiting.\n"
  exit 1
fi

if [ $? -eq 0 ]; then
  if [[ ! -f ${bfile}.bim_old  && -f ${bfile}.bim_new ]]; then # add a somewhat redundant check to be sure that .bim_new exists
    mv ${bfile}.bim ${bfile}.bim_old && \
      ln -s -f ${bfile}.bim_new ${bfile}.bim # create symlink instead of changing filename to indicate that we're using the .bim_new file as a replacement
  else
    echo "Error: Check that ${bfile}.bim_old does not exist and that ${bfile}.bim_new exists"
  fi
fi

time_check "chr${chrom} VCF to PLINK script complete"
