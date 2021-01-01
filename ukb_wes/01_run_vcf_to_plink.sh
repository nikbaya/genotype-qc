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
#$ -q short.qe
#$ -l h_vmem=14G
#$ -l h_rt=29:00:00
#$ -V
#$ -P lindgren.prjc
#$ -t 1-22

chrom=${SGE_TASK_ID}

wd=/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink

# input paths
vcf=/well/lindgren/UKBIOBANK/saskia/vcf_Filter/ukbb-wes-oqfe-pvcf-chr${chrom}_filtered.vcf.gz
fasta=/well/lindgren/saskia/WES/trial/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# output paths
prefix=ukb_wes_200k_chr${chrom}
tmp_vcf=${wd}/vcf/tmp-${prefix}.vcf.gz # intermediate VCF before going to PLINK
bfile=${wd}/plink/${prefix} # PLINK bfile prefix

timecheck() {
  echo -e "\n########\n$1 ($(date))\n########"
}

# check that VCFs are not truncated
if [ $( bcftools view ${vcf} 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
  >&2 echo -e "WARNING: ${vcf} may be truncated\nSkipping this VCF\n"
  exit 1
fi

timecheck "chr${chrom} start"

if [ ! -f ${tmp_vcf} ]; then
  bcftools norm -Ou -m -any ${vcf} | # split sites with multiple alleles
    bcftools norm -Oz -f ${fasta} -o ${tmp_vcf} # left-align and normalise indels using fasta file

  timecheck "chr${chrom} temporary VCF finished writing"

else
  echo "${tmp_vcf} already exists, skipping multiallelic split and indel normalisation"
fi

# NOTE: We only check for the bed file for simplicity
if [ ! -f ${bfile}}.bed ]; then
  plink2 --vcf ${tmp_vcf} \
    --keep-allele-order \
    --double-id \
    --allow-extra-chr \
    --memory 13500 \
    --make-bed \
    --out ${bfile}

  timecheck "chr${chrom} PLINK files finished writing"

else
  echo "${bfile}.bed already exists, skipping VCF to PLINK conversion"
fi

# update variant IDs
python3 ${wd}/scripts/make_variant_ids.py --chr ${chrom}

if [ $? -eq 0 ]; then
  if [[ ! -f ${bfile}.bim_old  && -f ${bfile}.bim_new ]]; then # add a somewhat redundant check to be sure that .bim_new exists
    mv ${bfile}.bim ${bfile}.bim_old && \
      ln -s -f ${bfile}.bim_new ${bfile}.bim # create symlink instead of changing filename to indicate that we're using the .bim_new file as a replacement
  else
    echo "WARNING: Check that ${bfile}.bim_old does not exist and that ${bfile}.bim_new exists"
  fi
fi

timecheck "chr${chrom} VCF to PLINK script complete"
