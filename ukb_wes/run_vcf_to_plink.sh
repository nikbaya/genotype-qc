#!/usr/bin/env bash
#
# This script is for splitting multi-nucleotide polymorphisms,
# then converting the VCF to PLINK format.
#
# Author: Nik Baya (2021-1-1)
#
#$ -N test_vcf_to_plink
#$ -o test_vcf_to_plink.log
#$ -e test_vcf_to_plink.errors.log
#$ -q short.qe
#$ -wd /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test
#$ -l h_vmem=14G
#$ -l h_rt=29:00:00
#$ -V
#$ -P lindgren.prjc
#$ -t 20-22


chrom=${SGE_TASK_ID}

# print date to log
echo -e "\n########"  && \
date  && \
echo "chr${chrom} start"  && \
echo -e "########\n"

vcf=ukb_test_data_chr${chrom}.vcf.gz
fasta=/well/lindgren/saskia/WES/trial/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
out=test_ukb_wes_chr${chrom}

if [ ! -f ${out}.vcf.gz ]; then

bcftools norm -Ou -m -any ${vcf} |
  bcftools norm -Ou -f ${fasta} -o ${out}.vcf.gz

echo -e "\n########"  && \
date  && \
echo "chr${chrom} VCF finished"  && \
echo -e "########\n"

fi

tabix -f -p vcf ${out}.vcf.gz

echo -e "\n########"  && \
date  && \
echo "chr${chrom} tabix finished"  && \
echo -e "########\n"


if [ ! -f ${out}.bed ]; then

plink2 --vcf ${out}.vcf.gz \
    --keep-allele-order \
    --double-fid \
    --allow-extra-chr  \
    --memory 13000 \
    --make-bed \
    --out ${out}

echo -e "\n########" && \
date && \
echo "chr${chrom} PLINK finished" && \
echo -e "########\n"

fi
