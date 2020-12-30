#!/usr/bin/env bash
#
# This script is for splitting multi-nucleotide polymorphisms,
# then converting the VCF to PLINK format.
#
# Author: Nik Baya (2020-12-23)
#
#$ -N vcf_to_plink
#$ -o vcf_to_plink.log
#$ -e vcf_to_plink.errors.log
#$ -q short.qe
#$ -wd /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test
#$ -l h_vmem=14G
#$ -l h_rt=20:00:00
#$ -V
#$ -P lindgren.prjc
#$ -t 1-22


chrom=${SGE_TASK_ID}

# print date to log
echo -e "\n########"  && \
date  && \
echo "chr${chrom} start"  && \
echo -e "########\n"


vcf=ukbb-wes-oqfe-pvcf-chr${chrom}_filtered.vcf.gz
out=test_ukb_wes_chr${chrom}

if [ ! -f ${out}.vcf.gz ]; then

bcftools norm -Ou -m -any ${vcf} |
  #bcftools norm -Ou -f human_g1k_v37.fasta |
  bcftools annotate -Oz -x ID \
    -I +'%CHROM:%POS:%REF:%ALT' > ${out}.vcf.gz

echo -e "\n########"  && \
date  && \
echo "chr${chrom} VCF finished"  && \
echo -e "########\n"

fi

if [ ! -f ${out}.vcf.gz.tbi ]; then

tabix -p vcf ${out}.vcf.gz

echo -e "\n########"  && \
date  && \
echo "chr${chrom} tabix finished"  && \
echo -e "########\n"

fi

if [ ! -f ${out}.bed ]; then

plink2 --vcf ${out}.vcf.gz \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr  \
    --memory 13000 \
    --make-bed \
    --out ${out}

echo -e "\n########" && \
date && \
echo "chr${chrom} PLINK finished" && \
echo -e "########\n"

fi
