#!/usr/bin/env bash
#
#
# Author: Nik Baya (2021-01-26)
#
#$ -N gq_stats
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/test/gq_stats.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/test/gq_stats.errors.log
#$ -q short.qf
#$ -V
#$ -P lindgren.prjc
#$ -t 1-24

CHR=${SGE_TASK_ID}
TYPE=$1 # either "ref" (hom ref) , "alt" (hom alt), or "het"

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vqsr/test" # TODO

IN=/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr${CHR}.vcf.gz

hom_alt=${WD}/test_hom_alt_gq_chr${CHR}.txt
hom_ref=${WD}/test_subset_hom_ref_gq_chr${CHR}.txt
het=${WD}/test_subset_het_gq_chr${CHR}.txt

if [[ ! -f ${hom_alt} && "$TYPE" == "alt" ]]; then
  bcftools query \
    ${IN}  \
    -i 'GT="AA"' \
    -f '[ %CHROM %POS %SAMPLE %GQ\n ]' > $hom_alt
else
  subset=${WD}/subset_chr${CHR}.txt
  if [ ! -f $subset ]; then
    bcftools view -H -G ${IN} | awk '{ print $1"\t"$2 }' | shuf -n 5000 | sort -k2 > $subset
  fi
  if [[ ! -f ${hom_ref} && "$TYPE" == "ref" ]]; then
    bcftools query \
      ${IN} \
      -R $subset \
      -i 'GT="RR"' \
      -f '[ %CHROM %POS %SAMPLE %GQ\n ]' > $hom_ref
  elif [[ ! -f ${het} && "$TYPE" == "het" ]]; then
    bcftools query \
      ${IN} \
      -R $subset \
      -i 'GT="het"' \
      -f '[ %CHROM %POS %SAMPLE %GQ\n ]' > $het
  fi
fi
