#!/usr/bin/env bash
#
# This script is for calculating genotype missingness
#
# Author: Nik Baya (2021-01-07)
#
#$ -N calc_missing_rate
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/calc_missing_rate.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/calc_missing_rate.errors.log
#$ -q long.qf
#$ -l h_rt=50:00:00
#$ -pe shmem 10
#$ -V
#$ -P lindgren.prjc
#$ -t 13,18

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink"
readonly VCF="${WD}/vcf/tmp-ukb_wes_200k_chr${CHR}.vcf.gz"
readonly BFILE="${WD}/plink/ukb_wes_200k_chr${CHR}" # bfile prefix for per-chrom files
readonly OUT="${WD}/test/missing_rate_chr${CHR}"

plink_files_exist() {
  [ ` ls -1 $1.{bed,bim,fam} 2> /dev/null | wc -l ` -eq 3 ]
}

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

SECONDS=0
if plink_files_exist ${BFILE} ; then
  plink --bfile ${BFILE} \
    --missing \
    --out ${OUT}
fi
duration=${SECONDS}
echo "finished calculating plink missingness, $( elapsed_time ${duration} ) ($( date ))"

SECONDS=0
# from https://gist.github.com/darencard/a91e2bedf9f81296d201cfb612e0337d
paste \
<(bcftools query -f '[%SAMPLE\t]\n' ${VCF} | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' ${VCF} | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2) > ${OUT}.txt
duration=${SECONDS}
echo "finished calculating vcf missingness, $( elapsed_time ${duration} ) ($( date ))"
