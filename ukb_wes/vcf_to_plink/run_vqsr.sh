#!/usr/bin/env bash
#
# Author: Nik Baya (2021-01-16)
#
#$ -N vqsr
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/vqsr.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/vqsr.errors.log
#$ -q short.qf
#$ -pe shmem 40
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
#readonly FASTA="/well/lindgren/saskia/WES/trial/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
readonly IN="${WD}/test_chr${CHR}.vcf.gz"
readonly OUT="${WD}/test_vqsr_chr${CHR}"

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

make_index() {
  if [ ! -f $1.tbi ]; then
    tabix -p vcf $1
  fi
}

SECONDS=0

make_index ${IN}

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

gatk --java-options "-Xmx100g -Xms100g" VariantAnnotator \
  -V ${IN} \
  -A ExcessHet -A MappingQualityRankSumTest -A ReadPosRankSumTest \
  -A QualByDepth -A StrandOddsRatio -A QualByDepth -A FisherStrand -A MappingQuality \
  -O ${OUT}_annot.vcf.gz

#make_index ${OUT

gatk --java-options "-Xmx100g -Xms100g" VariantFiltration \
  -V ${IN}_annot.vcf.gz \
  --filter-expression "ExcessHet > 54.69" \
  --filter-name ExcessHet \
  -O ${OUT}_filter_excesshet.vcf.gz

#make_index

gatk MakeSitesOnlyVcf \
  -I ${OUT}_filter_excesshet.vcf.gz \
  -O ${OUT}_sitesonly.vcf.gz

readonly REF="/well/lindgren/UKBIOBANK/nbaya/resources/ref"

gatk --java-options "-Xmx100g -Xms100g" VariantRecalibrator \
   -V cohort_sitesonly.vcf.gz \
   --trust-all-polymorphic \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
   -mode SNP \
   --max-gaussians 6 \
   -resource:hapmap,known=false,training=true,truth=true,prior=15:${REF}/hapmap_3.3.hg38.vcf.gz \
   -resource:omni,known=false,training=true,truth=true,prior=12:${REF}/1000G_omni2.5.hg38.vcf.gz \
   -resource:1000G,known=false,training=true,truth=false,prior=10:${REF}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   -resource:dbsnp,known=true,training=false,truth=false,prior=7:${REF}/Homo_sapiens_assembly38.dbsnp138.vcf \
   -O ${OUT}_recal_snp \
   --tranches-file ${OUT}_recal_snp.tranches


duration=${SECONDS}
echo "chr${CHR} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
