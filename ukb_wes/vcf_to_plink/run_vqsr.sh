#!/usr/bin/env bash
#
# Author: Nik Baya (2021-01-16)
#
#$ -N vqsr
#$ -o /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/test_vqsr.log
#$ -e /well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test/test_vqsr.errors.log
#$ -q short.qf
#$ -pe shmem 40
#$ -V
#$ -P lindgren.prjc
#$ -t 20-22

CHR=${SGE_TASK_ID}

if [ ${CHR} -eq 23 ]; then
  CHR="X"
elif [ ${CHR} -eq 24 ]; then
  CHR="Y"
fi
readonly CHR

readonly WD="/well/lindgren/UKBIOBANK/nbaya/wes_200k/vcf_to_plink/test"
#readonly FASTA="/well/lindgren/saskia/WES/trial/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
#readonly IN="${WD}/tmp-ukb_wes_200k_chr${CHR}.vcf.gz"
readonly IN=${WD}/test_ukb_wes_chr${CHR}.vcf.gz
#readonly OUT="${WD}/test_full_200k_vqsr_chr${CHR}"
readonly OUT=${WD}/test_ukb_chr${CHR}

raise_error() {
  >&2 echo -e "Error: $1. Exiting."
  exit 1
}

elapsed_time() {
  echo "elapsed time: $( echo "scale=2; $1/3600" | bc -l ) hrs"
}

vcf_check() {
  if [ ! -f $1 ]; then
    >&2 "Error: $1 does not exist. Exiting.\n"
    exit 1
  elif [ $( bcftools view $1 2>&1 | head | grep "No BGZF EOF marker" | wc -l ) -gt 0 ]; then
    >&2 "Error: $1 may be truncated. Exiting.\n"
    exit 1
  fi
}

make_index() {
  if [ ! -f $1.tbi ]; then
    tabix -p vcf $1
  fi
}

SECONDS=0

#if [ ! -f ${IN} ]; then
#bcftools view -Oz \
#  -r chr13:19177404-19179298 \
#  ${WD}/tmp-ukb_wes_200k_chr13.vcf.gz \
#  --threads 39 \
#  -o ${IN}
#fi

vcf_check ${IN}
make_index ${IN}

module load GATK/4.1.7.0-GCCcore-8.3.0-Java-11

#if [ ! -f ${OUT}_annot.vcf.gz ]; then
#gatk --java-options "-Xmx${MEM}g -Xms${MEM}g" VariantAnnotator \
#  -V ${IN} \
#  -A ExcessHet -A MappingQualityRankSumTest -A ReadPosRankSumTest -A RMSMappingQuality \
#  -A StrandOddsRatio -A QualByDepth -A FisherStrand -A MappingQuality -A DepthPerSampleHC \
#  -O ${OUT}_annot.vcf.gz
#fi

readonly MEM=120 # memory in gb used for gatk

if [ ! -f ${OUT}_annot.vcf.gz ]; then
gatk --java-options "-Xmx${MEM}g -Xms${MEM}g" VariantAnnotator \
  -V ${IN} \
  -A ExcessHet -A StrandOddsRatio -A QualByDepth -A FisherStrand \
  -O ${OUT}_annot.vcf.gz
fi

vcf_check ${OUT}_annot.vcf.gz
make_index ${OUT}_annot.vcf.gz

if [ ! -f ${OUT}_filter_excesshet.vcf.gz ]; then
gatk --java-options "-Xmx${MEM}g -Xms${MEM}g" VariantFiltration \
  -V ${OUT}_annot.vcf.gz \
  --filter-expression "ExcessHet > 54.69" \
  --filter-name ExcessHet \
  -O ${OUT}_filter_excesshet.vcf.gz
fi

vcf_check ${OUT}_filter_excesshet.vcf.gz
make_index ${OUT}_filter_excesshet.vcf.gz

if [ ! -f ${OUT}_sitesonly.vcf.gz ]; then
gatk MakeSitesOnlyVcf \
  -I ${OUT}_filter_excesshet.vcf.gz \
  -O ${OUT}_sitesonly.vcf.gz
fi

vcf_check ${OUT}_sitesonly.vcf.gz
make_index ${OUT}_sitesonly.vcf.gz

readonly REF="/well/lindgren/UKBIOBANK/nbaya/resources/ref"

if [ ! -f ${OUT}_recal_snp.tranches ]; then
gatk --java-options "-Xmx${MEM}g -Xms${MEM}g" VariantRecalibrator \
   -V ${OUT}_sitesonly.vcf.gz \
   --trust-all-polymorphic \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
   -an QD -an FS -an SOR \
   -mode SNP \
   --max-gaussians 4 \
   -resource:hapmap,known=false,training=true,truth=true,prior=15 ${REF}/hapmap_3.3.hg38.vcf.gz \
   -resource:omni,known=false,training=true,truth=true,prior=12 ${REF}/1000G_omni2.5.hg38.vcf.gz \
   -resource:1000G,known=false,training=true,truth=false,prior=10 ${REF}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${REF}/Homo_sapiens_assembly38.dbsnp138.vcf \
   -O ${OUT}_recal_snp \
   --tranches-file ${OUT}_recal_snp.tranches
fi

duration=${SECONDS}
echo "chr${CHR} finished, $( elapsed_time ${duration} ) (job id: ${JOB_ID}.${SGE_TASK_ID} $( date ))"
