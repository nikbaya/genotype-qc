#!/bin/bash
# Simple QC using PLINK
# Intended to recreate Ricopili QC pipeline

BFILE=${1?ERROR: No bfile given}

QC_WD=simple_qc_${BFILE}  #${BFILE}_simple_qc

if [ ! -d ${QC_WD} ]; then
	mkdir $QC_WD
fi

## QC thresholds (uses same naming convention as Ricopili)

PRE_GENO=0.05 # include only SNPs with missing-rate < $PRE_GENO (before ID filter), important for post merge of multiple platforms
MIND=0.02     # include only IDs with missing-rate < $MIND
FHET_TH=0.2   # include only samples with -$FHET_TH < FHET < $FHET_TH
GENO=0.02     # include only SNPs with missing-rate < $GENO
HWE_TH=1e-10  # include only SNPs with HWE p-val > $HWE_TH
MAF=0.005     # include only SNPs with minor allele frequency < $MAF


count () {
	# TODO: Only display change in samples/variants
	echo -e "\tsamples:  $( cat ${1}.fam | wc -l )"
	echo -e "\tvariants: $( cat ${1}.bim | wc -l )"
}

## SNP call rate 1st pass filtering
echo "...SNP call rate 1st pass filtering ($PRE_GENO)..."
plink --bfile ${BFILE} --geno $PRE_GENO --make-bed --out ${QC_WD}/${BFILE}_pregeno >> ${QC_WD}/simple_qc_${BFILE}.log
count ${QC_WD}/${BFILE}_pregeno

## Sample call rate filtering
echo "...Sample call rate filtering ($MIND)..."
plink --bfile ${QC_WD}/${BFILE}_pregeno --mind ${MIND} --make-bed --out ${QC_WD}/${BFILE}_mind >> ${QC_WD}/simple_qc_${BFILE}.log
count ${QC_WD}/${BFILE}_mind

## Inbreeding coefficient (Fhet)
echo "...Sample inbreeding coefficient filtering ($FHET_TH)..."
plink --bfile ${QC_WD}/${BFILE}_mind --het --out ${QC_WD}/${BFILE} >> ${QC_WD}/simple_qc_${BFILE}.log
awk -v FHET_TH=${FHET_TH} '{ if ($6 < -FHET_TH || $6> FHET_TH) print $1, $2, $6 }' ${QC_WD}/${BFILE}.het > ${QC_WD}/${BFILE}.remove.fhet.txt
plink --bfile ${QC_WD}/${BFILE}_mind --remove ${QC_WD}/${BFILE}.remove.fhet.txt --make-bed --out ${QC_WD}/${BFILE}_fhet >> ${QC_WD}/simple_qc_${BFILE}.log
count ${QC_WD}/${BFILE}_fhet

## SNP call rate 2nd pass filtering
echo "...SNP call rate 2nd pass filtering ($GENO)..."
plink --bfile ${QC_WD}/${BFILE}_fhet --geno $GENO --make-bed --out ${QC_WD}/${BFILE}_geno >> ${QC_WD}/simple_qc_${BFILE}.log
count ${QC_WD}/${BFILE}_geno

## MAF filtering
echo "...SNP MAF filtering ($MAF)..."
plink --bfile ${QC_WD}/${BFILE}_geno --maf $MAF --make-bed --out ${QC_WD}/${BFILE}_maf >> ${QC_WD}/simple_qc_${BFILE}.log
count ${QC_WD}/${BFILE}_maf

## HWE filtering
echo "...SNP HWE filtering ($HWE_TH)..."
plink --bfile ${QC_WD}/${BFILE}_maf --hwe $HWE_TH --make-bed --out ${QC_WD}/${BFILE}_postqc >> ${QC_WD}/simple_qc_${BFILE}.log
count ${QC_WD}/${BFILE}_postqc
