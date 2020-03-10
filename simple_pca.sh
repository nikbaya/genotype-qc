#!/usr/bin/env bash
# This script runs PLINK PCA on a dataset
#
# Overview:
# 1) Basic variant and sample QC to prepare for LD pruning
# 2) LD prune
# 3) Remove related samples
# 3) Run PCA
#
# Author: Nikolas Baya, Mar 2020

# TODO: FID/IIDs with asterisk do not work with grep
# To fix this, could sed all asterisks with random string, then sed after grep step to return the asterisks

bfile=${1?Missing param: No bfile prefix provided} # PLINK bfile prefix

# TODO: Set working directory

# Check that input files exist
input_files=( ${bfile}.{bed,bim,fam} )
for file in ${input_files[@]}; do
  test ! -f ${file} && echo "Error: ${file} does not exist" && exit 1
done

log=${bfile}.simple_pca.log # for redirecting all stdout, stderr

echo -e "...Performing basic variant & sample QC..." | tee -a ${log}

plink --bfile ${bfile} --silent --geno 0.05 --make-bed --out tmp0_${bfile} >> ${log} 2>&1

plink --bfile tmp0_${bfile} --silent --mind 0.02 --make-bed --out tmp1_${bfile} >> ${log} 2>&1

plink --bfile tmp1_${bfile} --silent --geno 0.02 --make-bed --out tmp2_${bfile} >> ${log} 2>&1

plink --bfile tmp2_${bfile} --silent --hwe 1e-6 --make-bed --out tmp3_${bfile} >> ${log} 2>&1

plink --bfile tmp3_${bfile} --silent --maf 0.05 --make-bed --out tmp4_${bfile} >> ${log} 2>&1

tmp4_files=( tmp4_${bfile}.{bed,bim,fam} )
for file in ${tmp4_files[@]}; do
  test ! -f ${file} && echo "Error: Basic variant and sample QC step failed." | tee -a ${log} && exit 1
done

echo -e "...LD pruning..." | tee -a ${log}

# LD pruning
ld_th=0.2
ld_wind=200
ld_move=$(($ld_wind/2))
maxiter=100 # maximum number of pruning iterations before timing out
i=1

plink --bfile tmp4_${bfile} --indep-pairwise $ld_wind $ld_move $ld_th  --silent --allow-no-sex --make-founders require-2-missing --out tmp4_${bfile}_prune$i >> ${log} 2>&1

nprune_old=$( wc -l < tmp4_${bfile}.bim )
nprune_new=$( wc -l < tmp4_${bfile}_prune$i.prune.in )

while [ ${nprune_old} -gt ${nprune_new} ]; do
  i=$(( i+1 ))
  echo "Pruning pass $i" >> ${log} 2>&1
  plink --bfile tmp4_${bfile} --extract tmp4_${bfile}_prune$((i-1)).prune.in --indep-pairwise $ld_wind $ld_move $ld_th --silent --allow-no-sex --make-founders require-2-missing --out tmp4_${bfile}_prune$i >> ${log} 2>&1

  nprune_old=$nprune_new
  nprune_new=$( wc -l < tmp4_${bfile}_prune$i.prune.in )
  echo "nprune_old: ${nprune_old}, nprune_new: ${nprune_new}" >> ${log} 2>&1
  
  if [[ $i -eq ${maxiter} ]]; then
    echo "Error: LD pruning step failed to converge after {maxiter} iterations." | tee -a ${log}  && exit 1
  fi
done

test ! -f tmp4_${bfile}_prune$i.prune.in  && echo "Error: LD pruning step failed." | tee -a ${log}  && exit 1

# Run PLINK IBD
echo -e "...Extracting LD pruned set and running IBD..." | tee -a ${log}
pihat_thresh=0.2

plink --bfile tmp4_${bfile} --extract tmp4_${bfile}_prune$i.prune.in --silent --genome --min ${pihat_thresh} --make-bed --allow-no-sex --out tmp4_${bfile}_finalpruned >> ${log} 2>&1

test ! -f tmp4_${bfile}_finalpruned.genome  && echo "Error: PLINK IBD (--genome) step failed." | tee -a ${log} && exit 1

# Removing related samples
echo -e "...Removing related samples..." | tee -a ${log}

ibd=tmp4_${bfile}_finalpruned.genome
cat <( awk '{print $1,$2}' ${ibd} ) <(awk '{print $3,$4}' ${ibd} ) \
  | grep -v "FID.*IID" \
  | sort \
  | uniq -c \
  | sort -rk 1,1 > tmp_${bfile}.related_ct.txt

related_fname=${bfile}.remove.related.txt

if [[ $( cat tmp_${bfile}.related_ct.txt | wc -l ) -eq 0 ]]; then
  touch ${related_fname}
  echo "No related samples identified by PLINK IBD" >> ${log}
  
else
  i=0
  ibd_idx=0
  ibd_curr=${ibd}
  old_grep_str=""
  old_line_ct=$(grep -v 'FID1.*IID1.*FID2.*IID2' ${ibd_curr} | wc -l)
  while read line; do
    arr=(${line})
    new_grep_str="${arr[1]}.*${arr[2]}"
  
    if [[ ${old_grep_str} == "" ]]; then
      grep_str=${new_grep_str}
    else
      grep_str="${old_grep_str}\|${new_grep_str}"
    fi
  
    old_grep_str=${grep_str}
    new_line_ct=$(grep -v ${grep_str} ${ibd_curr} | grep -v 'FID1.*IID1.*FID2.*IID2' | wc -l)
    if [[ ${new_line_ct} -eq 0 ]]; then
      break
    elif [[ ${new_line_ct} -lt ${old_line_ct} ]]; then
      echo -e "${arr[1]} ${arr[2]}" >> ${related_fname}
      old_line_ct=${new_line_ct}
    fi
  
    i=$(($i+1))
  
    if [[ $i -eq 100 ]]; then # every 100 samples, create new filtered .genome file and re-initialize grep str to speed up process
      ibd_new=tmp${ibd_idx}_${bfile}.removerelateds.genome
      grep -v ${grep_str} ${ibd_curr} > ${ibd_new}
      ibd_curr=${ibd_new}
      ibd_idx=$((${ibd_idx}+1))
      old_grep_str=""
      i=0
    fi

  done < tmp_${bfile}.related_ct.txt

  test ! -f ${related_fname}  && echo "Error: Obtaining list of related samples failed." | tee -a ${log} && exit 1
  
fi

echo -e "...Running PCA..." | tee -a ${log}
n_pcs=20

plink --bfile tmp4_${bfile} --remove ${related_fname} --pca ${n_pcs} header --silent --out tmp_${bfile}.simple_pca >> ${log} 2>&1

pca_files=( tmp_${bfile}.simple_pca.eigenv{al,ec} )
for file in ${pca_files[@]}; do
  test ! -f ${file} && echo "Error: PCA failed." | tee -a ${log} && exit 1
done

# remove intermediate temporary files (separated for readability, redirect stderr to /dev/null)
rm tmp[0-5]_${bfile}.{bed,bim,fam,hh,log,irem,nosex} 2> /dev/null
rm tmp4_${bfile}_prune*.{hh,log,prune.{in,out},nosex} 2> /dev/null
rm tmp4_${bfile}_finalpruned.{bed,bim,fam,hh,log,genome,nosex} 2> /dev/null
rm tmp_${bfile}.related_ct.txt 2> /dev/null
rm tmp*_${bfile}.removerelateds.genome 2> /dev/null

mv tmp_${bfile}.simple_pca.eigenval ${bfile}.simple_pca.eigenval
mv tmp_${bfile}.simple_pca.eigenvec ${bfile}.simple_pca.eigenvec
rm tmp_${bfile}.simple_pca.{hh,log,nosex} 2> /dev/null

echo "PCA complete." | tee -a ${log}
