#! /bin/bash
# Used to filter by PCA results
#
# Example:
#	./filter_by_pc.sh alc_stpe1_eur_nb-qc2 stpe1_pca2 2 false true
# Uses PCs from:
#	stpe1_pca2.menv.mds
#	stpe1_pca2_1kg.menv.mds
# Uses removed IDs from:
#	stpe1_pca2.mepr.famex
# Creates the following files (without Ricopili prefix removed, REMOVE_PREFIX=false):
#	alc_stpe1_eur_nb-qc2.remove.withoutref.txt
#	alc_stpe1_eur_nb-qc2.remove.withref.txt
#	alc_stpe1_eur_nb-qc2.remove.txt
# Using created files, performs PLINK filter on:
#	alc_stpe1_eur_nb-qc2.bed
#	alc_stpe1_eur_nb-qc2.bim
#	alc_stpe1_eur_nb-qc2.fam
# Output of PLINK filter (suffix=$STUDY.post_pca${PCA_ROUND}):
#	alc_stpe1_eur_nb-qc2.post_pca2.bed
#	alc_stpe1_eur_nb-qc2.post_pca2.bim
#	alc_stpe1_eur_nb-qc2.post_pca2.fam


# Arguments
STUDY=${1?ERROR: No study name given} # name of study for PCA results (*.menv.mds files)
BFILE=${2?ERROR: No bfile prefix given} # bfile from which PC outliers are removed (often the same as $STUDY)

# Arguments with defaults
DEF_PCA_ROUND=1 #set default PCA round
PCA_ROUND=${3:-${DEF_PCA_ROUND}} # which PCA round this is from, used for output (default: 1)
REMOVE_PREFIX=${4:true} # if true, removes prefix from .remove.txt files before filtering bfile (default: true)
PLINK_REMOVE=${5:true} # if true, uses PLINK to filter $BFILE using the $STUDY.remove* files created by this script (default: true)

PCA=${STUDY}.menv.mds
PCA_W_REF=${STUDY}_1kg.menv.mds
FAMEX=${STUDY}.mepr.famex
OUT=${BFILE}.post_pca${PCA_ROUND} # bfile prefix for filtered PLINK files output by this script

# Filtering arguments for PCA w/o ref, PCA w/ ref
# NOTE: For min_cutoffs or max_cutoffs use "NA" if there is no cutoff for that PC
# If either no min cutoffs or max cutoffs are required, the array can be made empty

# Cutoffs for PCA w/o ref (EDIT THESE VALUES)
# NOTE: These array variables cannot have "_" in their name in order to expand them inside functions
declare mincutoffs=( -0.11 "NA" "NA" -0.1) # min PC cutoffs, exclude individuals with PC1 < min_cutoffs[0], PC2 < min_cutoffs[1], etc.
declare maxcutoffs=( "NA" "NA" 0.1 0.1 ) # max PC cutoffs, exclude individuals with PC1 > max_cutoffs[0], PC2 > max_cutoffs[1], etc.
# Cutoffs for PCA w/o ref (EDIT THESE VALUES)
declare  mincutoffswref=( "NA" "NA" "NA" 0.002 ) # min PC cutoffs for PCA w/ ref
declare  maxcutoffswref=( -0.011 )


# ---------------------- NO NEED TO EDIT BELOW THIS LINE ----------------------

files=( $PCA $PCA_W_REF $FAMEX )
for file in ${files[@]}; do
  test ! -f $file && echo "ERROR: File $file not found." && exit 1
done

# TODO: Print arguments/input files used
# TODO: Check the max length of min_cutoffs/max_cutoffs arrays

# Start filtering
get_max () {
  # get max of two numbers
  if [[ $1 -gt $2 ]]; then
    echo $1
  else
    echo $2
  fi
}

awk_filter () {
  # used for performing awk filter on a PC in a PCA results file

  # load arguments
  local pca_fname=$1
  local tmp_fname=$2
  local pc=$3 # which PC to filter on
  local cutoff_type=$4 # used for print statements
  declare cutoffs=("${!5}") # array of cutoffs

  # initialize variables
  local -i cutoffs_len=${#cutoffs[@]} # length of array
  local pc_idx=$(( ${pc} - 1 )) # index of PC in cutoffs array
  local awk_pc_idx=$(( ${pc} + 3)) # PC1 starts is at column 4 in .mds file, PC2 at col 5, etc

  if [[ ${cutoffs_len} -ge ${pc} ]]; then
    cutoff=${cutoffs[$pc_idx]} # threshold for a given PC

    if [[ ${cutoff} == "NA" ]]; then
      echo -e "\tNo ${cutoff_type} cutoff"
    else
      local wc_old=$( cat ${tmp_fname} | wc -l ) # count lines previous to adding new samples to remove

      if [[ ${cutoff_type} == "min" ]]; then
        awk -v awk_pc_idx="${awk_pc_idx}" -v min_cutoff="${cutoff}" '{ if ( $awk_pc_idx < min_cutoff && $1 != "FID" && $2 != "IID" ) print $1, $2 }' ${pca_fname} \
          | grep -v "mis_pop"  >> ${tmp_fname}
      elif [[ ${cutoff_type} == "max" ]]; then
        awk -v awk_pc_idx="${awk_pc_idx}" -v max_cutoff="${cutoff}" '{ if ( $awk_pc_idx > max_cutoff && $1 != "FID" && $2 != "IID" ) print $1, $2 }' ${pca_fname} \
          | grep -v "mis_pop"  >> ${tmp_fname}
      fi

      local wc_new=$( cat ${tmp_fname} | wc -l )
      local n_to_remove=$(( ${wc_new} - ${wc_old} ))

      if [[ ${cutoff_type} == "min" ]]; then
        echo -e "\t${n_to_remove} samples < ${cutoff}"
      elif [[ ${cutoff_type} == "max" ]]; then
        echo -e "\t${n_to_remove} samples > ${cutoff}"
      fi
    fi
  else
    echo -e "\tNo ${cutoff_type} cutoff"
  fi
}

filter_by_pc () {
  # filters by PCs based on min/max cutoff arrays

  # load arguments
  local pca_fname=$1 # file name of PCA file (could be the file with or without ref samples)
  local with_ref=$2 # whether PCA version was done with or without reference samples
  declare -a min_cutoffs=("${!3}") # array of cutoffs (the ! in ${...} expands the argument)
  declare -a max_cutoffs=("${!4}") # array of cutoffs

  declare -a mincutoffs=(${min_cutoffs[@]})
  declare -a maxcutoffs=(${max_cutoffs[@]})

  # initialize variables dependent on function arguments
  local pca_version="PCA w/$( if ! ${with_ref}; then echo "o" ; fi ) ref" # used for later print statements
  local -i min_cutoffs_len=${#mincutoffs[@]} # length of array
  local -i max_cutoffs_len=${#maxcutoffs[@]} #length of array
  local -i max_pc=$( get_max ${min_cutoffs_len} ${max_cutoffs_len} ) # get maximum PC from min/mix cutoffs, required for iteration

  if ${with_ref}; then
    local remove_fname=${OUT}.remove.withref.txt
  else
    local remove_fname=${OUT}.remove.withoutref.txt
  fi

  local tmp_fname=tmp.${remove_fname}
  touch ${tmp_fname} # create empty file so that wc_old can be calculated even if file is empty

  if [[ ${max_pc} -gt 0 ]]; then
    echo -e "\nFiltering on ${pca_version}"
    for pc in `seq 1 ${max_pc}`; do
      echo "PC${pc}:"
      awk_filter ${pca_fname} ${tmp_fname} ${pc} "min" mincutoffs[@]
      awk_filter ${pca_fname} ${tmp_fname} ${pc} "max" maxcutoffs[@]
    done
  else
    echo "No PC filtering specified using ${pca_version}"
  fi

  sort  ${tmp_fname} | uniq > ${remove_fname}
  rm ${tmp_fname}
  echo "Total samples to remove from ${pca_version}: $(cat ${remove_fname} | wc -l)"
}

# run PCA filter
filter_by_pc ${PCA} false mincutoffs[@] maxcutoffs[@] # call filter_by_pc function for PCA w/o ref
filter_by_pc ${PCA_W_REF} true mincutoffswref[@] maxcutoffswref[@] # call filter_by_pc function for PCA w/ ref

echo "Total samples to remove from PCA: $( cat ${OUT}.remove.with{,out}ref.txt | sort | uniq | wc -l )"

echo -e "\nFiltering samples from ${FAMEX}"
famex_to_remove=$( cat ${FAMEX} | wc -l )
# NOTE: Assumes that FIDs have Ricopili "tags" with case/control info
if [[ ${famex_to_remove} -gt 0 ]]; then
  echo -e "\tCases: $( cut -d '_' -f1 ${FAMEX} | grep "cas" | wc -l)"
  echo -e "\tControls: $( cut -d '_' -f1 ${FAMEX} | grep "con" | wc -l)"
fi
echo "Total samples to remove using ${FAMEX}: ${famex_to_remove}"

# join all ${STUDY}.remove* files
# NOTE: It is possible for there to be samples in the famex file that also appear
#	in the *remove.*.txt files because the famex file comes from the PCA
#	w/o reference. The relatedness filter is not consistent between PCA w/o
#	reference vs. PCA w/ reference.
if ${REMOVE_PREFIX}; then
  # OPTION 1: removes QC prefix from FIDs
  # (assumes that the only '*' in each line occurs between the preimp_qc prefix and the original FID)
  cat ${OUT}.remove.with{,out}ref.txt <(cut -f 1-2 $STUDY.mepr.famex) \
    | column -t \
    | sort \
    | uniq \
    | cut -d'*' -f 2 > ${OUT}.remove.txt
else
  # OPTION 2: Does not remove QC prefix from FIDs
  cat ${OUT}.remove.with{,out}ref.txt <(cut -f 1-2 $STUDY.mepr.famex) \
    | column -t \
    | sort \
    | uniq  > ${OUT}.remove.txt
fi

echo -e "\nTotal samples to remove: $(cat ${OUT}.remove.txt | wc -l)"

if ${PLINK_REMOVE}; then

  files=($BFILE.{bed,bim,fam})
  for file in ${files[@]}; do
    test ! -f $file && echo "ERROR: File $file not found." && exit 1
  done

  echo # add carriage return before PLINK call

  # remove individuals outside of PC bounds
  plink --bfile ${BFILE} --remove ${OUT}.remove.txt --make-bed --out ${OUT}

  # check if any samples were removed
  ct_old=$(cat ${OUT}.fam | wc -l)
  ct_new=$(cat ${BFILE}.fam | wc -l)
  if [[ ${ct_old} -eq ${ct_new} ]]; then
    echo "WARNING: No samples were removed. All output files will be deleted."
    rm ${OUT}.remove{,.with{,out}ref}.txt  ${OUT}.{bed,bim,fam,log,hh}
  fi
fi
