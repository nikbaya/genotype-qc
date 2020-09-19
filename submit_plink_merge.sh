#!/usr/bin/env bash
#
# wrapper script for job submission on Broad UGER cluster
#
# The -V below above will provoke a warning that
# LD_LIBRARY_PATH won't be used for security reasons;
# this warning can be safely ignored
#
# Example: Submit plink_merge.sh to run on the Broad cluster and
# merge two datasets with bfile prefixes "data" and "ref"
#
# qsub submit_plink_merge.sh "data" "ref"
#
# NOTE: See notes in plink_merge.sh
#
# Author: Nikolas Baya (2020-09-19)

#$ -j y
#$ -cwd
#$ -V
#$ -N plink_merge
#$ -o ./plink_merge.log
#$ -q broad
#$ -l h_rt=1:00:00 # runtime in hours (currently set to one hour)
#$ -l m_mem_free=10g,h_vmem=10g # memory requested (currently set to 10 GB)

./plink_merge.sh $1 $2
