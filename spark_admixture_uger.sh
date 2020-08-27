#!/bin/bash
#$ -l h_rt=0:00:01:00
#$ -j y
#$ -l h_vmem=10g
#$ -l nodes=1
#$ -cwd
#$ -o /stanley/genetics/users/nbaya/spark/array_May2019/spark_ricopili/preimputation/preimp7_imus/admixture_hgdp/UGERtest.log
#$ -t 1-5

source /broad/software/scripts/useuse
reuse -q R-3.4
reuse -q Anaconda3
reuse -q .google-cloud-sdk

myArray=( 4 5 6 7 8 )
k=${myArray[${SGE_TASK_ID}]}

run_admixture() {
	cd /stanley/genetics/users/nbaya/spark/array_May2019/spark_ricopili/preimputation/preimp7_imus/admixture_hgdp
	
	echo $k

	./admixture_linux-1.3.0/admixture
	
	# ./admixture_linux-1.3.0/admixture preimp7.founders.imus.hgdp_v3.menv.bed $k
	#
	# cp preimp7.founders.imus.hgdp_v3.menv.$k.P spark.hgdp.admixture_tmp2.P.in
	#
	# ./admixture_linux-1.3.0/admixture -P spark.hgdp.admixture_tmp2.bed $k

}


run_admixture ${k}


/psych/genetics_data/ripke/references_outdated/hapmap_ref/impute2_ref/1KG_Aug12/ALL_1000G_phase1integrated_v3_impute_macGT1/4pops/qc/pop_4pop_mix_SEQ
