# cd /stanley/genetics/users/nbaya/spark/array_May2019/spark_ricopili/preimputation/spark_preimp1
# python

import pandas as pd

genome = pd.read_csv('SPARK.27K.genotype.20190501.genome',delim_whitespace=True)
fam = pd.read_csv('SPARK.27K.genotype.20190501_preimp1.fam',delim_whitespace=True,header=None)
fam = fam.rename(columns={0:'FID',1:'IID',2:'PAT',3:'MAT'})

w_both_parents = fam[(fam.PAT!='0')&(fam.MAT!='0')]

p1 = w_both_parents.PAT.values.tolist()
p2 = w_both_parents.MAT.values.tolist()

related_subset = genome[(genome.IID1.isin(p1)&genome.IID2.isin(p2))|(genome.IID1.isin(p2)&genome.IID2.isin(p1))]

related_pairs = list(zip(related_subset.IID1.values.tolist(),related_subset.IID2.values.tolist()))

related_parents = []



for pair in related_pairs:
	if len(w_both_parents[(w_both_parents.PAT==pair[0])&(w_both_parents.MAT==pair[1])])>0 or len(w_both_parents[(w_both_parents.PAT==pair[1])&(w_both_parents.MAT==pair[0])])>0:
        related_parents.append(pair[0])
        related_parents.append(pair[1])
        
parents = fam[(fam.IID.isin(fam.PAT))|(fam.IID.isin(fam.MAT))]
unrelated_parents = parents[~parents.IID.isin(related_parents)]
related_parents = parents[parents.IID.isin(related_parents)]

related_parents.to_csv('SPARK.related_parents.fam',sep=' ',index=False)
unrelated_parents.to_csv('SPARK.unrelated_parents.fam',sep=' ',index=False)