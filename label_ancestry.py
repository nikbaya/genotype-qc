# plink --bfile aut_spark_mix_nb-qc1 --bmerge pop_4pop_mix_SEQ --make-bed --out merged
# plink --bfile merged --geno 1e-6 --make-bed --out merged_v2

# conda init bash
# conda create --name py3 python=3.6.7
# conda activate py3
# conda install -yf pandas


import pandas as pd

df = pd.read_csv('merged_v2.fam',delim_whitespace=True,header=None)

# df = df.rename(columns={6:'pop'})

df['pop'] = df[0].str.split('_',expand=True)[3]
df.loc[df['pop'] == 'mix', 'pop'] = float('nan')

df['pop'].to_csv('SPARK.parents.IMUS.menv.pop',sep=' ',header=False,index=False,na_rep='-')


awk -F' ' '{ print $7 }' merged_v2.fam > merged_v2.pop