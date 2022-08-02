import os
import pandas as pd
os.chdir('matrix/')
df1=pd.read_csv('bhmatrix.txt',sep='\t',index_col='transcript')
df5=pd.read_csv('bhfpkm.txt',sep='\t',index_col='transcript')

df3=pd.DataFrame()
for i in df1.columns:
    df5[i]=df5[i]*df1[i].sum()
    df3[i]=(df5[i]/df5[i].sum())*1e6

df3.to_csv('bh_tpm2.txt',sep='\t')