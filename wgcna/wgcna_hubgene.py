
import os
import pandas as pd

def filt_mm():
    f1=open('datkme.txt','r')
    list1=f1.readlines()
    f2=open('mm_filt.txt','a+')
    f2.write('GeneName'+'\t'+'mm'+'\t'+'block'+'\n')
    for i in list1[1:]:
        i=i.strip()
        len1=len(i.split('\t'))
        #灰色不要
        len1=len1-1
        for n in range(1,len1):
            float1=i.split('\t')[n]
            if abs(float(float1)) > 0.8:
                label=list1[0].split('\t')[n-1]
                id=(i.split('\t')[0])
                f2.write(id + '\t'+float1+'\t'+label+'\n')
    f1.close()
    f2.close()

def filt_gs(file):
    df1=pd.read_csv(file,sep=',',index_col=0)
    df2=df1.loc[(abs(df1['PearsonCorrelation'])>= 0.6) &  (df1['q.Standard']<= 0.05)]
    df2.to_csv('gs.txt',sep='\t',mode='a+',index=0)

def duplicate_merge():
    df1=pd.read_csv('mm_filt.txt',sep='\t')
    df2=pd.read_csv('gs.txt',sep='\t')
    df3=df1.drop_duplicates(subset='GeneName',keep='first', inplace=False)
    df4=df2.drop_duplicates(subset='GeneName',keep='first', inplace=False)
    df5=pd.merge(df4,df3,on='GeneName',how='inner')
    #df4['sort'] = df4.mm.abs()
    #df4=df4.sort_values(by='sort',ascending=True)
    #df4=df4.drop_duplicates(subset='sort',keep='first', inplace=False)
    df5.to_csv('finalhub005.txt',sep='\t',mode='a+',index=0)

def main():
    dir='21es/'
    os.chdir(dir)
    #datakme (p value)
    filt_mm()
    for file in os.listdir(dir):
        if os.path.isfile(file) and file.endswith('bhexpr.csv'):
            filt_gs(file)
    duplicate_merge()


main()
print('finish')