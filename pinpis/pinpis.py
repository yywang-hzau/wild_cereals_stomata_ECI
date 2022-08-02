import os
import pandas as pd

def get_pinpis(dir,file):
    os.chdir(dir+file)
    df1=pd.read_csv('product_results.txt',sep='\t')
    #df1['pinpis']=df1['piN'].map(float)/df1['piS'].map(float)
    if 'Bh' in file:
        df1.to_csv('../bhpinpis.txt',sep='\t',index=None,mode='a')
    else:
        df1.to_csv('../bspinpis.txt',sep='\t',index=None,mode='a')
    os.chdir(dir)

def group_by_gene(dir):
    os.chdir(dir)
    df1=pd.read_csv('n_s_final.txt',sep='\t')
    df2=df1.groupby(df1['product']).mean()

    df2.to_csv('group_by_gene.txt',sep='\t')

def merge(dir):
    os.chdir(dir)
    df1 = pd.read_csv('degs_all.txt', sep='\t')
    df2=pd.read_csv('pinpis_morethan1.txt', sep='\t')
    df3=pd.merge(df1,df2,how='right',on='gene_id')
    df3.to_csv('piNpiS_plot.txt', sep='\t')

def main():
    dir='/Users/yuanyuan/Downloads/zhc/pinpis/'
    #for file in os.listdir(dir):
        #if file.endswith('snpgenie'):
            #print(file)
            #get_pinpis(dir,file)

    group_by_gene(dir)
    #merge(dir)
main()